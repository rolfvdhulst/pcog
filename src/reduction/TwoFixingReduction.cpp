//
// Created by rolf on 18-7-23.
//

#include "pcog/reduction/TwoFixingReduction.hpp"
#include "pcog/reduction/ReductionStack.hpp"

namespace pcog {
std::vector<DenseSet> attemptTwoFixing(Node first, Node second,
                                       const DenseReductionGraph& graph){
   DenseSet firstNeighbours = graph.complementNeighbourhood(first);
   DenseSet secondNeighbours = graph.complementNeighbourhood(second);
   DenseSet firstMinSecond = firstNeighbours.difference(secondNeighbours);
   if(firstMinSecond.empty()){
      return {};//Other reductions handle this case
   }
   for(Node node : firstMinSecond){
      if(graph.neighbourhood(node).intersects(firstMinSecond)){
         return {};
      }
   }
   DenseSet secondMinFirst = secondNeighbours.difference(firstNeighbours);
   if(secondMinFirst.empty()){
      return {};//Other reductions handle this case
   }
   for(Node node : secondMinFirst){
      if(graph.neighbourhood(node).intersects(secondMinFirst)){
         return {};
      }
   }
   firstMinSecond.add(first);
   secondMinFirst.add(second);


   //Use the DSATUR heuristic to divide the remaining nodes in the intersection
   //Since DSATUR is exact for bipartite graphs, this always finds a bipartition (if one exists)
   //Or shows that the subgraph formed by the non-neighbours and the two given nodes is not bipartite
   DenseSet firstSecondIntersection = firstNeighbours.intersection(secondNeighbours);
   if(firstSecondIntersection.empty()){
      return {}; //Case is handled by simplicial reduction
   }
   //Set up the statistics necessary for dsatur. We create a mapping from main graph to subgraph first
   struct NodeData{
      bool hasFirstColorNeighbour = false;
      bool hasSecondColorNeighbour = false;
      int saturationDegree = 0;
      std::size_t subgraphDegree = 0;
   };
   std::vector<NodeData> data;
   std::vector<Node> newToOld;
   for(Node node : firstSecondIntersection){
      bool firstNeighbour = graph.neighbourhood(node).intersects(firstMinSecond);
      bool secondNeighbour = graph.neighbourhood(node).intersects(secondMinFirst);
      if(firstNeighbour && secondNeighbour){
         //clique size 3 found, we terminate
         return {};
      }
      std::size_t subgraphDegreeNode = graph.neighbourhood(node).intersection(firstSecondIntersection).size();
      data.push_back(NodeData{
          .hasFirstColorNeighbour = firstNeighbour,
          .hasSecondColorNeighbour = secondNeighbour,
          .saturationDegree = (firstNeighbour ? 1 : 0) + (secondNeighbour ? 1 : 0),
          .subgraphDegree = subgraphDegreeNode,
      });
      newToOld.push_back(node);
   }

   std::vector<Node> oldToNew(firstSecondIntersection.capacity(),INVALID_NODE);
   for(std::size_t index = 0; index < newToOld.size(); ++index){
      oldToNew[newToOld[index]] = index;
   }

   std::size_t numQueueNodes = firstSecondIntersection.size();
   while(numQueueNodes != 0){
      //Find node with highest saturation degree
      Node pickNode = INVALID_NODE;
      int pickSaturation = -1;
      std::size_t pickDegree = 0;
      for(Node oldNode : firstSecondIntersection){
         Node newNode = oldToNew[oldNode];
         if(data[newNode].saturationDegree > pickSaturation ||
             (data[newNode].saturationDegree == pickSaturation && data[newNode].subgraphDegree > pickDegree)){
            pickNode = newNode;
            pickSaturation = data[newNode].saturationDegree;
            pickDegree = data[newNode].subgraphDegree;
         }
      }
      numQueueNodes -= 1;
      //Color it

      if(pickSaturation >= 2){
         assert(data[pickNode].hasFirstColorNeighbour && data[pickNode].hasSecondColorNeighbour);
         return {};// we need 3 or more colors
      }
      bool assignedToFirstColor = !data[pickNode].hasFirstColorNeighbour;
      if(assignedToFirstColor){
         firstMinSecond.add(newToOld[pickNode]);
      }else{
         secondMinFirst.add(newToOld[pickNode]);
      }
      for(Node neighbour : graph.neighbourhood(newToOld[pickNode]).intersection(firstSecondIntersection)){
         Node newNeighbour = oldToNew[neighbour];
         if(assignedToFirstColor){
            if(!data[newNeighbour].hasFirstColorNeighbour){
               data[newNeighbour].saturationDegree += 1;
            }
            data[newNeighbour].hasFirstColorNeighbour = true;
         }else{
            if(!data[newNeighbour].hasSecondColorNeighbour){
               data[newNeighbour].saturationDegree += 1;
            }
            data[newNeighbour].hasSecondColorNeighbour = true;
         }
      }
      firstSecondIntersection.remove(newToOld[pickNode]);
   }
   return {firstMinSecond,secondMinFirst};
}

bool twoFixingReduceNode(Node node, DenseReductionGraph &graph,
                         ReductionStack &stack, ReductionVertexQueue &queue) {
   for(Node other : graph.neighbourhood(node)){
      auto fixedSets = attemptTwoFixing(node,other,graph);
      if(!fixedSets.empty()){
         //TODO: double check that found sets are indeed stable in the graph
         TwoFixingReduction reduction{
            .firstNode = node,
            .secondNode = other,
            .firstSet = fixedSets.front(),
            .secondSet = fixedSets.back(),
         };
         stack.push(reduction);
         graph.removeStableSet(reduction.firstSet);
         graph.removeStableSet(reduction.secondSet);
         //TODO: think about what nodes to push
         for(Node graphNode : graph.nodes()){
            queue.push(graphNode,graph.lowerBoundNodes().contains(graphNode));
         }
         return true;
      }
   }
   return false;
}
}