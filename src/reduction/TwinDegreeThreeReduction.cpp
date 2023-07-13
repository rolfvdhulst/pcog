//
// Created by rolf on 12-7-23.
//

#include "pcog/reduction/TwinDegreeThreeReduction.hpp"
#include "pcog/reduction/ReductionStack.hpp"

namespace pcog {
bool twinDegreeThreeReduction(Node node, DenseReductionGraph &graph,
                              ReductionStack &stack,
                              ReductionVertexQueue &queue) {
   if(graph.nodeDegree(node) != (graph.numNodes() - 4)){
      return false;
   }
   DenseSet nonNeighbours = graph.complementNeighbourhood(node);

   Node w = nonNeighbours.first();
   Node x = nonNeighbours.find_next(w);
   Node y = nonNeighbours.find_next(x);
   assert(w != INVALID_NODE && x != INVALID_NODE && y != INVALID_NODE);
   assert(nonNeighbours.find_next(y) == INVALID_NODE);

   DenseSet twinCandidates = graph.neighbourhood(node).difference(graph.neighbourhood(w));
   twinCandidates.inplaceDifference(graph.neighbourhood(x));
   twinCandidates.inplaceDifference(graph.neighbourhood(y));

   Node twin = INVALID_NODE;
   for(Node candidate : twinCandidates){
      if(graph.nodeDegree(candidate) == graph.numNodes()-4){
         twin = candidate;
         break;
      }
   }
   if(twin == INVALID_NODE){
      return false;
   }
   bool edgeWX = graph.neighbourhood(w).contains(x);
   bool edgeWY = graph.neighbourhood(w).contains(y);
   bool edgeXY = graph.neighbourhood(x).contains(y);
   if(!(edgeWX && edgeWY && edgeXY)){
      //There is at least one non-edge, we can remove node twin and wxy from the graph by fixing two stable sets
      DenseSet firstSet(graph.nodes().capacity());
      firstSet.add(node);
      firstSet.add(w);
      if(!edgeWX){
         firstSet.add(x);
      }
      if(!edgeWY){
         firstSet.add(y);
      }
      DenseSet secondSet(graph.nodes().capacity());
      secondSet.add(twin);
      if(!edgeWX){
         secondSet.add(y);
      }else if(!edgeWY){
         secondSet.add(x);
      }else{
         assert(!edgeXY);
         secondSet.add(x);
         secondSet.add(y);
      }

      assert(firstSet.setUnion(secondSet).size() >= 5);
      graph.removeStableSet(firstSet);
      graph.removeStableSet(secondSet);

      stack.push(TwinDegreeThreeReduction{.u= node,.v = twin,.w = w, .x = x,.y = y,.fold = false});
      //TODO: what nodes to push to the queue? (perhaps intersection of neighbourhoods? or non-neighbourhoods)
      for(const auto& addNode : graph.nodes()){
         queue.push(addNode,graph.lowerBoundNodes().contains(addNode));
      }

      return true;
   }
   //Perform a folding reduction. Similar to degree 2 reduction, pick which nodes to remove to preserve clique lb;
   //If u and/or v are in clique, then no problems; folding does not destroy the clique. We arbitrarily pick from wxy
   //If one or two nodes of wxy are in the clique, then remove these clique nodes first (otherwise, clique lb might be invalidated)
   //If wxy are all in the clique, we can arbitrarily remove two of them without any problems.

   Node keepNode,removeNode1,removeNode2;
   if(!graph.lowerBoundNodes().contains(x)){
      //Fold neighbourhoods, keeping x
      keepNode = x;
      removeNode1 = w;
      removeNode2 = y;
   }else if(!graph.lowerBoundNodes().contains(y)) {
      //Fold neighbourhoods keeping y
      keepNode = y;
      removeNode1 = w;
      removeNode2 = x;
   }else{
      //keep w and fold neighbourhoods
      keepNode = w;
      removeNode1 = x;
      removeNode2 = y;
   }
   DenseSet keepRemoveNeighbours = graph.neighbourhood(keepNode);
   DenseSet possibleRemoveEdges = graph.complementNeighbourhood(removeNode1).setUnion(graph.complementNeighbourhood(removeNode2));
   keepRemoveNeighbours.inplaceIntersection(possibleRemoveEdges);
   graph.removeNode(node);
   graph.removeNode(twin);
   graph.removeNode(removeNode1);
   graph.removeNode(removeNode2);
   graph.removeNodeEdges(keepNode,keepRemoveNeighbours);

   stack.push(TwinDegreeThreeReduction{.u= node,.v = twin,.w = w, .x = x,.y = y,.fold = true});
   //TODO: what nodes to push to the queue? (perhaps intersection of neighbourhoods? or non-neighbourhoods)
   for(const auto& addNode : graph.nodes()){
      queue.push(addNode,graph.lowerBoundNodes().contains(addNode));
   }

   return true;
}
}