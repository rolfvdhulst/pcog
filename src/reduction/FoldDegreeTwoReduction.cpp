//
// Created by rolf on 10-7-23.
//

#include "pcog/reduction/FoldDegreeTwoReduction.hpp"
#include "pcog/reduction/ReductionStack.hpp"
namespace pcog {

//Checks if two vertex sets form a biClique, e.g. if each edge between them exists
bool formBiClique(const DenseReductionGraph& graph, const DenseSet& firstSet, const DenseSet& secondSet){
   for(Node node : firstSet){
      if(!secondSet.isSubsetOf(graph.neighbourhood(node))){
         return false;
      }
   }
   //Since we have an undirected graph, also all second nodes are connected
   //to all first nodes by the above loop
   return true;
}
bool foldDegreeTwoReduceNode(Node node, DenseReductionGraph &graph,
                             ReductionStack &stack,
                             ReductionVertexQueue &queue) {
   if(graph.nodeDegree(node) != (graph.numNodes() - 3)){
      return false;
   }
   //Get the two non-neighbouring nodes
   DenseSet nonNodes = graph.neighbourhood(node);
   nonNodes.complement().inplaceIntersection(graph.nodes());
   nonNodes.remove(node);
   Node first = nonNodes.first();
   Node second = nonNodes.find_next(first);
   assert(nonNodes.find_next(second) == INVALID_NODE);
   //We can only apply the reduction if these nodes are connected with an edge
   //(otherwise, this is a simplicial vertex)
   if(!graph.neighbourhood(first).contains(second)){
      return false;
   }
   DenseSet firstNonNeighbours = graph.neighbourhood(first);
   firstNonNeighbours.complement().inplaceIntersection(graph.nodes());
   firstNonNeighbours.remove(first);

   DenseSet secondNonNeighbours = graph.neighbourhood(second);
   secondNonNeighbours.complement().inplaceIntersection(graph.nodes());
   secondNonNeighbours.remove(second);

   DenseSet firstMinSecond = firstNonNeighbours.difference(secondNonNeighbours);
   DenseSet secondMinFirst = secondNonNeighbours.difference(firstNonNeighbours);
   if(!formBiClique(graph,firstMinSecond,secondMinFirst)){
      return false;
   }
   //We can validly do the reduction.
   //We do this by removing u,v,w and replacing them by a single vertex that is
   //connected to all nodes except firstNonNeighbours and secondNonNeighbours

   //We pick the vertices in a way in which the computed clique lower bound remains a clique,
   //even after removing edges.
   //Concretely, we explicitly remove the original node and one other node.
   //If none of the two nodes is in the clique, v should be in the clique
   //If exactly one of the two nodes is in the clique, then we remove this node to ensure that the remaining graph's clique is valid
   //If both nodes are in the clique, none of their nonNeighbours are and thus we can safely remove these edges from one of these after removing the other.


   bool removeFirst = graph.lowerBoundNodes().contains(first) && !graph.lowerBoundNodes().contains(second);
   Node keepNode = removeFirst ?  second : first;
   Node removeNode = removeFirst ? first : second;
   const DenseSet& keepRemoveNeighbours = removeFirst ? firstMinSecond : secondMinFirst;
   const DenseSet& keepNonNeighbours = removeFirst ? secondNonNeighbours : firstNonNeighbours;
   const DenseSet& removeNonNeighbours = removeFirst ? firstNonNeighbours : secondNonNeighbours;
   stack.push(FoldDegreeTwoReduction(node,keepNode,removeNode,keepNonNeighbours,removeNonNeighbours));
   //TODO: what nodes to push to queue? Nodes in keepRemoveNeighbours get smaller degrees w.r.t lower bound
   queue.push(keepNode,graph.lowerBoundNodes().contains(keepNode));
   //TODO: think/test about if more/less nodes need to be pushed to the queue
   DenseSet queueNodes = firstNonNeighbours.setUnion(secondNonNeighbours);
   for(const auto& queueNode : queueNodes){
      if(queueNode == removeNode || queueNode == node ) continue;
      queue.push(queueNode,graph.lowerBoundNodes().contains(queueNode));
   }

   graph.removeNode(node);
   graph.removeNode(removeNode);
   graph.removeNodeEdges(keepNode,keepRemoveNeighbours);

   return true;
}
void FoldDegreeTwoReduction::transformStableSet(DenseSet &set) const {
   assert(!set.contains(degree2Node));
   assert(!set.contains(removedNode));
   if(set.contains(keptNode)){
      //If set is stable w.r.t removed neighbourhood of the removed node, we change it to be that
      if(set.isSubsetOf(removeNonNeighbours)){
         set.remove(keptNode);
         set.add(removedNode);
      }else{
         assert(set.isSubsetOf(keepNonNeighbours));
      }
   }
}
void FoldDegreeTwoReduction::newToOldColoring(NodeColoring &coloring) const {
   std::size_t count = coloring.numColors();

   assert(coloring[degree2Node] == INVALID_COLOR);
   assert(coloring[removedNode] == INVALID_COLOR);
   assert(coloring[keptNode] != INVALID_COLOR);

   Color keptNodeColor = coloring[keptNode];

   assert(keptNodeColor != INVALID_COLOR);

   bool inRemoveNonNeighbours = true;
   for(std::size_t i = 0; i < coloring.numNodes(); ++i){
      if(i == keptNode) continue;
      if(keptNodeColor == coloring[i] && !removeNonNeighbours.contains(i)){
         inRemoveNonNeighbours = false;
         break;
      }
   }

   coloring[degree2Node] = count;
   if(!inRemoveNonNeighbours){
      coloring[removedNode] = count;
   }else{
      coloring[removedNode] = coloring[keptNode];
      coloring[keptNode] = count;
   }

   ++count;
   coloring.setNumColors(count);
   //TODO?
}
}