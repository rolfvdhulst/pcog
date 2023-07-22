//
// Created by rolf on 4-7-23.
//

#include "pcog/reduction/DominatedReduction.hpp"
#include "pcog/reduction/ReductionStack.hpp"

namespace pcog{
DominatedReduction::DominatedReduction(Node node, Node dominatingNode) : node{node},dominatedBy{dominatingNode}{

}
void DominatedReduction::transformStableSet(DenseSet &set) const {
   assert(!set.contains(node));
   set.set(node,set.contains(dominatedBy));
}
void DominatedReduction::newToOldColoring(NodeColoring &coloring) const {
   assert(coloring[dominatedBy] != INVALID_COLOR);
   assert(coloring[node] == INVALID_COLOR);
   coloring[node] = coloring[dominatedBy];
}

//This function is faster for very dense graphs e.g. graphs where we are effectively
//solving the vertex clique cover on the sparse complement graph
bool dominatedReduceNodeDense(Node node,
                              DenseReductionGraph& graph,
                              ReductionStack& stack,
                              ReductionVertexQueue& queue){
   DenseSet nonAdjacent = graph.neighbourhood(node);
   nonAdjacent.complement();
   nonAdjacent.remove(node);
   for(Node candidateDominating : nonAdjacent){
      if(graph.nodeDegree(candidateDominating) < graph.nodeDegree(node)) continue;
      if(graph.neighbourhood(node).isSubsetOf(graph.neighbourhood(candidateDominating))){
         assert(candidateDominating != INVALID_NODE);
         assert(graph.containsNode(candidateDominating));
         assert(graph.neighbourhood(node).isSubsetOf(graph.neighbourhood(candidateDominating)));

         DominatedReduction reduction(node,candidateDominating);
         stack.push(reduction);
         for(const auto& neighbour : graph.neighbourhood(node)){
            queue.push(neighbour,graph.lowerBoundNodes().contains(neighbour));
         }
         graph.removeNode(node);
         return true;
      }
   }

   return false;
}
bool dominatedReduceNode(Node node,
                         DenseReductionGraph& graph,
                         ReductionStack& stack,
                         ReductionVertexQueue& queue
){
   const DenseSet& nodeNeighbourhood = graph.neighbourhood(node);
   Node iterNode = nodeNeighbourhood.first();
   if(iterNode == INVALID_NODE){
      //Node has degree 0; it will be removed by degree rule
      //but strictly speaking it is also 'dominated' by any other node
      return false;
   }
   //The node itself or any node in its neighbourhood cannot dominate it, so we exclude them already
   DenseSet dominatedCandidates = graph.neighbourhood(iterNode).difference(nodeNeighbourhood);
   dominatedCandidates.remove(node);
   iterNode = nodeNeighbourhood.find_next(iterNode);
   if(dominatedCandidates.empty()){
      return false;
   }
   while(iterNode != INVALID_NODE){
      dominatedCandidates.inplaceIntersection(graph.neighbourhood(iterNode));
      if(dominatedCandidates.empty()){
         return false;
      }
      iterNode = nodeNeighbourhood.find_next(iterNode);
   }
   Node dominatingNode = dominatedCandidates.first();
   assert(dominatingNode != INVALID_NODE);
   assert(graph.containsNode(dominatingNode));
   assert(graph.neighbourhood(node).isSubsetOf(graph.neighbourhood(dominatingNode)));

   DominatedReduction reduction(node,dominatingNode);
   stack.push(reduction);
   for(const auto& neighbour : nodeNeighbourhood){
      queue.push(neighbour,graph.lowerBoundNodes().contains(neighbour));
   }
   graph.removeNode(node);
   return true;
}

}
