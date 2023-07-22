//
// Created by rolf on 4-7-23.
//

#include "pcog/reduction/LowDegreeReduction.hpp"
#include "pcog/reduction/ReductionStack.hpp"

namespace pcog{
bool lowDegreeReduceNode(Node node,
                         DenseReductionGraph& graph,
                         ReductionStack& stack,
                         ReductionVertexQueue& queue
){
   if(!graph.hasLowerBound()){
      return false;
   }
   std::size_t degree = graph.nodeDegree(node);
   std::size_t lb = graph.lowerBound();

   bool inClique = graph.lowerBoundNodes().contains(node);
   if(inClique){
      return false;
   }
   if(degree >= lb){ //TODO: what to do with clique nodes with 'small' degree?
                     // For now, ignore but think about this case when considering other subgraphs as lower bounds
      return false;
   }

   assert((!graph.lowerBoundNodes().contains(node) && degree < lb) || degree+1 < lb);
   LowDegreeReduction reduction(node,graph.neighbourhood(node));
   stack.push(reduction);
   for(const auto& neighbour : graph.neighbourhood(node)){
      queue.push(neighbour,graph.lowerBoundNodes().contains(node));
   }
   graph.removeNode(node);

   return true;
}

void LowDegreeReduction::transformStableSet(DenseSet &set) const {
   assert(!set.contains(node));
   if(!set.intersects(neighbours)){
      set.add(node);
   }
}
void LowDegreeReduction::newToOldColoring(NodeColoring &coloring) const {
   DenseSet unusedColors(coloring.numColors(),true);
   for(Node neighbour : neighbours){
      assert(coloring[neighbour] != INVALID_NODE);
      unusedColors.remove(coloring[neighbour]);
   }
   Color nodeColor = unusedColors.first();
   assert(nodeColor != INVALID_NODE);
   coloring[node] = nodeColor;

}
}