//
// Created by rolf on 5-7-23.
//

#include "pcog/reduction/SimplicialReduction.hpp"
#include "pcog/reduction/ReductionStack.hpp"
namespace pcog{

bool simplicialReduceNode(Node node,
                          DenseReductionGraph& graph,
                          ReductionStack& stack,
                          ReductionVertexQueue& queue
){
   DenseSet set = graph.neighbourhood(node);
   set.complement();
   set.inplaceIntersection(graph.nodes());

   DenseSet disallowed_nodes(set.capacity());
   for(const auto& setNode : set){
      if(disallowed_nodes.contains(setNode)){
         return false;
      }
      disallowed_nodes.add(setNode);
      //Technically we would need to intersect the neighbourhood in the next line with present_nodes
      //however, we only check for nodes within the present set of nodes, so this is fine anyways
      disallowed_nodes.inplaceUnion(graph.neighbourhood(setNode));
   }
   SimplicialReduction reduction(node,set);
   stack.push(reduction);

   //As a by-result, disallowed_nodes already exactly contains all of the affected neighbours,
   for(const auto& neighbour : disallowed_nodes){
      assert(graph.containsNode(neighbour));
      queue.push(neighbour,graph.lowerBoundNodes().contains(neighbour));
   }

   graph.removeStableSet(set);
   return true;
}
SimplicialReduction::SimplicialReduction(Node t_node, DenseSet t_set) : node{t_node}, set{t_set}{}
void SimplicialReduction::transformStableSet(DenseSet &) const {
   //no-op, removed
}
void SimplicialReduction::newToOldColoring(NodeColoring &coloring) const {
   std::size_t count = coloring.numColors();
   for (Node fixed : set) {
      assert(coloring[fixed] == INVALID_COLOR);
      coloring[fixed] = count;
   }
   ++count;

   coloring.setNumColors(count);
}

}