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

   bool setIsStable = true;
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
   assert(setIsStable);
   SimplicialReduction reduction(node,set);
   stack.push(reduction);

   //As a by-result, disallowed_nodes already exactly contains all of the affected neighbours,
   for(const auto& neighbour : disallowed_nodes){
      assert(graph.containsNode(neighbour));
      queue.push(neighbour);
   }

   for(Node setNode : set){
      graph.removeNode(setNode);
   }
   return true;
}
SimplicialReduction::SimplicialReduction(Node node, DenseSet set) : node{node}, set{set}{}

}