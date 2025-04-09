//
// Created by rolf on 4-7-23.
//

#ifndef PCOG_INCLUDE_PCOG_REDUCTION_LOWDEGREEREDUCTION_HPP
#define PCOG_INCLUDE_PCOG_REDUCTION_LOWDEGREEREDUCTION_HPP
#include "DenseReductionGraph.hpp"
#include "ReductionVertexQueue.hpp"
#include "pcog/utilities/Coloring.hpp"

namespace pcog {
class ReductionStack;
struct LowDegreeReduction {
   LowDegreeReduction(Node t_node, DenseSet t_neighbours) : node{t_node}, neighbours{std::move(t_neighbours)}{};
   Node node;
   DenseSet neighbours;

   void transformStableSet(DenseSet& set) const;
   void newToOldColoring(NodeColoring &coloring) const;
};

bool lowDegreeReduceNode(Node node,
                   DenseReductionGraph& graph,
                   ReductionStack& stack,
                   ReductionVertexQueue& queue
                   );
}
#endif // PCOG_INCLUDE_PCOG_REDUCTION_LOWDEGREEREDUCTION_HPP
