//
// Created by rolf on 4-7-23.
//

#ifndef PCOG_INCLUDE_PCOG_REDUCTION_LOWDEGREEREDUCTION_HPP
#define PCOG_INCLUDE_PCOG_REDUCTION_LOWDEGREEREDUCTION_HPP
#include "DenseReductionGraph.hpp"
#include "ReductionVertexQueue.hpp"

namespace pcog {
class ReductionStack;
struct LowDegreeReduction {
   LowDegreeReduction(Node node, DenseSet neighbours) : node{node},neighbours{std::move(neighbours)}{};
   Node node;
   DenseSet neighbours;
};

void lowDegreeReduceNode(Node node,
                   DenseReductionGraph& graph,
                   ReductionStack& stack,
                   ReductionVertexQueue& queue
                   );
}
#endif // PCOG_INCLUDE_PCOG_REDUCTION_LOWDEGREEREDUCTION_HPP
