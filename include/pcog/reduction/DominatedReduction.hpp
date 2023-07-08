//
// Created by rolf on 4-7-23.
//

#ifndef PCOG_INCLUDE_PCOG_REDUCTION_DOMINATEDREDUCTION_HPP
#define PCOG_INCLUDE_PCOG_REDUCTION_DOMINATEDREDUCTION_HPP

#include "DenseReductionGraph.hpp"
#include "ReductionVertexQueue.hpp"

namespace pcog {
class ReductionStack;
struct DominatedReduction {
   DominatedReduction(Node node, Node dominatingNode);
   Node node;
   Node dominatedBy;
};

bool dominatedReduceNodeDense(Node node,
                         DenseReductionGraph& graph,
                         ReductionStack& stack,
                         ReductionVertexQueue& queue
);
bool dominatedReduceNode(Node node,
                         DenseReductionGraph& graph,
                         ReductionStack& stack,
                         ReductionVertexQueue& queue
);
}
#endif // PCOG_INCLUDE_PCOG_REDUCTION_DOMINATEDREDUCTION_HPP
