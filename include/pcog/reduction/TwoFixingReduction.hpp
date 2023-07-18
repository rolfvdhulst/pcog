//
// Created by rolf on 18-7-23.
//

#ifndef PCOG_TWOFIXINGREDUCTION_HPP
#define PCOG_TWOFIXINGREDUCTION_HPP

#include "DenseReductionGraph.hpp"
#include "ReductionVertexQueue.hpp"

namespace pcog {
class ReductionStack;
struct TwoFixingReduction {
   Node firstNode;
   Node secondNode;
   DenseSet firstSet;
   DenseSet secondSet;
};

bool twoFixingReduceNode(Node node,
                          DenseReductionGraph& graph,
                          ReductionStack& stack,
                          ReductionVertexQueue& queue
);
}

#endif // PCOG_TWOFIXINGREDUCTION_HPP
