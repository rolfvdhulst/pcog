//
// Created by rolf on 18-7-23.
//

#ifndef PCOG_TWOFIXINGREDUCTION_HPP
#define PCOG_TWOFIXINGREDUCTION_HPP

#include "DenseReductionGraph.hpp"
#include "ReductionVertexQueue.hpp"
#include "pcog/utilities/Coloring.hpp"

namespace pcog {
class ReductionStack;
struct TwoFixingReduction {
   Node firstNode;
   Node secondNode;
   DenseSet firstSet;
   DenseSet secondSet;
   void transformStableSet(DenseSet& set) const;
   void newToOldColoring(NodeColoring &coloring) const;
};

bool twoFixingReduceNode(Node node,
                          DenseReductionGraph& graph,
                          ReductionStack& stack,
                          ReductionVertexQueue& queue
);

bool attemptApplyTwoFixing(Node node, Node other,
                           DenseReductionGraph& graph,
                           ReductionStack& stack,
                           ReductionVertexQueue& queue);

bool findTwoFixings(DenseReductionGraph& graph,
                    ReductionStack& stack,
                    ReductionVertexQueue& queue);
}

#endif // PCOG_TWOFIXINGREDUCTION_HPP
