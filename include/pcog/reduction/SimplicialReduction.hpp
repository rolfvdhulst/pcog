//
// Created by rolf on 5-7-23.
//

#ifndef PCOG_SRC_REDUCTION_SIMPLICIALREDUCTION_HPP
#define PCOG_SRC_REDUCTION_SIMPLICIALREDUCTION_HPP

#include "DenseReductionGraph.hpp"
#include "ReductionVertexQueue.hpp"
#include "pcog/utilities/Coloring.hpp"

namespace pcog {
class ReductionStack;
struct SimplicialReduction {
   SimplicialReduction(Node node, DenseSet set);
   Node node;
   DenseSet set;
   void transformStableSet(DenseSet& set) const;
   void newToOldColoring(NodeColoring &coloring) const;
};

bool simplicialReduceNode(Node node,
                         DenseReductionGraph& graph,
                         ReductionStack& stack,
                         ReductionVertexQueue& queue
);
}

#endif // PCOG_SRC_REDUCTION_SIMPLICIALREDUCTION_HPP
