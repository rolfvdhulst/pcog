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

//If a node is contained in a unique maximum stable set, we can fix this color
//and remove the stable set from the graph. For a given node v, we check this by
//checking if V \ N(v) is a stable set.
bool simplicialReduceNode(Node node,
                         DenseReductionGraph& graph,
                         ReductionStack& stack,
                         ReductionVertexQueue& queue
);
}

#endif // PCOG_SRC_REDUCTION_SIMPLICIALREDUCTION_HPP
