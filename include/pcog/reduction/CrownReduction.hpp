//
// Created by rolf on 15-7-23.
//

#ifndef PCOG_CROWNREDUCTION_HPP
#define PCOG_CROWNREDUCTION_HPP

#include "DenseReductionGraph.hpp"
#include "ReductionVertexQueue.hpp"
#include "pcog/utilities/Coloring.hpp"

namespace pcog {
class ReductionStack;
///This reduction handles some nodes which have degree 2 in the complement graph,
/// e.g. which only have two non-neighbours in the current graph
struct CrownReduction {
   std::vector<DenseSet> fixedSets;

   void transformStableSet(DenseSet& set) const;
   void newToOldColoring(NodeColoring &coloring) const;
};

bool findCrownReductions(DenseReductionGraph& graph,
                         ReductionStack& stack,
                         ReductionVertexQueue& queue
);
}

#endif // PCOG_CROWNREDUCTION_HPP
