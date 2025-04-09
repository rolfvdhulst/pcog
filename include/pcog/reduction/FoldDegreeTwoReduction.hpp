//
// Created by rolf on 10-7-23.
//

#ifndef PCOG_SRC_REDUCTION_FOLDDEGREETWOREDUCTION_HPP
#define PCOG_SRC_REDUCTION_FOLDDEGREETWOREDUCTION_HPP

#include <utility>

#include "DenseReductionGraph.hpp"
#include "ReductionVertexQueue.hpp"
#include "pcog/utilities/Coloring.hpp"

namespace pcog {
class ReductionStack;
///This reduction handles some nodes which have degree 2 in the complement graph,
/// e.g. which only have two non-neighbours in the current graph
struct FoldDegreeTwoReduction {
   FoldDegreeTwoReduction(Node node, Node kept, Node removed,
                          DenseSet t_keepNonNeighbours, DenseSet t_removeNonNeighbours)
       : degree2Node{node}, keptNode{kept}, removedNode{removed},
         keepNonNeighbours{std::move(t_keepNonNeighbours)}, removeNonNeighbours{std::move(t_removeNonNeighbours)}
         {};
   Node degree2Node;
   Node keptNode;
   Node removedNode;

   DenseSet keepNonNeighbours;
   DenseSet removeNonNeighbours;

   void transformStableSet(DenseSet& set) const;
   void newToOldColoring(NodeColoring &coloring) const;
};

bool foldDegreeTwoReduceNode(Node node,
                         DenseReductionGraph& graph,
                         ReductionStack& stack,
                         ReductionVertexQueue& queue
);
}

#endif // PCOG_SRC_REDUCTION_FOLDDEGREETWOREDUCTION_HPP
