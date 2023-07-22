//
// Created by rolf on 12-7-23.
//

#ifndef PCOG_TWINDEGREETHREEREDUCTION_HPP
#define PCOG_TWINDEGREETHREEREDUCTION_HPP

#include "DenseReductionGraph.hpp"
#include "ReductionVertexQueue.hpp"
#include "pcog/utilities/Coloring.hpp"
namespace pcog{
class ReductionStack;
///This reduction handles a pair of nodes which have degree 3 in the complement graph,
struct TwinDegreeThreeFoldReduction {
   //degree 3 nodes
   Node u;
   Node v;

   //Non-neighbours
   Node keep;
   Node remove1;
   Node remove2;
   DenseSet keepNonNeighbours;
   DenseSet remove1NonNeighbours;
   DenseSet remove2NonNeighbours;

   void transformStableSet(DenseSet& set) const;
   void newToOldColoring(NodeColoring &coloring) const;
};

struct TwinDegreeThreeReduction {

   //degree 3 nodes
   Node u;
   Node v;

   DenseSet uSet;
   DenseSet vSet;
   //Non-neighbours

   void transformStableSet(DenseSet& set) const;
   void newToOldColoring(NodeColoring &coloring) const;
};

bool twinDegreeThreeReduction(Node node,
                             DenseReductionGraph& graph,
                             ReductionStack& stack,
                             ReductionVertexQueue& queue
);
}

#endif // PCOG_TWINDEGREETHREEREDUCTION_HPP
