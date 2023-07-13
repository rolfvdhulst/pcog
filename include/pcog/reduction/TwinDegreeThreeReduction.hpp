//
// Created by rolf on 12-7-23.
//

#ifndef PCOG_TWINDEGREETHREEREDUCTION_HPP
#define PCOG_TWINDEGREETHREEREDUCTION_HPP

#include "DenseReductionGraph.hpp"
#include "ReductionVertexQueue.hpp"
namespace pcog{
class ReductionStack;
///This reduction handles a pair of nodes which have degree 3 in the complement graph,
struct TwinDegreeThreeReduction {
   //degree 3 nodes
   Node u;
   Node v;

   //Non-neighbours
   Node w;
   Node x;
   Node y;
   bool fold;
};

bool twinDegreeThreeReduction(Node node,
                             DenseReductionGraph& graph,
                             ReductionStack& stack,
                             ReductionVertexQueue& queue
);
}

#endif // PCOG_TWINDEGREETHREEREDUCTION_HPP
