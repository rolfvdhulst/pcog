//
// Created by rolf on 6-7-23.
//

#ifndef PCOG_INCLUDE_PCOG_REDUCTION_REDUCER_HPP
#define PCOG_INCLUDE_PCOG_REDUCTION_REDUCER_HPP

#include "ReductionStack.hpp"
#include "ReductionVertexQueue.hpp"
#include "DenseReductionGraph.hpp"
namespace pcog{

void reduceGraph(const DenseGraph& graph);
}

#endif // PCOG_INCLUDE_PCOG_REDUCTION_REDUCER_HPP
