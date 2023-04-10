//
// Created by rolf on 28-2-23.
//

#ifndef PCOG_INCLUDE_PCOG_STABLESETMAXIMIZER_HPP
#define PCOG_INCLUDE_PCOG_STABLESETMAXIMIZER_HPP

#include "pcog/utilities/DenseGraph.hpp"
#include <random>
namespace pcog{
class StableSetMaximizer {
 public:
   explicit StableSetMaximizer(std::size_t seed);
   void maximizeRandomly(DenseSet& set, const DenseGraph& graph);
   void maximizeDegree(DenseSet& set, const DenseGraph& graph);
   void setSeed(std::size_t seed);
 private:
   std::mt19937 random_engine;
};
}
#endif // PCOG_INCLUDE_PCOG_STABLESETMAXIMIZER_HPP
