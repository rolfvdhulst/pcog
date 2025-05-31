//
// Created by rolf on 26-4-23.
//

#ifndef PCOG_SRC_MWSS_AUGMENTINGSEARCH_HPP
#define PCOG_SRC_MWSS_AUGMENTINGSEARCH_HPP

#include "pcog/utilities/DenseGraph.hpp"
#include "pcog/utilities/Numerics.hpp"
#include "pcog/utilities/SparseSet.hpp"

namespace pcog {
class AugmentingSearch {
 public:
   AugmentingSearch(const std::vector<SafeWeight> &greedyWeights,
                    const DenseGraph &graph);
   using Weight = SafeWeight;

   void setSolution(const DenseSet &set);
   void clearSolution();

   void updateWeights(const std::vector<SafeWeight> &greedyWeights);
   [[nodiscard]] bool solutionIsMaximal() const;
   [[nodiscard]] DenseSet getDenseSolution() const;
   bool searchImprovements();
   [[nodiscard]] Weight totalWeight() const;
   [[nodiscard]] std::size_t currentSolutionSize() const;

   [[nodiscard]] bool nodeFree(Node node) const;
   void addNodeToSolution(Node node);
   void removeNodeFromSolution(Node node);
   [[nodiscard]] std::vector<Node> decreasingWeightOrdering() const;
   [[nodiscard]] std::vector<Node> decreasingRatioOrdering() const;
   [[nodiscard]] std::vector<Node> decreasingMaxRatioOrdering() const;

   enum class GreedyStrategy { ByWeight, ByDynamicSurplus, BySurplus, ByRatio, ByMaxRatio};

   std::size_t greedyByWeight(long rotateOrderBy);
   std::size_t greedyByDynamicSurplus(long rotateOrderBy);
   std::size_t greedyByStaticSurplus(long rotateOrderBy);
   std::size_t greedyByRatio(long rotateOrderBy);
   std::size_t greedyByMaxRatio(long rotateOrderBy);
   template <AugmentingSearch::GreedyStrategy strategy>
   std::size_t greedyAlgorithm(long rotateOrderBy = 0) {
      if constexpr (strategy == AugmentingSearch::GreedyStrategy::ByWeight) {
         return greedyByWeight(rotateOrderBy);
      }
      if constexpr (strategy ==
                    AugmentingSearch::GreedyStrategy::ByDynamicSurplus) {
         return greedyByDynamicSurplus(rotateOrderBy);
      }
      if constexpr (strategy == AugmentingSearch::GreedyStrategy::ByRatio) {
         return greedyByRatio(rotateOrderBy);
      }
      if constexpr (strategy == AugmentingSearch::GreedyStrategy::ByMaxRatio) {
         return greedyByMaxRatio(rotateOrderBy);
      }
      return greedyByStaticSurplus(rotateOrderBy);
   }

   template <AugmentingSearch::GreedyStrategy strategy>
   std::size_t doTwoImprovementsWithGreedy() {
      std::size_t totalSwaps = 0;
      std::size_t iterSwaps;
      do {
         iterSwaps = 0;
         for (Node node = 0; node < graph.numNodes(); ++node) {
            if (nodeInSolution(
                    node)) { // TODO: iterate over a sparse set rather than over
                             // all nodes with if NodeInSolution()?
               if (searchTwoImprovement(node)) {
                  --node; // rerun iteration on the same node again
                  iterSwaps++;
               }
            }
         }
         if (iterSwaps) {
            std::size_t greedyExtra = greedyAlgorithm<strategy>();
            iterSwaps += greedyExtra;
            totalSwaps += iterSwaps;
         }
      } while (iterSwaps);
      return totalSwaps;
   }
   std::size_t doTwoKImprovementsWithGreedy();
   template <AugmentingSearch::GreedyStrategy strategy>
   std::vector<DenseSet>
   repeatedGreedyLocalSearch(AugmentingSearch::Weight cutoff,
                             bool reduceNumStarts = false) {
      std::vector<DenseSet> solutions;
      long numStarts = long(table.size());
      if (reduceNumStarts) {
         numStarts = std::max(numStarts / 20l, 1l);
      }
      // long lastImprovingStart = -1;
      long firstValidStart = -1;

      for (long rotateBy = 0; rotateBy < numStarts; rotateBy++) {
         clearSolution();
         std::size_t changes = greedyAlgorithm<strategy>(rotateBy);
         while (changes) {
            changes = doTwoImprovementsWithGreedy<strategy>();

            std::size_t change = doTwoKImprovements();
            changes += change;
            if (change) {
               greedyAlgorithm<strategy>(0);
            }
         }
         if (totalWeight() > cutoff) {
            solutions.push_back(getDenseSolution());
            // lastImprovingStart = rotateBy;
            if (firstValidStart == -1) {
               firstValidStart = rotateBy;
            }
            numStarts *= 2;
            numStarts /= 3;
         }
      }
      return solutions;
   }

 private:
   bool doTwoImprovements();
   bool doTwoKImprovements();
   bool searchTwoKImprovement(Node node);
   bool searchTwoImprovement(Node node);
   [[nodiscard]] bool nodeInSolution(Node node) const;
   [[nodiscard]] std::size_t getTightness(Node node) const;
   [[nodiscard]] Weight getWeight(Node node) const;

   enum class GreedySolutionStatus { FREE = 0, SOLUTION = 1, NONFREE = 2 };
   struct GreedySearchTableNode {
      Weight greedy_weight; // TODO: listen to weightfunction
      std::size_t tightness;
      GreedySolutionStatus status;
   };
   std::vector<GreedySearchTableNode> table;

   std::size_t num_solution_nodes;
   std::size_t num_free_nodes;
   const DenseGraph &graph;
   SparseSet temporary_set;
};

} // namespace pcog

#endif // PCOG_SRC_MWSS_AUGMENTINGSEARCH_HPP
