//
// Created by rolf on 4-7-23.
//

#ifndef PCOG_SRC_REDUCTION_REDUCTIONGRAPH_HPP
#define PCOG_SRC_REDUCTION_REDUCTIONGRAPH_HPP

#include "pcog/utilities/DenseGraph.hpp"
namespace pcog {
class DenseReductionGraph {
 public:
   explicit DenseReductionGraph(const DenseGraph& graph);
   void removeNode(Node node);
//   void removeStableSet(const DenseSet& set);
   [[nodiscard]] const DenseSet& neighbourhood(Node node) const;
   [[nodiscard]] const DenseSet& nodes() const;
   [[nodiscard]] std::size_t nodeDegree(Node node) const;
   [[nodiscard]] bool setIsStable(const DenseSet& set) const;
   [[nodiscard]] bool setIsStableMaximal(const DenseSet& set) const;
   [[nodiscard]] bool containsNode(Node node) const;

   [[nodiscard]] bool hasLowerBound() const;
   [[nodiscard]] std::size_t lowerBound() const;
   [[nodiscard]] const DenseSet& lowerBoundNodes() const;
 private:
   std::vector<DenseSet> adjacencyMatrix;
   std::vector<std::size_t> degrees;
   DenseSet presentNodes;

   DenseSet cliqueNodes;
   std::size_t lowerbound = 0;
};
}
#endif // PCOG_SRC_REDUCTION_REDUCTIONGRAPH_HPP
