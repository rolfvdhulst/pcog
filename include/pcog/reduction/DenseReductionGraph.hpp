//
// Created by rolf on 4-7-23.
//

#ifndef PCOG_SRC_REDUCTION_REDUCTIONGRAPH_HPP
#define PCOG_SRC_REDUCTION_REDUCTIONGRAPH_HPP

#include "pcog/utilities/DenseGraph.hpp"
namespace pcog {
struct InducedGraph{
   DenseGraph graph;
   NodeMap newToOld;
   NodeMap oldToNew;
};
class DenseReductionGraph {
 public:
   explicit DenseReductionGraph(const DenseGraph& graph);
   void removeNode(Node node);
   void removeStableSet(const DenseSet& set);
   void removeNodeEdges(Node node, const DenseSet& otherNodes);
   [[nodiscard]] DenseSet complementNeighbourhood(Node node) const;
   [[nodiscard]] const DenseSet& neighbourhood(Node node) const;
   [[nodiscard]] const DenseSet& nodes() const;
   [[nodiscard]] std::size_t numNodes() const;
   [[nodiscard]] std::size_t nodeDegree(Node node) const;
   [[nodiscard]] bool setIsStable(const DenseSet& set) const;
   [[nodiscard]] bool setIsStableMaximal(const DenseSet& set) const;
   [[nodiscard]] bool containsNode(Node node) const;
   [[nodiscard]] InducedGraph currentGraph() const;


   [[nodiscard]] bool hasLowerBound() const;
   [[nodiscard]] std::size_t lowerBound() const;
   [[nodiscard]] const DenseSet& lowerBoundNodes() const;
   void setLowerBoundClique(const DenseSet& clique, std::size_t numCliqueNodes);
 private:
   void removeGraphNode(Node node);
   std::vector<DenseSet> adjacencyMatrix;
   std::vector<std::size_t> degrees;
   DenseSet presentNodes;
   std::size_t nodeCount;

   DenseSet cliqueNodes;
   std::size_t lowerbound = 0;
};
}
#endif // PCOG_SRC_REDUCTION_REDUCTIONGRAPH_HPP
