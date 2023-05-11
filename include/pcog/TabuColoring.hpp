//
// Created by rolf on 6-5-23.
//

#ifndef PCOG_SRC_TABUCOLORING_HPP
#define PCOG_SRC_TABUCOLORING_HPP
#include "utilities/Coloring.hpp"
#include "utilities/DenseGraph.hpp"

namespace pcog {
class TabuColoring {

 public:
   explicit TabuColoring(const DenseGraph &graph);
   void setMaxIterations(std::size_t iters);
   void setGamma(double gamma);
   void setTabuBase(std::size_t numBase);
   std::optional<NodeColoring> run(NodeColoring initialColoring);
   [[nodiscard]] std::size_t numIterations() const;
   static long numViolatedEdges(const DenseGraph &graph,
                                const NodeColoring &coloring);

 private:
   const DenseGraph &graph;
   std::size_t max_iterations = 100'000;
   std::size_t tabu_base = 50;
   double tabu_gamma = 0.9;

   std::size_t done_iterations = 0;
};
}
#endif // PCOG_SRC_TABUCOLORING_HPP
