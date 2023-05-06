//
// Created by rolf on 6-5-23.
//

#ifndef PCOG_SRC_GREEDYCOLORING_HPP
#define PCOG_SRC_GREEDYCOLORING_HPP

#include "utilities/DenseGraph.hpp"
#include "utilities/Coloring.hpp"

namespace pcog {

class GreedyColoring {
 public:
   explicit GreedyColoring(const DenseGraph &graph);
   [[nodiscard]] SetColoring run_degree_node() const;
   [[nodiscard]] SetColoring run_saturation_degree() const;

   [[nodiscard]] NodeColoring
   run_partial_sequential_coloring(const std::vector<Node> &ordering,
                                   std::size_t colors) const;
   [[nodiscard]] NodeColoring
   run_partial_saturation_coloring(const std::vector<Node> &ordering,
                                   std::size_t colors) const;

 private:
   const DenseGraph &graph;
};
}
#endif // PCOG_SRC_GREEDYCOLORING_HPP
