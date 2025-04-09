//
// Created by rolf on 1-8-24.
//

#ifndef PCOG_SPARSEGRAPH_HPP
#define PCOG_SPARSEGRAPH_HPP

#include <vector>
class SparseGraph {
 public:
   SparseGraph() = default;

 private:
   std::vector<std::vector<std::size_t>> m_adjacentNodes;
};

#endif // PCOG_SPARSEGRAPH_HPP
