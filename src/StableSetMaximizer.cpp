//
// Created by rolf on 28-2-23.
//

#include "pcog/StableSetMaximizer.hpp"
namespace  pcog {
StableSetMaximizer::StableSetMaximizer(std::size_t seed)
    : random_engine(seed) {}

void StableSetMaximizer::maximizeDegree(DenseSet &set,
                                        const DenseGraph &graph) {
   assert(set.capacity() == graph.numNodes());
   std::vector<std::pair<Node, std::size_t>> node_degrees(graph.numNodes());
   for (Node i = 0; i < graph.numNodes(); ++i) {
      node_degrees[i].first = i;
      node_degrees[i].second = graph.neighbourhood(i).size();
   }
   std::sort(node_degrees.begin(), node_degrees.end(),
             [](const std::pair<Node, std::size_t> &a,
                const std::pair<Node, std::size_t> &b) {
                return a.second <
                       b.second; // TODO: is this ascending or descending
             });
   for (const auto &node_degree_pair : node_degrees) {
      if (!set.contains(node_degree_pair.first)) {
         if (graph.neighbourhood(node_degree_pair.first)
                 .intersection(set)
                 .empty()) {
            set.add(node_degree_pair.first);
         }
      }
   }
}

void StableSetMaximizer::maximizeRandomly(DenseSet &set,
                                          const DenseGraph &graph) {
   std::vector<Node> nodes(graph.numNodes(), 0);
   std::iota(nodes.begin(), nodes.end(), 0);
   std::shuffle(nodes.begin(), nodes.end(), random_engine);
   for (const auto &node : nodes) {
      if (!set.contains(node)) {
         if (graph.neighbourhood(node).intersection(set).empty()) {
            set.add(node);
         }
      }
   }
}

void StableSetMaximizer::setSeed(std::size_t seed) { random_engine.seed(seed); }
}