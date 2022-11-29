//
// Created by rolf on 27-11-22.
//

#include "pcog/Preprocessing.hpp"
#include <utility>

using namespace pcog;
PreprocessingResult::PreprocessingResult(DenseGraph graph, PreprocessedMap map)
    : graph{std::move(graph)}, map{std::move(map)} {}

PreprocessingResult preprocessOriginalGraph(const DenseGraph &graph,
                                            const DenseSet &clique) {

   std::vector<PreprocessedNode> removed_nodes;
   DenseSet present_nodes(graph.numNodes(), true);
   assert(present_nodes.full());

   std::size_t previous_round_size = 1;
   while (previous_round_size != removed_nodes.size()) {
      previous_round_size = removed_nodes.size();
      while (true) {
         auto low_degree_removed =
             removeLowDegreeVerticesClique(graph, present_nodes, clique);
         if (low_degree_removed.empty()) {
            break;
         }
         removed_nodes.insert(removed_nodes.end(), low_degree_removed.begin(),
                              low_degree_removed.end());
      }

      while (true) {
         auto dominated_removed =
             removeDominatedVerticesClique(graph, present_nodes, clique);
         if (dominated_removed.empty()) {
            break;
         }
         removed_nodes.insert(removed_nodes.end(), dominated_removed.begin(),
                              dominated_removed.end());
      }
   }

   auto [newGraph, newToOld] = graph.nodeInducedSubgraph(present_nodes);
   NodeMap oldToNew = NodeMap::inverse(newToOld, graph.numNodes());

   return PreprocessingResult(
       newGraph, PreprocessedMap(removed_nodes, newToOld, oldToNew));
}

std::vector<PreprocessedNode>
removeLowDegreeVerticesClique(const DenseGraph &graph, DenseSet &present_nodes,
                              const DenseSet &clique) {
   std::size_t lower_bound = clique.size();
   DenseSet iterate_nodes = present_nodes;
   iterate_nodes.inplaceDifference(clique);

   DenseSet removed_nodes(graph.numNodes());
   assert(removed_nodes.empty());

   for (const auto &node : iterate_nodes) {
      if (graph.neighbourhood(node).intersection(present_nodes).size() <
          lower_bound) {
         removed_nodes.add(node);
      }
   }

   present_nodes.inplaceDifference(removed_nodes);
   std::vector<PreprocessedNode> removed;
   for (const auto &node : removed_nodes) {
      removed.emplace_back(
          PreprocessedNode(node, PreprocessedReason::LOW_DEGREE));
   }
   return removed;
}

std::vector<PreprocessedNode>
removeDominatedVerticesClique(const DenseGraph &graph, DenseSet &present_nodes,
                              const DenseSet &clique) {
   DenseSet iterate_nodes = present_nodes;
   iterate_nodes.inplaceDifference(clique);
   std::vector<PreprocessedNode> removed_nodes;
   for (const auto &node_a : iterate_nodes) {
      for (const auto &node_b : present_nodes) {
         if (node_a == node_b) {
            continue;
         }
         const DenseSet masked_neighbourhood_a =
             graph.neighbourhood(node_a).intersection(present_nodes);
         const DenseSet masked_neighbourhood_b =
             graph.neighbourhood(node_b).intersection(present_nodes);
         if (masked_neighbourhood_a.isSubsetOf(masked_neighbourhood_b)) {
            removed_nodes.emplace_back(
                node_a, PreprocessedReason::DOMINATED_NODE, node_b);
            present_nodes.remove(node_a);
            break;
         }
      }
   }
   return removed_nodes;
}

PreprocessedNode::PreprocessedNode(Node t_node, PreprocessedReason t_reason,
                                   Node t_dominatedByNode)
    : m_node{t_node}, m_reason{t_reason}, m_dominated_by{t_dominatedByNode} {}

Node PreprocessedNode::removedNode() const { return m_node; }

PreprocessedReason PreprocessedNode::removedReason() const { return m_reason; }

Node PreprocessedNode::dominatingNode() const { return m_dominated_by; }

PreprocessedMap::PreprocessedMap(std::vector<PreprocessedNode> removed_nodes,
                                 NodeMap newToOld, NodeMap oldToNew)
    : removed_nodes{std::move(removed_nodes)}, newToOldIDs{std::move(newToOld)},
      oldToNewIDs{std::move(oldToNew)} {}

std::vector<PreprocessedNode> removeLowDegreeVertices(const DenseGraph &graph,
                                                      DenseSet &present_nodes,
                                                      std::size_t lower_bound) {

   DenseSet removed_nodes(graph.numNodes());
   assert(removed_nodes.empty());

   for (const auto &node : present_nodes) {
      if (graph.neighbourhood(node).intersection(present_nodes).size() <
          lower_bound - 1) { // Note the -1, otherwise we could remove clique
                             // nodes by accident!
         removed_nodes.add(node);
      }
   }

   present_nodes.inplaceDifference(removed_nodes);
   std::vector<PreprocessedNode> removed;
   for (const auto &node : removed_nodes) {
      removed.emplace_back(
          PreprocessedNode(node, PreprocessedReason::LOW_DEGREE));
   }
   return removed;
}

std::vector<PreprocessedNode>
removeStrictlyDominatedVertices(const DenseGraph &graph,
                                DenseSet &present_nodes) {

   std::vector<PreprocessedNode> removed_nodes;
   for (const auto &node_a : present_nodes) {
      for (const auto &node_b : present_nodes) {
         if (node_a == node_b) {
            continue;
         }
         const DenseSet masked_neighbourhood_a =
             graph.neighbourhood(node_a).intersection(present_nodes);
         const DenseSet masked_neighbourhood_b =
             graph.neighbourhood(node_b).intersection(present_nodes);
         if (masked_neighbourhood_a.isProperSubsetOf(masked_neighbourhood_b)) {
            removed_nodes.emplace_back(
                node_a, PreprocessedReason::DOMINATED_NODE, node_b);
            break;
         }
      }
   }
   for (const auto &processed_node : removed_nodes) {
      present_nodes.remove(processed_node.removedNode());
   }
   return removed_nodes;
}

NodeColoring extendColoring(const NodeColoring& coloring,
                            const PreprocessedMap &map,
                            const DenseGraph &originalGraph){
   NodeColoring newColoring(map.oldToNewIDs.size());
   newColoring.setNumColors(coloring.numColors());

   for (std::size_t node = 0; node < coloring.numNodes(); ++node) {
      newColoring[map.newToOldIDs[node]] = coloring[node];
   }
   if (!map.removed_nodes.empty()) {
      DenseSet unused_colors(coloring.numColors());
      for (auto it = map.removed_nodes.rbegin(); it != map.removed_nodes.rend();
           it++) {
         const PreprocessedNode &preprocessedNode = *it;
         if (preprocessedNode.removedReason() == PreprocessedReason::LOW_DEGREE) {
            assert(newColoring[preprocessedNode.removedNode()] >= coloring.numColors());
            unused_colors.setAll();
            Node removed_node = preprocessedNode.removedNode();
            for (const auto &node : originalGraph.neighbourhood(removed_node)) {
               if (newColoring[node] < coloring.numColors()) { // the removed node may have some other removed
                                 // nodes in its neighbourhood which were
                                 // removed earlier
                  unused_colors.remove(newColoring[node]);
               }
            }
            std::size_t assigned_color = unused_colors.first();
            newColoring[removed_node] = assigned_color;
            assert(newColoring[preprocessedNode.removedNode()] <
                   coloring.numColors());
         } else if (preprocessedNode.removedReason() ==
                    PreprocessedReason::DOMINATED_NODE) {
            // in case the node was dominated, we can extend it by using the
            // same color
            assert(newColoring[preprocessedNode.removedNode()] >=
                   coloring.numColors());
            assert(newColoring[preprocessedNode.dominatingNode()] <
                   coloring.numColors());
            newColoring[preprocessedNode.removedNode()] =
                newColoring[preprocessedNode.dominatingNode()];
            assert(newColoring[preprocessedNode.removedNode()] <
                   coloring.numColors());
         }
      }
   }

   assert(newColoring.hasNoInvalidNodes());
   return newColoring;
}
