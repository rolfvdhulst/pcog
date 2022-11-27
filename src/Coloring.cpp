//
// Created by rolf on 26-11-22.
//

#include "pcog/Coloring.hpp"
#include "pcog/DenseGraph.hpp"
// #include "pcog/Preprocessing.h"
namespace pcog {

SetColoring::SetColoring(const NodeColoring &t_nodeColoring)
    : m_colors(t_nodeColoring.num_colors,
               DenseSet(t_nodeColoring.color.size())) {
   for (Node node = 0; node < t_nodeColoring.color.size(); node++) {
      if (t_nodeColoring.color[node] != INVALID_COLOR) {
         m_colors[t_nodeColoring.color[node]].add(node);
      }
   }
}

bool SetColoring::isValid(const DenseGraph &t_graph) const {
   for (const auto &color : m_colors) {
      assert(t_graph.numNodes() == color.capacity());
      if (!t_graph.setIsStable(color)) {
         return false;
      }
   }
   // Check if all nodes are colored
   DenseSet covered(t_graph.numNodes());
   for (const auto &color : m_colors) {
      covered.inplaceUnion(color);
   }
   return covered.full();
}
void SetColoring::addColor(const DenseSet &t_set) { m_colors.push_back(t_set); }
std::size_t SetColoring::numColors() const { return m_colors.size(); }
const std::vector<DenseSet> &SetColoring::colors() const { return m_colors; }

NodeColoring::NodeColoring(std::size_t num_nodes, const SetColoring &coloring)
    : color(num_nodes, INVALID_COLOR), num_colors{coloring.m_colors.size()} {
   std::size_t index = 0;
   for (const auto &setColor : coloring.m_colors) {
      assert(setColor.capacity() == num_nodes);
      for (const Node &node : setColor) {
         color[node] = index;
      }
      index++;
   }
}
NodeColoring::NodeColoring(std::size_t num_nodes)
    : color(num_nodes, INVALID_COLOR), num_colors{0} {}

// TODO: move to preprocessing
/*
NodeColoring NodeColoring::extend(const PreprocessedMap &map,
                                  const DenseGraph &originalGraph) const {
   NodeColoring newColoring(map.oldToNewIDs.size());
   newColoring.num_colors = num_colors;
   for (std::size_t node = 0; node < color.size(); ++node) {
      newColoring.color[map.newToOldIDs[node]] = color[node];
   }
   if (!map.removed_nodes.empty()) {
      DenseSet unused_colors(num_colors);
      for (auto it = map.removed_nodes.rbegin(); it != map.removed_nodes.rend();
           it++) {
         const PreprocessedNode &preprocessedNode = *it;
         if (preprocessedNode.removedReason() ==
             PreprocessedReason::LOW_DEGREE) {
            assert(newColoring.color[preprocessedNode.removedNode()] >=
                   num_colors);
            unused_colors.setAll();
            Node removed_node = preprocessedNode.removedNode();
            for (const auto &node : originalGraph.neighbourhood(removed_node)) {
               if (newColoring.color[node] <
                   num_colors) { // the removed node may have some other removed
                                 // nodes in its neighbourhood which were
                                 // removed earlier
                  unused_colors.remove(newColoring.color[node]);
               }
            }
            std::size_t assigned_color = unused_colors.first();
            newColoring.color[removed_node] = assigned_color;
            assert(newColoring.color[preprocessedNode.removedNode()] <
                   num_colors);
         } else if (preprocessedNode.removedReason() ==
                    PreprocessedReason::DOMINATED_NODE) {
            // in case the node was dominated, we can extend it by using the
            // same color
            assert(newColoring.color[preprocessedNode.removedNode()] >=
                   num_colors);
            assert(newColoring.color[preprocessedNode.dominatingNode()] <
                   num_colors);
            newColoring.color[preprocessedNode.removedNode()] =
                newColoring.color[preprocessedNode.dominatingNode()];
            assert(newColoring.color[preprocessedNode.removedNode()] <
                   num_colors);
         }
      }
   }

   assert(hasNoInvalidNodes());
   return newColoring;
}
*/

bool NodeColoring::hasNoInvalidNodes() const {
   return std::all_of(color.begin(), color.end(),
                      [&](const Color &node_color) -> bool {
                         return node_color < num_colors;
                      }); // TODO: weird way to check invalid colors?
}

std::vector<std::size_t> NodeColoring::getStableSetSizes() const {
   std::vector<std::size_t> sizes(num_colors, 0);
   for (unsigned long node_color : color) {
      if (node_color < num_colors) {
         sizes[node_color]++;
      }
   }
   return sizes;
}

std::size_t NodeColoring::numUncolored() const {
   std::size_t sum = 0;
   for (unsigned long node_color : color) {
      if (node_color >= num_colors) { // TODO: why this instead of checking
                                      // equality with INVALID_COLOR?
         sum++;
      }
   }
   return sum;
}
const Color &NodeColoring::operator[](Node node) const { return color[node]; }
Color &NodeColoring::operator[](Node node) { return color[node]; }
void NodeColoring::setNumColors(std::size_t numColors) {
   num_colors = numColors;
}
std::size_t NodeColoring::numColors() const {
   return num_colors;
}
degree_type NodeColoring::numNodes() const { return color.size(); }
} // namespace pcog