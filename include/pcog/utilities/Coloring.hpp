//
// Created by rolf on 26-11-22.
//

#ifndef PCOG_SRC_COLORING_HPP
#define PCOG_SRC_COLORING_HPP

#include "DenseSet.hpp"
#include "Definitions.hpp"
///\file We use two classes to store colorings. NodeColoring is a class which
/// stores a single color assigned to each node, whereas SetColoring can assign
/// multiple colors to a node and may be more efficient/practical to use when
/// comparing with other sets

namespace pcog {

class DenseGraph;
class NodeColoring;


using Color = Node;
static constexpr Color INVALID_COLOR = INVALID_NODE;

class SetColoring {
   friend class NodeColoring;
 public:
   ///Constructs an empty coloring with no color classes.
   SetColoring() = default;
   ///Constructs a set coloring from the given node coloring
   explicit SetColoring(const NodeColoring &t_nodeColoring);

   /// Adds a color (independent set) to the coloring
   /// \param t_set The nodes given the same color.
   void addColor(const DenseSet& t_set);

   /// Returns a vector with the color classes
   /// \return
   [[nodiscard]] const std::vector<DenseSet>& colors() const;
   std::vector<DenseSet>& colors();

   /// \return The number of colors in the current coloring
   [[nodiscard]] std::size_t numColors() const;

   /// Checks if the set coloring is valid, e.g. each set is stable and all nodes are colored.
   /// \param t_graph Graph to check the coloring against
   /// \return True if the coloring is valid, false otherwise.
   [[nodiscard]] bool isValid(const DenseGraph &t_graph) const;
 private:
   std::vector<DenseSet> m_colors;
};

class NodeColoring {
   friend class SetColoring;
 public:
   /// Constructs a node coloring for the given number of nodes.
   /// Initially no node is colored and each node is given INVALID_COLOR
   /// \param num_nodes the number of nodes this coloring has.
   explicit NodeColoring(std::size_t num_nodes);

   /// Constructs a node coloring from a set coloring.
   /// The number of nodes passed should be equal to the capacity of all colors in coloring
   /// \param num_nodes The number of nodes the coloring has
   /// \param coloring The set coloring
   NodeColoring(std::size_t num_nodes, const SetColoring &coloring);

   /// Accesses the color of the given node. Allows it to be overwritten.
   /// \param node The node to get the color for.
   /// \return A reference to the color of the node.
   Color& operator[](Node node);
   /// Read the color of a given node.
   /// \param node The node to get the color for.
   /// \return A constant reference to the color of the node.
   const Color& operator[](Node node) const;

   /// Sets the number of colors of this coloring. This is mostly important for this class to efficiently interoperate with other classes.
   /// \param numColors The number of colors in the current coloring.
   void setNumColors(std::size_t numColors);

   /// \return The number of colors in the current coloring
   [[nodiscard]] std::size_t numColors() const;


   /// \return the number of nodes which this coloring has.
   [[nodiscard]] degree_type numNodes() const;

   /// Checks if the coloring has any nodes with INVALID_COLOR
   /// \return true if the coloring has any node with INVALID_COLOR
   [[nodiscard]] bool hasNoInvalidNodes() const;
   /// Returns a vector with the sizes of each stable set of only the valid colors
   /// \return
   [[nodiscard]] std::vector<std::size_t> getStableSetSizes() const;
   ///
   /// \return The number of uncolored nodes in this coloring
   [[nodiscard]] std::size_t numUncolored() const;
 private:
   std::vector<std::size_t> color; // color[node] gives color of each node
   std::size_t num_colors;
};
} // namespace pcog
#endif // PCOG_SRC_COLORING_HPP
