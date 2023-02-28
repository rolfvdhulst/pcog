//
// Created by rolf on 23-11-22.
//

#ifndef PCOG_SRC_NODEMAP_HPP
#define PCOG_SRC_NODEMAP_HPP

#include "Definitions.hpp"
#include <vector>
namespace pcog {

class DenseSet;

/// A class which stores a function which maps nodes to other nodes in a dense
/// manner.
class NodeMap {
   std::vector<Node> permutation_map;

 public:
   NodeMap() = default;
   explicit NodeMap(degree_type size);
   /// Constructs the inverse of this function, where the size of the input
   /// domain is given by size.
   static NodeMap inverse(const NodeMap &map, degree_type size);

   static NodeMap identity(degree_type size);
   /// Gets the number of inputs allowed for this function/ the maximal node
   /// this function has as input
   [[nodiscard]] std::size_t size() const;

   /// Transform the given set given the node map and stores the result
   ///  \param set Set to transform
   ///  \param toStoreIn Out parameter to store the transformed set in. It is
   ///  assumed enough space is allocated to hold the result, user should ensure
   ///  this.
   void transform(const DenseSet &set, DenseSet &toStoreIn) const;

   ///Empties the node-map
   void clear();

   /// Query/modify the function for a single node
   ///  \param node
   ///  \return A reference to the queried output.
   Node &operator[](Node node);
   const Node &operator[](Node node) const;
   /// Resets this function to an identity function
   void setIdentity();
   /// \return True if this function is an identity function
   [[nodiscard]] bool isIdentity() const;
   /// Increase the size by one and adds a new output to the current map
   /// \param node To append to the map
   void extend(Node node);
};
} // namespace pcog

#endif // PCOG_SRC_NODEMAP_HPP
