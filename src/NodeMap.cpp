//
// Created by rolf on 23-11-22.
//

#include "pcog/NodeMap.hpp"
#include "pcog/DenseSet.hpp"

#include <cassert>
#include <numeric>
namespace pcog {
NodeMap::NodeMap(degree_type size) : permutation_map(size, INVALID_NODE) {}

void NodeMap::setIdentity() {
   std::iota(permutation_map.begin(), permutation_map.end(), 0);
}

bool NodeMap::isIdentity() const {
   Node val = 0;
   for (const auto &node : permutation_map) {
      if (node != val) {
         return false;
      }
      val++;
   }
   return true;
}

Node &NodeMap::operator[](Node node) { return permutation_map[node]; }

NodeMap NodeMap::inverse(const NodeMap &map, degree_type size) {
   NodeMap inverse;
   inverse.permutation_map = std::vector(size, INVALID_NODE);
   for (std::size_t i = 0; i < map.permutation_map.size(); i++) {
      Node value = map.permutation_map[i];
      assert(value < size);
      inverse.permutation_map[value] = i;
   }
   return inverse;
}

const Node &NodeMap::operator[](Node node) const {
   return permutation_map[node];
}

void NodeMap::transform(const DenseSet &set, DenseSet &toStoreIn) const {
   for (const auto &elem : set) {
      if (permutation_map[elem] != INVALID_NODE) {
         assert(permutation_map[elem] < toStoreIn.capacity());
         toStoreIn.add(permutation_map[elem]);
      }
   }
}

std::size_t NodeMap::size() const { return permutation_map.size(); }

void NodeMap::extend(Node node) { permutation_map.push_back(node); }
} // namespace pcog