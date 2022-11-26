//
// Created by rolf on 22-11-22.
//

#ifndef PCOG_INCLUDE_PCOG_DEFINITIONS_HPP
#define PCOG_INCLUDE_PCOG_DEFINITIONS_HPP

#include <cstddef>
#include <limits>

namespace pcog {
using Node = std::size_t;
using degree_type = std::size_t;

constexpr static Node INVALID_NODE = std::numeric_limits<Node>::max();
} // namespace pcog
#endif // PCOG_INCLUDE_PCOG_DEFINITIONS_HPP
