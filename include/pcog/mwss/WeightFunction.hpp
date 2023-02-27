//
// Created by rolf on 27-2-23.
//

#ifndef PCOG_INCLUDE_PCOG_MWSS_WEIGHTFUNCTION_HPP
#define PCOG_INCLUDE_PCOG_MWSS_WEIGHTFUNCTION_HPP

#include "pcog/DenseSet.hpp"
#include <concepts>

namespace pcog {
enum class WeightFunctionType {
   UNIFORM_COST, // all nodes have the same cost/weight
   LINEAR, // The weight function is linear over the nodes (e.g. every node gets
           // one weight)
   SUBADDITIVE // The weight function is nonlinear but subadditive over the
               // nodes
};

template <typename T>
concept WeightFunction = requires(T a, const DenseSet &set, const Node &node) {
                            typename T::weight_type;
                            {
                               T::function_type
                            } -> std::convertible_to<WeightFunctionType>;
                            {
                               a.nodeCost(node)
                            } -> std::convertible_to<typename T::weight_type>;
                            {
                               a.setCost(set)
                            } -> std::convertible_to<typename T::weight_type>;
                            {
                               a.dualNodeCost(node)
                            } -> std::convertible_to<typename T::weight_type>;
                            {
                               a.fixedDualCost()
                            } -> std::convertible_to<typename T::weight_type>;
                         };
} // namespace pcog
#endif // PCOG_INCLUDE_PCOG_MWSS_WEIGHTFUNCTION_HPP
