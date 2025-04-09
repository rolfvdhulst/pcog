//
// Created by rolf on 2-3-23.
//

#include "pcog/SafeDualWeights.hpp"
namespace pcog {
SafeDualWeights::weight_type
SafeDualWeights::fixedDualCost() {
   return 0;
}

SafeDualWeights::weight_type
SafeDualWeights::setCost(const DenseSet &set) const {
   weight_type sum = 0;
   for (const auto &node : set) {
      assert(sum + node_dual_value[node] >= sum); //checking for overflow
      sum += node_dual_value[node];
   }
   return sum;
}

SafeDualWeights::weight_type
SafeDualWeights::nodeCost(Node node) const {
   return node_dual_value[node];
}
SafeDualWeights::weight_type
SafeDualWeights::dualNodeCost(Node node) const {
   return nodeCost(node);
}

SafeDualWeights::SafeDualWeights(
    const std::vector<double> &t_weights) {

   auto result = toSafeWeights(t_weights);
   node_dual_value = result.first;
   scalar = result.second;
}

SafeDualWeights::weight_type
SafeDualWeights::getOne() const {
   return toSafeWeightUpperbound(1.0, scalar); //TODO: upper bound here or lower bound?
}

std::vector<SafeDualWeights::weight_type>
SafeDualWeights::weights() const {
   return node_dual_value;
}

SafeDualWeights::weight_type
SafeDualWeights::setCostUpperBound(const DenseSet &set) const {
   // Every double value can be just below the safe value n+1
   // so adding these doubles we need to add 1 for every element of the set = set.size()
   return setCost(set) + static_cast<SafeDualWeights::weight_type>(set.size());
}

SafeDualWeights::weight_type
SafeDualWeights::fullCost() const {
   weight_type sum = 0;
   for (const auto &dual_value : node_dual_value) {
      assert(sum + dual_value >= sum);
      sum += dual_value;
   }
   return sum;
}

double SafeDualWeights::scale() const { return scalar; }

std::vector<SafeDualWeights::weight_type> &
SafeDualWeights::mutable_weights() {
   return node_dual_value;
}
}