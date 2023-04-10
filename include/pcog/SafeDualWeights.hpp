//
// Created by rolf on 2-3-23.
//

#ifndef PCOG_INCLUDE_PCOG_SAFEDUALWEIGHTS_HPP
#define PCOG_INCLUDE_PCOG_SAFEDUALWEIGHTS_HPP

#include "pcog/utilities/Numerics.hpp"
#include "mwss/WeightFunction.hpp"

namespace pcog {
class SafeDualWeights {
 public:
   SafeDualWeights(const std::vector<double> &t_weights);
   using weight_type = SafeWeight;
   static constexpr WeightFunctionType function_type =
       WeightFunctionType::LINEAR;
   [[nodiscard]] weight_type nodeCost(Node node) const;
   [[nodiscard]] weight_type setCost(const DenseSet &set) const;
   [[nodiscard]] weight_type setCostUpperBound(const DenseSet &set) const;
   [[nodiscard]] weight_type fullCost() const;
   [[nodiscard]] double scale() const;
   static weight_type fixedDualCost();
   [[nodiscard]] weight_type dualNodeCost(Node node) const;
   [[nodiscard]] weight_type getOne() const;

   [[nodiscard]] std::vector<weight_type> weights() const;
   std::vector<weight_type> &mutable_weights();

 private:
   std::vector<weight_type> node_dual_value;
   double scalar;
};
static_assert(WeightFunction<SafeDualWeights>);

}

#endif // PCOG_INCLUDE_PCOG_SAFEDUALWEIGHTS_HPP
