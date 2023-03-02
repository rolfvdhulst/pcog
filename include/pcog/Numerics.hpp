//
// Created by rolf on 2-3-23.
//

#ifndef PCOG_SRC_NUMERICS_HPP
#define PCOG_SRC_NUMERICS_HPP

#include <vector>
#include <cfenv>

namespace pcog {
using SafeWeight = int; // TODO: consider long?
using RoundingMode = int;
SafeWeight toSafeWeightUpperbound(double value, double scale);
std::pair<std::vector<SafeWeight>, double>
toSafeWeights(const std::vector<double> &t_weights);
} // namespace pcog

#endif // PCOG_SRC_NUMERICS_HPP
