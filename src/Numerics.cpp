//
// Created by rolf on 2-3-23.
//

#include "pcog/Numerics.hpp"
#include <cfenv>
#include <limits>
#include <cmath>

namespace pcog {
//TODO: consider Boost.NumericConversion for rounding double to SafeWeight
SafeWeight toSafeWeightUpperbound(double value, double scale){
   RoundingMode mode = fegetround();

   fesetround(FE_UPWARD);
   double result = value * scale;
   SafeWeight weight = std::ceil(result);

   fesetround(mode);
   return weight;
}


std::pair<std::vector<SafeWeight>, double>
toSafeWeights(const std::vector<double> &t_weights) {
   RoundingMode mode = fegetround();

   //Adjust floating point rounding mode so that any floating point errors in the following calculations can only occur downwards
   fesetround(FE_UPWARD);
   double maxSetWeight = static_cast<double>(t_weights.size());

   fesetround(FE_DOWNWARD);
   double scale = std::numeric_limits<SafeWeight>::max();
   scale /= maxSetWeight;

   std::vector<SafeWeight> safeWeights(t_weights.size());
   for (std::size_t i = 0; i < t_weights.size(); ++i) {
      double result = t_weights[i] * scale;
      SafeWeight weight = std::floor(result); //TODO: check if this cast is okay or not/dive into the stdlib or use Boost.NumericConversion
      safeWeights[i] = weight;
   }
   fesetround(mode); //reset the rounding mode again
   return {safeWeights,scale};
}
}