//
// Created by rolf on 5-7-23.
//

#ifndef PCOG_SRC_REDUCTION_REDUCTIONSTACK_HPP
#define PCOG_SRC_REDUCTION_REDUCTIONSTACK_HPP
#include <variant>
#include "LowDegreeReduction.hpp"
#include "DominatedReduction.hpp"
#include "SimplicialReduction.hpp"
#include "FoldDegreeTwoReduction.hpp"
#include "TwinDegreeThreeReduction.hpp"

namespace pcog{
using Reduction = std::variant<LowDegreeReduction,DominatedReduction,SimplicialReduction,
                               FoldDegreeTwoReduction,TwinDegreeThreeReduction>;
class ReductionStack {
 public:
   void push(const Reduction& reduction);
 private:
   std::vector<Reduction> reductions;

   std::size_t removedNodes = 0;
};
}

#endif // PCOG_SRC_REDUCTION_REDUCTIONSTACK_HPP
