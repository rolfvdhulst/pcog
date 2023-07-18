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
#include "CrownReduction.hpp"
#include "TwoFixingReduction.hpp"

namespace pcog{
using Reduction = std::variant<LowDegreeReduction,DominatedReduction,SimplicialReduction,
                               FoldDegreeTwoReduction,TwinDegreeThreeReduction,CrownReduction,
                               TwoFixingReduction>;
class ReductionStack {
 public:
   void push(const Reduction& reduction);
   std::size_t numFixedColors() const;
 private:
   std::vector<Reduction> reductions;
   std::size_t nFixedColors = 0;

};
}

#endif // PCOG_SRC_REDUCTION_REDUCTIONSTACK_HPP
