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
#include "pcog/utilities/Coloring.hpp"

namespace pcog{
using Reduction = std::variant<LowDegreeReduction,DominatedReduction,SimplicialReduction,
                               FoldDegreeTwoReduction,TwinDegreeThreeReduction,TwinDegreeThreeFoldReduction,CrownReduction,
                               TwoFixingReduction>;
class ReductionStack {
 public:
   void push(const Reduction& reduction);
   [[nodiscard]] std::size_t numFixedColors() const;
   void transformStableSet(DenseSet& set) const;
   void clear();
   void newToOldColoring(NodeColoring& oldColoring) const;
 private:
   std::vector<Reduction> reductions;
   std::size_t nFixedColors = 0;

};
}

#endif // PCOG_SRC_REDUCTION_REDUCTIONSTACK_HPP
