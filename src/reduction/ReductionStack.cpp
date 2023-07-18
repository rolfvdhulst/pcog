//
// Created by rolf on 5-7-23.
//

#include "pcog/reduction/ReductionStack.hpp"
namespace pcog {
void ReductionStack::push(const Reduction& reduction) {
   if(std::holds_alternative<SimplicialReduction>(reduction)){
      nFixedColors += 1;
   }else if(std::holds_alternative<FoldDegreeTwoReduction>(reduction)){
      nFixedColors += 1;
   }else if(std::holds_alternative<TwinDegreeThreeReduction>(reduction)){
      nFixedColors += 2;
   }else if(std::holds_alternative<CrownReduction>(reduction)){
      const auto& crownReduction = std::get<CrownReduction>(reduction);
      nFixedColors += crownReduction.fixedSets.size();
   }else if(std::holds_alternative<TwoFixingReduction>(reduction)){
      nFixedColors += 2;
   }
   reductions.push_back(reduction);
}
std::size_t ReductionStack::numFixedColors() const {
   return nFixedColors;
}
}