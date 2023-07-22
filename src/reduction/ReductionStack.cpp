//
// Created by rolf on 5-7-23.
//

#include <ranges>

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
void ReductionStack::transformStableSet(DenseSet &set) const{
   for(const auto & reduction : reductions){
      std::visit([&set](const auto &arg) { arg.transformStableSet(set); },
                 reduction);
   }

}
void ReductionStack::clear() {
   reductions.clear();
   nFixedColors = 0;
}
void ReductionStack::newToOldColoring(NodeColoring &coloring) const {
   std::size_t i = reductions.size() -1;
   for(auto it = reductions.rbegin(); it != reductions.rend(); ++it){
      std::visit([&coloring](const auto& red){
         red.newToOldColoring(coloring);
      },*it);
      --i;
   }
}
}