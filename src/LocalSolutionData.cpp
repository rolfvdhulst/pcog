//
// Created by rolf on 2-6-23.
//

#include "pcog/LocalSolutionData.hpp"

namespace pcog{
LocalSolutionData::LocalSolutionData() :
                                         m_numStartGlobalVariables{0},
                                         m_lpIterations{0},
                                         m_pricingIterations{0}{

}
const std::vector<StableSetVariable> &LocalSolutionData::variables() const {
   return m_variables;
}
bool LocalSolutionData::isNewSet(const DenseSet &t_set) const {
   return std::all_of(m_variables.begin(), m_variables.end(),
                      [&](const StableSetVariable &variable) {
                         return variable.set() != t_set;
                      });
}

std::size_t LocalSolutionData::addStableSet(pcog::DenseSet t_set) {
   assert(isNewSet(t_set));
   std::size_t index = m_variables.size();
   m_variables.emplace_back(std::move(t_set));
   return index;
}
std::size_t LocalSolutionData::findOrAddStableSet(const DenseSet &t_set) {
   for (std::size_t i = 0; i < m_variables.size(); ++i) {
      if (m_variables[i].set() == t_set) {
         return i;
      }
   }
   return addStableSet(t_set);
}
void LocalSolutionData::addPricingIteration() {
   ++m_pricingIterations;
}
void LocalSolutionData::addLPIterations(std::size_t numIterations) {
   m_lpIterations += numIterations;
}
}