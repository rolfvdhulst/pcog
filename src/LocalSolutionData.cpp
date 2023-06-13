//
// Created by rolf on 2-6-23.
//

#include "pcog/LocalSolutionData.hpp"

namespace pcog{
LocalSolutionData::LocalSolutionData() :
                                         m_numStartGlobalVariables{0},
                                         m_lastGlobalAddedIndex{0},
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
std::size_t
LocalSolutionData::mapToGlobalIndex(std::size_t localVarIndex) const {
   if(localVarIndex < m_numStartGlobalVariables){
      return localVarIndex;
   }
   return m_variable_mapping[localVarIndex-m_numStartGlobalVariables];
}
void LocalSolutionData::addSolution(
    const std::vector<std::size_t>& t_stable_set_indices) {
   m_solutions.push_back(t_stable_set_indices);
}
GlobalToLocalMapping LocalSolutionData::getGlobalToLocalMapping() const {
   if(m_variable_mapping.empty()){
      return GlobalToLocalMapping{.m_sameUntilIndex = m_numStartGlobalVariables,.m_globalToLocalMapping={}};
   }
   std::size_t max_global = *std::max_element(m_variable_mapping.begin(),m_variable_mapping.end());
   if(max_global < m_numStartGlobalVariables){
      return GlobalToLocalMapping{.m_sameUntilIndex = m_numStartGlobalVariables,.m_globalToLocalMapping={}};
   }
   std::vector<size_t> globalToLocalMap(max_global-m_numStartGlobalVariables+1,std::numeric_limits<std::size_t>::max());
   for(std::size_t i = 0; i < m_variable_mapping.size(); ++i){
      std::size_t globalIndex = m_variable_mapping[i];
      if(globalIndex >= m_numStartGlobalVariables){
         globalToLocalMap[globalIndex-m_numStartGlobalVariables] = m_numStartGlobalVariables + i;
      }
   }

   return GlobalToLocalMapping{.m_sameUntilIndex = m_numStartGlobalVariables,
                               .m_globalToLocalMapping = globalToLocalMap};
}
std::size_t GlobalToLocalMapping::mapGlobalToLocal(std::size_t index) const {
   if(index < m_sameUntilIndex){
      return index;
   }
   if(index >= m_sameUntilIndex + m_globalToLocalMapping.size()){
      return std::numeric_limits<std::size_t>::max();
   }
   return m_globalToLocalMapping[index-m_sameUntilIndex];
}
}