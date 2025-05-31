//
// Created by rolf on 2-6-23.
//

#ifndef PCOG_SRC_LOCALSOLUTIONDATA_HPP
#define PCOG_SRC_LOCALSOLUTIONDATA_HPP

#include "utilities/DenseSet.hpp"

namespace pcog{
class StableSetVariable{
 public:
   StableSetVariable() = default;
   explicit StableSetVariable(DenseSet set) : m_set{std::move(set)}{};
   [[nodiscard]] const DenseSet& set() const { return m_set;}
 private:
   DenseSet m_set;
};

struct LowerBoundInfo{
   LowerBoundInfo(double t_frac, std::size_t t_lb) : m_fractional_lb{t_frac}, m_lb{t_lb}{};
   LowerBoundInfo() : m_fractional_lb{0.0},m_lb{0}{};
   double m_fractional_lb;
   std::size_t m_lb;
};

struct GlobalToLocalMapping{
   std::size_t m_sameUntilIndex;
   std::vector<std::size_t> m_globalToLocalMapping;
   std::size_t mapGlobalToLocal(std::size_t index) const;
};
///Local solution information which is synchronized with the global information.
class LocalSolutionData{
 public:
   friend class SolutionData;
   LocalSolutionData();

   [[nodiscard]] const std::vector<StableSetVariable>& variables() const;
   [[nodiscard]] bool isNewSet(const DenseSet& t_set) const;
   std::size_t addStableSet(DenseSet t_set);
   std::size_t findOrAddStableSet(const DenseSet& t_set);

   void addPricingIteration();
   void addLPIterations(std::size_t numIterations);
   void addSolution(const std::vector<std::size_t>& t_stable_set_indices);
   [[nodiscard]] std::size_t mapToGlobalIndex(std::size_t localVarIndex) const;
   GlobalToLocalMapping getGlobalToLocalMapping() const;

   std::size_t getPricingIterations() const;
 private:

   std::vector<StableSetVariable> m_variables;
   //Until this index, the indexing of global and local variables should be the same.
   //Afterwards, the local variables
   std::size_t m_numStartGlobalVariables;

   std::vector<std::size_t> m_variable_mapping;//indices the mapping of the variables local to global index
                                                // from m_numStartGlobalVariables to m_lastGlobalAddedIndex
   std::size_t m_lastGlobalAddedIndex{}; //the index of the last variable which was added to the global pool

   //Solutions which were found locally
   std::vector<std::vector<std::size_t>> m_solutions;

   std::size_t m_lpIterations;
   std::size_t m_pricingIterations;

   LowerBoundInfo info;
};

}

#endif // PCOG_SRC_LOCALSOLUTIONDATA_HPP
