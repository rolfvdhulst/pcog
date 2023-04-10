//
// Created by rolf on 7-4-23.
//

#include "pcog/SolutionData.hpp"

namespace pcog {

void SolutionData::addSolution(std::vector<std::size_t> t_stable_set_indices) {
#ifndef NDEBUG
   {
      // assert that coloring is indeed a valid coloring
      DenseSet coveredNodes(m_preprocessedGraph.numNodes());
      for (const auto &index : t_stable_set_indices) {
         assert(index < m_variables.size());
         coveredNodes.inplaceUnion(m_variables[index].set());
      }
      assert(coveredNodes.full());
   }
#endif

   std::size_t ub = t_stable_set_indices.size();
   if (ub < m_upperBound) {
      // prune redundant nodes from the tree
      m_tree.pruneUpperBound(ub);
      m_upperBound = ub;
      m_incumbent_index = m_colorings.size();
      // TODO: notify in interface
   }
   m_colorings.emplace_back(std::move(t_stable_set_indices));
}
std::size_t SolutionData::upperBound() const { return m_upperBound; }
std::size_t SolutionData::lowerBound() const { return m_lowerBound; }
double SolutionData::fractionalLowerBound() const {
   return m_fractionalLowerBound;
}
bool SolutionData::isNewSet(const DenseSet &t_set) const {
   return std::all_of(m_variables.begin(), m_variables.end(),
                      [&](const StableSetVariable &variable) {
                         return variable.set() != t_set;
                      });
}

std::size_t SolutionData::addStableSet(pcog::DenseSet t_set) {
   assert(isNewSet(t_set));
   assert(m_preprocessedGraph.setIsStable(t_set));
   std::size_t index = m_variables.size();
   m_variables.emplace_back(std::move(t_set));
   return index;
}
std::size_t SolutionData::findOrAddStableSet(const DenseSet &t_set) {
   for (std::size_t i = 0; i < m_variables.size(); ++i) {
      if (m_variables[i].set() == t_set) {
         return i;
      }
   }
   return addStableSet(t_set);
}
SolutionData::SolutionData(const DenseGraph &t_originalGraph)
    : m_originalGraph{t_originalGraph} {}

const DenseGraph &SolutionData::originalGraph() const {
   return m_originalGraph;
}

} // namespace pcog