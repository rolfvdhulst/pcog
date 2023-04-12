//
// Created by rolf on 25-2-23.
//

#ifndef PCOG_SRC_COLORNODEWORKER_HPP
#define PCOG_SRC_COLORNODEWORKER_HPP

#include "LPSolver.hpp"
#include "SolutionData.hpp"
#include "StableSetMaximizer.hpp"
#include "mwss/CombinatorialStableSet.hpp"
#include "pcog/SafeDualWeights.hpp"
#include "utilities/DenseGraph.hpp"

namespace pcog {


enum class PricingResult{
   FOUND_COLUMN,
   NO_COLUMNS,
   ABORT
};

/// This class is responsible for solving the branch-and-bound nodes e.g.
/// performing the pricing loop. It contains the information which is local to
/// each branch-and-bound node, such as the LP and graph representations.
class ColorNodeWorker {
 public:
   ColorNodeWorker() : m_colorSolver{nullptr},m_maximizer{42}, m_focusNode{INVALID_BB_NODE}{};//TODO: add seed to constructor

   void processNextNode(SolutionData& t_solData);
   LPBasis basis();
   NodeMap mapToFocus() const;
 private:
   void processNode(BBNode& node, SolutionData& t_solver);

   /// Performs 'node preprocessing', diminishing the size of the graphs
   void setupGraphs(BBNode& node, const SolutionData& t_solver);
   void setupLP(BBNode& node, const SolutionData& t_solver);
   void farkasPricing(BBNode& node, SolutionData& t_solver);
   void pricingLoop(BBNode& node, SolutionData& t_solver, bool duringDiving);
   void roundingHeuristic(BBNode& node, SolutionData& t_solver);
   void divingHeuristic(BBNode& t_node, SolutionData& t_solver);

   PricingResult priceColumn(BBNode& node, SolutionData& t_solver, bool duringDiving);
   void solutionCallback(const DenseSet &current_nodes, SafeWeight weight,
    void *user_data, bool first_solution, bool &stop_solving, bool &accepted_solution);
   /// Checks if the node is cut off. If not, then decides the branching vertices u and v
   void computeBranchingVertices(BBNode& node, const SolutionData& t_solver);

   void fixLPBasis(LPBasis& basis, const NodeMap& previous_nodemap);

   void maximizeStableSet(DenseSet& set, const DenseGraph& graph);
   void addColumns(const std::vector<DenseSet>& sets, SolutionData& t_solver);
   void addColumnsToLP(const std::vector<DenseSet>& sets, SolutionData& t_solver);
   std::vector<std::size_t> repairColoring(const std::vector<DenseSet>& coloring, SolutionData& t_solver);
   LPSolver m_lpSolver;
   std::unique_ptr<MaxWeightStableSetCombinatorial<SafeDualWeights>> m_mwssSolver;
   const SolutionData * m_colorSolver;
   std::vector<DenseSet> m_pricedVariables;
   StableSetMaximizer m_maximizer;
   node_id m_focusNode;

   DenseGraph m_focusGraph;
   DenseGraph m_completeFocusGraph;
   PreprocessedMap m_mapToPreprocessed;
};
} // namespace pcog

#endif // PCOG_SRC_COLORNODEWORKER_HPP
