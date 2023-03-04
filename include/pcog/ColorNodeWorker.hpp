//
// Created by rolf on 25-2-23.
//

#ifndef PCOG_SRC_COLORNODEWORKER_HPP
#define PCOG_SRC_COLORNODEWORKER_HPP

#include "BBNode.hpp"
#include "DenseGraph.hpp"
#include "LPSolver.hpp"
#include "StableSetMaximizer.hpp"
#include "mwss/CombinatorialStableSet.hpp"
#include "pcog/SafeDualWeights.hpp"

namespace pcog {

class ColorSolver;

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

   void processNode(BBNode& node, ColorSolver& t_solver);
 private:
   /// Performs 'node preprocessing', diminishing the size of the graphs
   void setupGraphs(BBNode& node, const ColorSolver& t_solver);
   void setupLP(BBNode& node, const ColorSolver& t_solver);
   void farkasPricing(BBNode& node, ColorSolver& t_solver);
   void pricingLoop(BBNode& node, ColorSolver& t_solver);
   void roundingHeuristic(BBNode& node, ColorSolver& t_solver);

   PricingResult priceColumn(BBNode& node, ColorSolver& t_solver);
   void solutionCallback(const DenseSet &current_nodes, SafeWeight weight,
    void *user_data, bool first_solution, bool &stop_solving, bool &accepted_solution);
   /// Checks if the node is cut off. If not, then decides the branching vertices u and v
   void computeBranchingVertices(BBNode& node, const ColorSolver& t_solver);


   void maximizeStableSet(DenseSet& set, const DenseGraph& graph);
   void addColumns(const std::vector<DenseSet>& set, ColorSolver& t_solver);
   LPSolver m_lpSolver;
   std::unique_ptr<MaxWeightStableSetCombinatorial<SafeDualWeights>> m_mwssSolver;
   const ColorSolver * m_colorSolver;
   std::vector<DenseSet> m_pricedVariables;
   StableSetMaximizer m_maximizer;
   node_id m_focusNode;

   DenseGraph m_focusGraph;
   DenseGraph m_completeFocusGraph;
   PreprocessedMap m_mapToPreprocessed;
};
} // namespace pcog

#endif // PCOG_SRC_COLORNODEWORKER_HPP
