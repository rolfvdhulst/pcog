//
// Created by rolf on 25-2-23.
//

#ifndef PCOG_SRC_COLORNODEWORKER_HPP
#define PCOG_SRC_COLORNODEWORKER_HPP

#include "LPSolver.hpp"
#include "StableSetMaximizer.hpp"
#include "mwss/CombinatorialStableSet.hpp"
#include "pcog/SafeDualWeights.hpp"
#include "pcog/BBNode.hpp"
#include "utilities/DenseGraph.hpp"
#include "LocalSolutionData.hpp"
#include "Settings.hpp"

namespace pcog {

class SolutionData;

enum class PricingResult { FOUND_COLUMN, NO_COLUMNS, ABORT };

/// This class is responsible for solving the branch-and-bound nodes e.g.
/// performing the pricing loop. It contains the information which is local to
/// each branch-and-bound node, such as the LP and graph representations.
class ColorNodeWorker {
 public:
   explicit ColorNodeWorker(std::size_t t_id)
       : m_worker_id{t_id},
         m_focusNode{INVALID_BB_NODE},
         m_random_engine{t_id},
         m_maximizer{t_id}{};
   std::size_t id() const {return m_worker_id;}
   void runLoop(SolutionData &t_soldata, std::atomic_bool& stop);
   std::size_t successiveChildrenProcessed() const;
   std::size_t pickChildNode(NodeChildSelectionStrategy strategy, const std::vector<std::size_t>& children);

   LPBasis basis();
   NodeMap focusToPreprocessed() const;

 private:
   void processNode(BBNode &node, SolutionData &t_solver);

   /// Performs 'node preprocessing', diminishing the size of the graphs

   void setupNode(BBNode &node, const SolutionData &t_soLData);
   void setupChildLP(BBNode &node, const SolutionData &t_solver,
                     const PreprocessedMap &t_childMap);
   void setupLPFromScratch(BBNode &node);
   void farkasPricing(BBNode &node);
   void pricingLoop(BBNode &node, SolutionData &t_solver, bool duringDiving);
   void roundingHeuristic(BBNode &node, SolutionData &t_solver);
   void divingHeuristic(BBNode &t_node, SolutionData &t_solver);
   void improvementHeuristic(const std::vector<std::size_t> &solVarindices,
                             BBNode &t_node, SolutionData &t_solver);

   void addSolution(const std::vector<std::size_t> &indices, BBNode &t_node,
               SolutionData &t_solData, bool checkImprovements);

   PricingResult priceColumn(BBNode &node, SolutionData &t_solver,
                             bool duringDiving);
   void solutionCallback(const DenseSet &current_nodes, SafeWeight weight,
                         void *user_data, bool first_solution,
                         bool &stop_solving, bool &accepted_solution);
   /// Checks if the node is cut off. If not, then decides the branching
   /// vertices u and v
   void computeBranchingVertices(BBNode &node, SolutionData &t_solver);

   LPBasis fixLPBasis(const SmallBasis &basis, const NodeMap &previous_nodemap);

   void maximizeStableSet(DenseSet &set, const DenseGraph &graph);
   void addColumns(const std::vector<DenseSet> &sets);
   void addColumnsToLP(const std::vector<DenseSet> &sets);
   std::vector<std::size_t> repairColoring(const std::vector<DenseSet> &coloring);
   std::size_t m_worker_id;
   LocalSolutionData m_localData;
   LPSolver m_lpSolver;

   node_id m_focusNode;

   //Local graph information
   DenseGraph m_focusGraph;
   DenseGraph m_completeFocusGraph;
   PreprocessedMap m_mapToPreprocessed;

   NodeMap m_nodeToLPRow;
   NodeMap m_LPRowToNode;

   std::vector<node_id> m_childNodes; // child nodes of current focus node
   std::size_t m_successiveChildNodesProcessed = 0;

   std::mt19937 m_random_engine;

   //pricing information
   std::unique_ptr<MaxWeightStableSetCombinatorial<SafeDualWeights>>
       m_mwssSolver;
   std::vector<DenseSet> m_pricedVariables;
   StableSetMaximizer m_maximizer;


   bool solveLP();

   void synchronizeStastistics(BBNode &node, SolutionData &t_data);
};
} // namespace pcog

#endif // PCOG_SRC_COLORNODEWORKER_HPP
