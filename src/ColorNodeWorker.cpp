//
// Created by rolf on 25-2-23.
//

#include "pcog/ColorNodeWorker.hpp"
#include "pcog/Branching.hpp"
#include "pcog/BranchingSelection.hpp"
#include "pcog/ColorSolver.hpp"
#include "pcog/mwss/CombinatorialStableSet.hpp"

// TODO: add error handling, particularly LP-related.
namespace pcog {
void ColorNodeWorker::processNode(BBNode &t_node, ColorSolver &t_solver) {
   // Node preprocessing to compute focus graph and complete focus graph
   setupGraphs(t_node, t_solver);
   // Set up LP representation
   setupLP(t_node, t_solver);
   m_focusNode = t_node.id();
   // Solve lp
   m_lpSolver.solve();
   // Farkas pricing: price in new columns until lp solution is feasible
   while (m_lpSolver.status() == LPSolverStatus::INFEASIBLE) {
      farkasPricing(t_node, t_solver);
   }
   // For graph coloring with ryan-foster branching the subproblem cannot be
   // infeasible, so we can ensure optimality of the lp here.
   if (m_lpSolver.status() != LPSolverStatus::OPTIMAL) {
      // TODO: error handling
      return;
   }

   // Price-lp loop
   pricingLoop(t_node, t_solver);
   computeBranchingVertices(t_node, t_solver);
}
void ColorNodeWorker::setupGraphs(BBNode &node, const ColorSolver &solver) {

   if (node.branchDecisions().size() == 1) {
      // For the root node; simply copy
      m_completeFocusGraph = solver.preprocessedGraph();
      m_focusGraph = solver.preprocessedGraph();
      m_mapToPreprocessed =
          PreprocessedMap({}, NodeMap::identity(m_focusGraph.numNodes()),
                          NodeMap::identity(m_focusGraph.numNodes()));
   } else {
      // Perform node preprocessing
      auto [completeGraph, result] =
          constructBranchedGraphs(solver.preprocessedGraph(),
                                  node.branchDecisions(), node.lowerBound());
      // TODO: check if moving is necesary here or if auto-detected.
      m_completeFocusGraph = completeGraph;
      m_focusGraph = result.graph;
      m_mapToPreprocessed = result.map;
   }
}
void ColorNodeWorker::setupLP(BBNode &bb_node, const ColorSolver &t_solver) {
   // TODO: adjust method below so that the lp is not reinitialized from scratch
   // every time
   m_lpSolver.clear();
   m_lpSolver.setObjectiveSense(ObjectiveSense::MINIMIZE);
   for (Node vertex = 0; vertex < m_mapToPreprocessed.newToOldIDs.size();
        ++vertex) {
      m_lpSolver.addRow({}, 1.0, std::nullopt);
   }
   assert(m_lpSolver.numRows() == m_focusGraph.numNodes());
   std::vector<std::vector<ColElem>> entries;
   std::vector<double> objectives;
   std::vector<double> lbs;
   std::vector<double> ubs;
   for (const auto &variable : t_solver.variables()) {
      std::vector<ColElem> columnEntries;
      for (const auto &node : variable.set()) {
         auto index = m_mapToPreprocessed
                          .oldToNewIDs[node]; // TODO: check if there is not a
                                              // 'double mapping' here
         if (index != INVALID_NODE) {
            columnEntries.push_back(
                ColElem(index, 1.0)); // TODO: fix conversion to int
         }
      }

      // Fix variables to zero from ryan-foster branching
      double upperBound = soplex::infinity; // TODO: experiment with finite upper bound
      const auto &set = variable.set();
      for (const auto &decision : bb_node.branchDecisions()) {
         // Check if the stable set is viable in the Subproblem
         if ((decision.type == BranchType::DIFFER &&
              set.contains(decision.first) && set.contains(decision.second)) ||
             (decision.type == BranchType::SAME &&
              (set.contains(decision.first) !=
               set.contains(decision.second)))) {
            upperBound = 0.0;
            break;
         }
      }
      entries.push_back(columnEntries);
      objectives.push_back(1.0);
      lbs.push_back(0.0);
      ubs.push_back(upperBound);
   }
   m_lpSolver.addColumns(entries,objectives,lbs,ubs);

   //set up LP Basis:
   if(bb_node.depth() != 0){
      LPBasis basis = bb_node.basis();
      NodeMap previousMapping = bb_node.previousNodeMap();
      fixLPBasis(basis,previousMapping);
      m_lpSolver.setBasis(basis);
   }

}
void ColorNodeWorker::farkasPricing(BBNode &t_node, ColorSolver &t_solver) {
   assert(m_lpSolver.status() == LPSolverStatus::INFEASIBLE);
   const auto &stable_sets = t_solver.variables();
   DenseSet colored_nodes(m_completeFocusGraph.numNodes());
   // color any nodes whom's constraints are not enforced in this subproblem
   // (e.g. have a SAME constraint)

   // TODO: is this equivalent to the SAME node thing below? Maybe need to
   // explicitly check LP.
   for (Node node = 0; node < m_mapToPreprocessed.oldToNewIDs.size(); ++node) {
      if (m_mapToPreprocessed.oldToNewIDs[node] == INVALID_NODE) {
         colored_nodes.add(node);
      }
   }

   // Infeasible in the root node, so we are just computing an initial covering
   // basically
   if (t_node.depth() == 0) {
      for (const auto &stable_set : stable_sets) {
         colored_nodes.inplaceUnion(stable_set.set());
      }
   } else {
      RowVector upperBounds = m_lpSolver.columnUpperBounds();
      for (const auto &bound : upperBounds) {
         if (bound.value != 0.0) {
            colored_nodes.inplaceUnion(
                stable_sets[bound.column].set()); // TODO: integer...
         }
      }
   }
   DenseSet uncolored_nodes = colored_nodes;
   uncolored_nodes.complement();

   DenseSet reduced_uncolored_nodes(m_focusGraph.numNodes());
   m_mapToPreprocessed.oldToNewIDs.transform(uncolored_nodes,
                                             reduced_uncolored_nodes);

   PartialUniformWeightFunction weightFunction(reduced_uncolored_nodes);
   MaxWeightStableSetCombinatorial<PartialUniformWeightFunction> mwss_algorithm(
       weightFunction, m_focusGraph);

   DenseSet computedSet(m_focusGraph.numNodes());

   auto lambda = [](const DenseSet &current_nodes, std::size_t size,
                    void *user_data, bool first_solution, bool &stop_solving,
                    bool &accepted_solution) {
      auto *set = reinterpret_cast<DenseSet *>(user_data);
      *set = current_nodes;
      accepted_solution = true;
      stop_solving = false;
   };
   std::size_t limit = std::max(1'000ul, m_focusGraph.numNodes());
   mwss_algorithm.setNodeLimit(limit); // run with a small limit to emulate
                                       // greedily finding a 'good' set.
   mwss_algorithm.setInfiniteUpperBound();
   mwss_algorithm.setTimeLimitCheckInterval(10'000);
   mwss_algorithm.setUserData(&computedSet);
   mwss_algorithm.setCallback(lambda);

   std::vector<DenseSet> pricedSets;
   while (!reduced_uncolored_nodes.empty()) {
      mwss_algorithm.setLowerBound(0);
      mwss_algorithm.run();

      reduced_uncolored_nodes.inplaceDifference(computedSet);
      mwss_algorithm.updateWeightFunction(
          PartialUniformWeightFunction(reduced_uncolored_nodes));

      // map set back to original graph and maximize it in the
      // completeFocusGraph
      DenseSet computedStableSet(m_completeFocusGraph.numNodes());
      m_mapToPreprocessed.newToOldIDs.transform(computedSet, computedStableSet);
      maximizeStableSet(computedStableSet, m_completeFocusGraph);
      assert(m_completeFocusGraph.setIsStable(computedStableSet));
      assert(m_completeFocusGraph.setIsStableMaximal(computedStableSet));
      assert(t_solver.isNewSet(computedStableSet));
      assert(t_node.verifyStableSet(computedStableSet));
      pricedSets.push_back(computedStableSet);
   }
   // Add new columns and resolve LP
   addColumns(pricedSets, t_solver);
   m_lpSolver.solve();
}
void ColorNodeWorker::pricingLoop(BBNode &node, ColorSolver &t_solver) {
   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   while (true) {
      // Run rounding heuristics / integrality check on the current LP solution
      roundingHeuristic(node,t_solver);
      // Solve pricing problem. If new columns can be found, add them and
      // resolve the LP. Otherwise, break away from loop as we have optimally
      // solved this node, and the lower bounds are updated.
      auto result = priceColumn(node, t_solver);
      if (result != PricingResult::FOUND_COLUMN) {
         // TODO: also account for other sources of pricing stopping, such as
         // time limits being hit.
         break;
      }
      // Resolve lp with newly added columns
      m_lpSolver.solve();
      if (m_lpSolver.status() != LPSolverStatus::OPTIMAL) {
         // TODO: error handling
         break;
      }
   }
}
void ColorNodeWorker::computeBranchingVertices(BBNode &node,
                                               const ColorSolver &t_solver) {
   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   if(node.lowerBound() >= t_solver.globalUpperBound()){
      node.setStatus(BBNodeStatus::CUT_OFF);
      return;
   }

   auto primalSol = m_lpSolver.getPrimalSolution();
   std::erase_if(primalSol,[](const RowElem& elem){
      return elem.value == 0.0;
   });
   //Get the candidate branching edges
   auto scoredEdges = getAllBranchingEdgesViolatedInBoth(
       m_focusGraph, primalSol, t_solver.variables(),
       m_mapToPreprocessed.oldToNewIDs);
   std::mt19937 random_device(42); //TODO: make member

   //Score them according to some function
   scoreBranchingCandidates(
       scoredEdges, BranchingStrategy::INTERSECTION_UNION_SIZE, m_focusGraph,
       m_lpSolver,t_solver.variables(),m_mapToPreprocessed.newToOldIDs,m_completeFocusGraph.numNodes());
   std::shuffle(scoredEdges.begin(),scoredEdges.end(),random_device);
   std::stable_sort(scoredEdges.begin(),scoredEdges.end(),[](const ScoredEdge& a, const ScoredEdge& b){return a.score > b.score;});

   //Select pair, checking if they are violated in both branches
   const NodeMap& focusToPreprocessed = m_mapToPreprocessed.newToOldIDs;
   auto best_it = selectBestPair(scoredEdges,SelectionStrategy::VIOLATED_IN_BOTH,
                  primalSol,focusToPreprocessed,t_solver.variables());
   assert(best_it != scoredEdges.end());
   assert(best_it->node1 != INVALID_NODE && best_it->node2 != INVALID_NODE);
   //TODO: how to sort/choose between node1 /node2?
   node.setBranchingNodes(focusToPreprocessed[best_it->node1],focusToPreprocessed[best_it->node2]);
   node.setStatus(BBNodeStatus::BRANCHED);
}
void ColorNodeWorker::maximizeStableSet(DenseSet &t_set,
                                        const DenseGraph &t_graph) {
   assert(t_graph.setIsStable(t_set));
   m_maximizer.maximizeRandomly(t_set, t_graph);
   assert(t_graph.setIsStable(t_set));
   assert(t_graph.setIsStableMaximal(t_set));
}
void ColorNodeWorker::addColumns(const std::vector<DenseSet> &sets,
                                 ColorSolver &t_solver) {
   assert(t_solver.variables().size() == m_lpSolver.numCols());

   // First add the columns to our internal storage
   for (const auto &set : sets) {
      assert(t_solver.isNewSet(set));
      assert(set.capacity() == m_completeFocusGraph.numNodes());
      t_solver.addStableSet(set);
   }
   // Add the columns to the LP
   const NodeMap& mapToFocussed = m_mapToPreprocessed.oldToNewIDs;
   for (const auto &set : sets) {
      ColVector colRows;
      for (const Node elem : set) {
         Node node = mapToFocussed[elem];
         if(node != INVALID_NODE){
            assert(node < m_lpSolver.numRows());
            colRows.emplace_back(node, 1.0);
         }
      }
      // TODO: add all columns at once in soplex (more efficient)
      m_lpSolver.addColumn(
          colRows, 1.0, 0.0,
          soplex::infinity); // TODO: infinity/weak upper bound here
   }
   assert(t_solver.variables().size() == m_lpSolver.numCols());
}
PricingResult ColorNodeWorker::priceColumn(BBNode &node,
                                           ColorSolver &t_solver) {

   // Get dual weights and make them numerically safe
   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   assert(m_lpSolver.numRows() == m_focusGraph.numNodes());
   auto dual_values = m_lpSolver.getDualSolution();
   std::vector<double> dualSolution(m_focusGraph.numNodes(), 0.0);
   for (const auto &row : dual_values) {
      dualSolution[row.column] = row.value;
   }
   double dualSum =
       std::accumulate(dualSolution.begin(), dualSolution.end(), 0.0);
   double test = fabs(dualSum - m_lpSolver.objective());
   assert(
       test <=
       1e-8); // TODO: check why the dual solution is not equal... how? How to
              // query dual variable bounds?
              //  Only breaks when variable upper bounds are 1.0 for DSJC125.9
   // round up weights to prevent negative dual_values
   for (auto &dual : dualSolution) {
      dual = std::max(dual, 0.0);
   }
   SafeDualWeights dualWeights(dualSolution);

   // Find column(s) with negative reduced cost in focusGraph
   bool solvedOptimally = false;
   assert(m_pricedVariables.empty());
   {
      // run a combinatorial algorithm
      // TODO: refactor a bit.
      auto lambda = [](const DenseSet &current_nodes, SafeWeight weight,
                       void *user_data, bool first_solution, bool &stop_solving,
                       bool &accepted_solution) {
         auto *worker = reinterpret_cast<ColorNodeWorker *>(user_data);
         worker->solutionCallback(current_nodes, weight, user_data,
                                  first_solution, stop_solving,
                                  accepted_solution);
      };

      m_mwssSolver.reset();
      m_mwssSolver =
          std::make_unique<MaxWeightStableSetCombinatorial<SafeDualWeights>>(
              dualWeights, m_focusGraph);
      m_mwssSolver->setLowerBound(dualWeights.getOne());
      m_mwssSolver->setInfiniteUpperBound();
      m_mwssSolver->setTimeLimitCheckInterval(10'000);
      m_mwssSolver->setTimeLimit(std::chrono::seconds(3'600)); // TODO:
      m_mwssSolver->setUserData(this);
      m_mwssSolver->setCallback(lambda);
      m_colorSolver = &t_solver;

      m_mwssSolver->run();
      solvedOptimally = !m_mwssSolver->stoppedBeforeOptimality();
   }
   // Add column(s) to the LP and resolve

   // Add variables which were priced
   SafeDualWeights::weight_type bestPricingValue = 0;
   const auto &map = m_mapToPreprocessed.oldToNewIDs;
   std::vector<DenseSet> newSets;
   for (auto &set : m_pricedVariables) {
      if (t_solver.isNewSet(set)) {
         DenseSet focusSet(m_focusGraph.numNodes());
         map.transform(set, focusSet);
         SafeDualWeights::weight_type maxPricingValue =
             dualWeights.setCostUpperBound(focusSet);
         assert(maxPricingValue > dualWeights.getOne());
         if (maxPricingValue > bestPricingValue) {
            bestPricingValue = maxPricingValue;
         }
         // Don't care too much about numerical problems here, as the effect of
         // numerical errors
         //  on the choice of new column should be very minimal in terms of
         //  performance
         double approximatePricingValue =
             double(maxPricingValue) / double(dualWeights.getOne());
         assert(m_completeFocusGraph.setIsStable(set));
         assert(m_completeFocusGraph.setIsStableMaximal(set));
         newSets.push_back(set);
      }
   }
   m_pricedVariables.clear();
   if (solvedOptimally) {
      if (newSets.empty()) {

         SafeDualWeights::weight_type weightSum = dualWeights.fullCost();
         RoundingMode mode = fegetround();
         fesetround(FE_DOWNWARD);
         double safeLowerBound = weightSum;
         safeLowerBound /= dualWeights.scale();
         fesetround(mode);
         assert(m_lpSolver.objective() >= safeLowerBound);
         node.updateLowerBound(safeLowerBound);
         return PricingResult::NO_COLUMNS;
      }
   }
   RoundingMode mode = fegetround();
   fesetround(FE_UPWARD);
   double pricingProblemUB =
       static_cast<double>(bestPricingValue) / dualWeights.getOne();
   assert(pricingProblemUB >= 1.0 - 1e-8); // TODO: make epsilon global?
   fesetround(FE_DOWNWARD);
   double derivedLowerBound = m_lpSolver.objective() / pricingProblemUB;
   fesetround(mode);
   assert(m_lpSolver.objective() >= derivedLowerBound);

   node.updateLowerBound(derivedLowerBound);

   addColumns(newSets, t_solver);
   return PricingResult::FOUND_COLUMN; // TODO: pick up on abort here
}
void ColorNodeWorker::solutionCallback(const DenseSet &current_nodes,
                                       SafeWeight weight, void *user_data,
                                       bool first_solution, bool &stop_solving,
                                       bool &accepted_solution) {

   accepted_solution = false;
   stop_solving = false;

   assert(m_focusGraph.setIsStable(current_nodes));
   const auto &map = m_mapToPreprocessed.newToOldIDs;
   DenseSet original_space_set(m_completeFocusGraph.numNodes());
   map.transform(current_nodes, original_space_set);
   assert(m_completeFocusGraph.setIsStable(original_space_set));
   maximizeStableSet(original_space_set, m_completeFocusGraph);
   assert(m_completeFocusGraph.setIsStable(original_space_set));
   assert(m_completeFocusGraph.setIsStableMaximal(original_space_set));

   bool is_new_set = m_colorSolver->isNewSet(original_space_set);
   if (!is_new_set) {
      return;
   }
   accepted_solution = true;
   m_pricedVariables.emplace_back(original_space_set);
   std::size_t search_extra =
       std::min(10 * m_completeFocusGraph.numNodes(), 1000ul);
   std::size_t numBBnodes = m_mwssSolver->numBranchAndBoundNodes();
   m_mwssSolver->setNodeLimit(numBBnodes + search_extra);
   assert(m_mwssSolver->getWeightFunction().setCostUpperBound(current_nodes) >
          m_mwssSolver->getWeightFunction().getOne());
}
void ColorNodeWorker::roundingHeuristic(BBNode &node, ColorSolver &t_solver) {
   //TODO: fix
   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   double cutOffBound = static_cast<double>(t_solver.globalUpperBound())-1.0;
   if(m_lpSolver.objective() >= (cutOffBound +1e-8) ){
      return;
   }
   const auto& focusGraph = m_focusGraph;
   const auto& preprocessedToFocus = m_mapToPreprocessed.oldToNewIDs;
   const auto& preprocessedGraph = t_solver.preprocessedGraph();
   const auto& variables = t_solver.variables();
   auto LPsol = m_lpSolver.getPrimalSolution();
   assert(!LPsol.empty());

   struct PartialSolVar{
      PartialSolVar(std::size_t n) : index{INVALID_NODE},value{0.0},projectedSet(n){};
      std::size_t index;
      double value;
      DenseSet projectedSet;
   };
   std::vector<PartialSolVar> projectedSolution(LPsol.size(),PartialSolVar(focusGraph.numNodes()));
   for(std::size_t i = 0; i < LPsol.size(); i++){
      projectedSolution[i].value = LPsol[i].value;
      projectedSolution[i].index = LPsol[i].column;
      preprocessedToFocus.transform(variables[LPsol[i].column].set(),projectedSolution[i].projectedSet);
      if(!focusGraph.setIsStable(projectedSolution[i].projectedSet)){
         if(!preprocessedGraph.setIsStable(variables[projectedSolution[i].index].set())){
            std::cout<<"set not stable in original??"<<std::endl;
         }
      }

   }
   DenseSet uncoloredNodes(focusGraph.numNodes(),true);

   std::vector<std::size_t> color_indices;
   while(!uncoloredNodes.empty()){
      auto best = std::max_element(projectedSolution.begin(),projectedSolution.end(),[&](const PartialSolVar& a, const PartialSolVar& b){
         return a.value*(a.projectedSet.intersection(uncoloredNodes)).size() < b.value*(b.projectedSet.intersection(uncoloredNodes)).size();
      });
      uncoloredNodes.inplaceDifference(best->projectedSet);
      color_indices.push_back(best->index);
      if(color_indices.size() > LPsol.size()){
         break;
      }
   }
   if(uncoloredNodes.empty() && color_indices.size() < t_solver.globalUpperBound()){
      //TODO: add setting for storing 'only best' solutions or also suboptimal settings
      t_solver.addSolution(color_indices); //TODO: create solution pool
   }

}
LPBasis ColorNodeWorker::basis() {
   return m_lpSolver.getLPBasis();
}
NodeMap ColorNodeWorker::mapToFocus() const {
   return m_mapToPreprocessed.oldToNewIDs;
}
void ColorNodeWorker::fixLPBasis(LPBasis &basis,
                                 const NodeMap &previous_nodemap) {
   //fix rows of the basis
   //Rows which were not in previous rowmaps are also not in this one
   std::vector<soplex::SPxSolverBase<double>::VarStatus> originalNodeStatus(previous_nodemap.size(),soplex::SPxSolverBase<double>::ON_LOWER);
   for(std::size_t i = 0; i < previous_nodemap.size();++i){
      if(previous_nodemap[i] != INVALID_NODE){
         originalNodeStatus[i] = basis.rowStatus[previous_nodemap[i]];
      }
   }
   assert(m_lpSolver.numRows() == m_mapToPreprocessed.newToOldIDs.size());
   basis.rowStatus.resize(m_lpSolver.numRows());
   for(std::size_t i = 0; i < m_mapToPreprocessed.newToOldIDs.size(); ++i){
      basis.rowStatus[i] = originalNodeStatus[m_mapToPreprocessed.newToOldIDs[i]];
   }

   //fix (extra) columns

   std::size_t previousSize = basis.colStatus.size();
   std::size_t nextSize = m_lpSolver.numCols();
   if(previousSize < nextSize){
      basis.colStatus.resize(nextSize);
      for(std::size_t i = previousSize; i < nextSize; ++i){
         basis.colStatus[i] = soplex::SPxSolverBase<double>::ON_LOWER;
      }
   }

   assert(basis.rowStatus.size() == m_lpSolver.numRows() && basis.colStatus.size() == m_lpSolver.numCols());

}

} // namespace pcog
