//
// Created by rolf on 25-2-23.
//

#include "pcog/ColorNodeWorker.hpp"
#include "pcog/Branching.hpp"
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

      // TODO: experiment with upper bound set to infinity (e.g. the 'virtual'
      // bounds scip uses) Fix variables to zero from ryan-foster branching
      double upperBound = soplex::infinity;
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
      m_lpSolver.addColumn(columnEntries, 1.0, 0.0, upperBound);
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
                                               const ColorSolver &t_solver) {}
void ColorNodeWorker::maximizeStableSet(DenseSet &t_set,
                                        const DenseGraph &t_graph) {
   assert(t_graph.setIsStable(t_set));
   m_maximizer.maximizeRandomly(t_set, t_graph);
   assert(t_graph.setIsStable(t_set));
   assert(t_graph.setIsStableMaximal(t_set));
}
void ColorNodeWorker::addColumns(const std::vector<DenseSet> &sets,
                                 ColorSolver &t_solver) {

   // First add the columns to our internal storage
   for (const auto &set : sets) {
      assert(t_solver.isNewSet(set));
      t_solver.addStableSet(set);
   }
   // Add the columns to the LP
   for (const auto &set : sets) {
      ColVector colRows;
      for (const auto &elem : set) {
         colRows.emplace_back(elem, 1.0); // TODO: fix conversion
      }
      //TODO: add all columns at once in soplex (more efficient)
      m_lpSolver.addColumn(colRows, 1.0, 0.0,
                           soplex::infinity); // TODO: infinity/weak upper bound here
   }
}
PricingResult ColorNodeWorker::priceColumn(BBNode &node,
                                           ColorSolver &t_solver) {

   // Get dual weights and make them numerically safe
   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   assert(m_lpSolver.numRows() == m_focusGraph.numNodes());
   auto dual_values = m_lpSolver.getDualSolution();
   std::vector<double> dualSolution(m_focusGraph.numNodes(),0.0);
   for (const auto &row : dual_values) {
      dualSolution[row.column] = row.value;
   }
   double dualSum = std::accumulate(dualSolution.begin(),dualSolution.end(),0.0);
   double test =fabs(dualSum-m_lpSolver.objective());
   assert(test <=1e-8); //TODO: check why the dual solution is not equal... how? How to query dual variable bounds?
                        // Only breaks when variable upper bounds are 1.0 for DSJC125.9
   //round up weights to prevent negative dual_values
   for(auto& dual : dualSolution){
      dual = std::max(dual,0.0);
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
         auto * worker = reinterpret_cast<ColorNodeWorker*>(user_data);
         worker->solutionCallback(current_nodes,weight,user_data,first_solution,stop_solving,accepted_solution);
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

   //Add variables which were priced
   SafeDualWeights::weight_type bestPricingValue = 0;
   const auto& map =m_mapToPreprocessed.oldToNewIDs;
   std::vector<DenseSet> newSets;
   for (auto &set : m_pricedVariables) {
      if (t_solver.isNewSet(set)) {
         DenseSet focusSet(m_focusGraph.numNodes());
         map.transform(set,focusSet);
         SafeDualWeights::weight_type maxPricingValue = dualWeights.setCostUpperBound(focusSet);
         assert(maxPricingValue > dualWeights.getOne());
         if(maxPricingValue > bestPricingValue){
            bestPricingValue = maxPricingValue;
         }
         //Don't care too much about numerical problems here, as the effect of numerical errors
         // on the choice of new column should be very minimal in terms of performance
         double approximatePricingValue = double(maxPricingValue) / double(dualWeights.getOne());
         assert(m_completeFocusGraph.setIsStable(set));
         assert(m_completeFocusGraph.setIsStableMaximal(set));
         newSets.push_back(set);
      }
   }
   m_pricedVariables.clear();
   if(solvedOptimally){
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
   double pricingProblemUB = static_cast<double>(bestPricingValue) / dualWeights.getOne();
   assert(pricingProblemUB >= 1.0 - 1e-8); //TODO: make epsilon global?
   fesetround(FE_DOWNWARD);
   double derivedLowerBound = m_lpSolver.objective() / pricingProblemUB;
   fesetround(mode);
   assert(m_lpSolver.objective() >= derivedLowerBound);

   node.updateLowerBound(derivedLowerBound);

   addColumns(newSets, t_solver);
   return PricingResult::FOUND_COLUMN; //TODO: pick up on abort here
}
void ColorNodeWorker::solutionCallback(const DenseSet &current_nodes, SafeWeight weight,
                      void *user_data, bool first_solution, bool &stop_solving,
                      bool &accepted_solution){

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
      assert(m_mwssSolver->getWeightFunction().setCostUpperBound(
                 current_nodes) >
             m_mwssSolver->getWeightFunction().getOne());
}

} // namespace pcog
