//
// Created by rolf on 25-2-23.
//

#include "pcog/ColorNodeWorker.hpp"
#include "pcog/Branching.hpp"
#include "pcog/BranchingSelection.hpp"
#include "pcog/ColorSolver.hpp"
#include "pcog/mwss/CombinatorialStableSet.hpp"
#include "pcog/mwss/AugmentingSearch.hpp"

// TODO: add error handling, particularly LP-related.
namespace pcog {
void ColorNodeWorker::processNode(BBNode &t_node, SolutionData &t_solData) {
   resetNodeStatistics();

   // Node preprocessing to compute focus graph and complete focus graph
   setupGraphs(t_node, t_solData);
   // Set up LP representation
   setupLP(t_node, t_solData);
   m_focusNode = t_node.id();
   // Solve lp
   solveLP();
   // Farkas pricing: price in new columns until lp solution is feasible
   while (m_lpSolver.status() == LPSolverStatus::INFEASIBLE) {
      farkasPricing(t_node, t_solData);
   }
   // For graph coloring with ryan-foster branching the subproblem cannot be
   // infeasible, so we can ensure optimality of the lp here.
   if (m_lpSolver.status() != LPSolverStatus::OPTIMAL) {
      // TODO: error handling
      return;
   }

   // Price-lp loop
   //TODO: pricing stores the basis, but ideally we only want to compute branching vertices after diving.
   //TODO: however, diving destroys LP information. Two options:
   // 1. reset LP manually before computing branching solutions (downside = restoring basis might be expensive)
   // 2. perform diving in a separate LP solver object (downside = lots of copying)
   pricingLoop(t_node, t_solData,false);
   computeBranchingVertices(t_node, t_solData);
   divingHeuristic(t_node,t_solData);

   writeNodeStatistics(t_node,t_solData);
}
void ColorNodeWorker::setupGraphs(BBNode &node, const SolutionData &solver) {

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
      m_completeFocusGraph = completeGraph;
      m_focusGraph = result.graph;
      m_mapToPreprocessed = result.map;
   }
}
void ColorNodeWorker::setupLP(BBNode &bb_node, const SolutionData &t_solData) {
   if(bb_node.depth() != 0 && bb_node.parent() == m_focusNode){
      //TODO: check if this actually saves any time at all

      if(m_focusGraph.numNodes() != m_lpSolver.numRows()){
         assert(m_focusGraph.numNodes() < m_lpSolver.numRows());

         std::vector<int> permutation;
         m_lpSolver.removeRows(permutation);
         //Save permutation of rows

      }
      //TODO: Check if row removal 'preserves' basis

      auto decisions = bb_node.branchDecisions();
      std::size_t numDecisions = decisions.size();
      std::size_t numAdded = bb_node.getNumAddedBranchingDecisions();

      auto bounds = m_lpSolver.columnUpperBounds();
      const auto& variables = t_solData.variables();
      assert(bounds.size() == variables.size());
      //Fix variables to zero for columns which do not fit the current subproblem
      for(std::size_t j = 0; j < bounds.size(); ++j){
         if(bounds[j].value == 0.0) continue;
         const auto &set = variables[j].set();
         for(std::size_t i = numDecisions-numAdded; i < numDecisions; ++i){
            const auto& decision = decisions[i];
            if ((decision.type == BranchType::DIFFER &&
                 set.contains(decision.first) && set.contains(decision.second)) ||
                (decision.type == BranchType::SAME &&
                 (set.contains(decision.first) !=
                  set.contains(decision.second)))) {
               //set upper bound to zero. We do this one at a time (and not all at once), so that SoPlex might be better able to preserve basis information
               m_lpSolver.changeBounds(j,0.0,0.0);
               break;
            }
         }
      }
      return;
   }
   // TODO: adjust method below so that the lp is not reinitialized from scratch
   // every time
   m_lpSolver.clear();

   m_lpSolver.setIntegralityPolishing(true);
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
   for (const auto &variable : t_solData.variables()) {
      std::vector<ColElem> columnEntries;
      for (const auto &node : variable.set()) {
         auto index = m_mapToPreprocessed
                          .oldToNewIDs[node]; // TODO: check if there is not a
                                              // 'double mapping' here
         if (index != INVALID_NODE) {
            columnEntries.push_back(
                ColElem(index, 1.0));
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
void ColorNodeWorker::farkasPricing(BBNode &t_node, SolutionData &t_solData) {
   assert(m_lpSolver.status() == LPSolverStatus::INFEASIBLE);
   const auto &stable_sets = t_solData.variables();
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
      assert(t_solData.isNewSet(computedStableSet));
      assert(t_node.verifyStableSet(computedStableSet));
      pricedSets.push_back(computedStableSet);
   }
   // Add new columns and resolve LP
   addColumns(pricedSets, t_solData);
   solveLP();
}
void ColorNodeWorker::pricingLoop(BBNode &node, SolutionData &t_solData, bool duringDiving) {
   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   while (true) {
      // Run rounding heuristics / integrality check on the current LP solution
      roundingHeuristic(node,t_solData);
      // Solve pricing problem. If new columns can be found, add them and
      // resolve the LP. Otherwise, break away from loop as we have optimally
      // solved this node, and the lower bounds are updated.
      auto result = priceColumn(node, t_solData,duringDiving);
      if (result != PricingResult::FOUND_COLUMN) {
         // TODO: also account for other sources of pricing stopping, such as
         // time limits being hit.
         break;
      }
      // Resolve lp with newly added columns
      solveLP();
      if (m_lpSolver.status() != LPSolverStatus::OPTIMAL) {
         // TODO: error handling
         break;
      }
   }
   node.setBasis(m_lpSolver.getLPBasis()); //TODO: only save here when diving, maybe?
}
void ColorNodeWorker::computeBranchingVertices(BBNode &node,
                                               const SolutionData &t_solData) {
   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   if(node.lowerBound() >= t_solData.upperBound()){
      node.setStatus(BBNodeStatus::CUT_OFF);
      return;
   }

   auto primalSol = m_lpSolver.getPrimalSolution();
   //Get the candidate branching edges
   auto scoredEdges = getAllBranchingEdgesViolatedInBoth(
       m_focusGraph, primalSol, t_solData.variables(),
       m_mapToPreprocessed.oldToNewIDs);
   std::mt19937 random_device(42); //TODO: make member

   //Score them according to some function
   scoreBranchingCandidates(
       scoredEdges, t_solData.settings().branchingStrategy(), m_focusGraph,
       m_lpSolver,t_solData.variables(),m_mapToPreprocessed.newToOldIDs,m_completeFocusGraph.numNodes());
   std::shuffle(scoredEdges.begin(),scoredEdges.end(),random_device);
   std::stable_sort(scoredEdges.begin(),scoredEdges.end(),[](const ScoredEdge& a, const ScoredEdge& b){return a.score > b.score;});

   //Select pair, checking if they are violated in both branches
   const NodeMap& focusToPreprocessed = m_mapToPreprocessed.newToOldIDs;
   auto best_it = selectBestPair(scoredEdges,t_solData.settings().candidateSelectionStrategy(),
                  primalSol,focusToPreprocessed,t_solData.variables());
   assert(best_it != scoredEdges.end());
   assert(best_it->node1 != INVALID_NODE && best_it->node2 != INVALID_NODE);
   //TODO: how to sort/choose between node1 /node2?
   //Do some tests on which node to eliminate when branching (e.g. simple heuristic like degree or so should do)

   constexpr bool neighbourHoodDisjunction = false;
   Node firstNode = focusToPreprocessed[best_it->node1];
   Node secondNode = focusToPreprocessed[best_it->node2];
   node.setBranchingNodes(firstNode,secondNode);

   if(!neighbourHoodDisjunction) {
      node.setBranchingData(
          {{BranchData(firstNode, secondNode, BranchType::SAME)},
           {BranchData(firstNode, secondNode, BranchType::DIFFER)}});
   }else{
      DenseSet differenceOne = m_focusGraph.neighbourhood(best_it->node1).difference(m_focusGraph.neighbourhood(best_it->node2));
      DenseSet differenceTwo = m_focusGraph.neighbourhood(best_it->node2).difference(m_focusGraph.neighbourhood(best_it->node1));

      bool firstSmaller = differenceOne.size() <= differenceTwo.size();

      std::vector<std::vector<BranchData>> data{{BranchData(firstNode,secondNode,BranchType::SAME)}};

      const auto& differChoice = firstSmaller ? differenceOne : differenceTwo;
      Node fixedNode = firstSmaller ? secondNode : firstNode;
      for(const Node neighbourOne : differChoice){
         std::vector<BranchData> branching;
         for(const Node diffOne : differChoice){
            if(diffOne == neighbourOne){
               break;
            }

            branching.emplace_back(focusToPreprocessed[diffOne],fixedNode,BranchType::DIFFER);
         }
         branching.emplace_back(focusToPreprocessed[neighbourOne],fixedNode,BranchType::SAME);
         data.emplace_back(branching);
      }

      node.setBranchingData(data);
   }

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
                                 SolutionData &t_solData) {
   assert(t_solData.variables().size() == m_lpSolver.numCols());

   // First add the columns to our internal storage
   for (const auto &set : sets) {
      assert(t_solData.isNewSet(set));
      assert(set.capacity() == m_completeFocusGraph.numNodes());
      t_solData.addStableSet(set);
   }
   addColumnsToLP(sets,t_solData);
}
PricingResult ColorNodeWorker::priceColumn(BBNode &node,
                                           SolutionData &t_solData,
                                           bool duringDiving) {
   ++m_numPricingIterations;
   // Get dual weights and make them numerically safe
   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   assert(m_lpSolver.numRows() == m_focusGraph.numNodes());
   auto dual_values = m_lpSolver.getDualSolution();
   std::vector<double> dualSolution(m_focusGraph.numNodes(), 0.0);
   for (const auto &row : dual_values) {
      dualSolution[row.column] = row.value;
   }
//  if(!duringDiving){
//         double dualSum = std::accumulate(dualSolution.begin(), dualSolution.end(), 0.0);
//         double test = fabs(dualSum - m_lpSolver.objective());
//         assert(test <= 1e-8);
//         // TODO: check why the dual solution is not equal... how? How to
//         // query dual variable bounds?
//         // Only breaks when variable upper bounds are 1.0 for DSJC125.9
//  }

   // round up weights to prevent negative dual_values
   for (auto &dual : dualSolution) {
      dual = std::max(dual, 0.0);
   }
   SafeDualWeights dualWeights(dualSolution);

   // Find column(s) with negative reduced cost in focusGraph
   bool solvedOptimally = false;
   assert(m_pricedVariables.empty());

   //First attempt to do so heuristically
   bool doHeuristics = false;
   if(doHeuristics){
      AugmentingSearch search(dualWeights.weights(),m_focusGraph);
      std::vector<DenseSet> solutions =
          search.repeatedGreedyLocalSearch<AugmentingSearch::GreedyStrategy::ByWeight>(dualWeights.getOne());
      //TODO: seems to be slower to do the following, test.
//      if(solutions.empty()){
//         auto sols = search.
//                     repeatedGreedyLocalSearch<AugmentingSearch::GreedyStrategy::BySurplus>(dualWeights.getOne(),!solutions.empty());
//         solutions.insert(solutions.end(),sols.begin(),sols.end());
//      }
      for(const auto& sol : solutions){
         DenseSet preprocessedSol(m_completeFocusGraph.numNodes());
         m_mapToPreprocessed.newToOldIDs.transform(sol,preprocessedSol);
         maximizeStableSet(preprocessedSol,m_completeFocusGraph);
         bool duplicate = false;
         for(const auto& other_sol : m_pricedVariables){
            if(other_sol == preprocessedSol){
               duplicate = true;
               break;
            }
         }
         if(!duplicate){
            m_pricedVariables.push_back(preprocessedSol);
         }
      }
   }

   //then exactly
   if(m_pricedVariables.empty() && (!duringDiving || (duringDiving && !doHeuristics))){
      // run an exact combinatorial algorithm
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
      m_mwssSolver->setTimeLimit(std::chrono::seconds(3'600)); // TODO: fix
      if(duringDiving){
         m_mwssSolver->setNodeLimit(10*m_focusGraph.numNodes()); //TODO: fix
      }
      m_mwssSolver->setUserData(this);
      m_mwssSolver->setCallback(lambda);
      m_colorSolver = &t_solData;

      m_mwssSolver->run();
      solvedOptimally = !m_mwssSolver->stoppedBeforeOptimality();
   }
   // Add column(s) to the LP and resolve

   // Add variables which were priced
   SafeDualWeights::weight_type bestPricingValue = 0;
   const auto &map = m_mapToPreprocessed.oldToNewIDs;
   std::vector<DenseSet> newSets;
   for (auto &set : m_pricedVariables) {
      if (t_solData.isNewSet(set)) {
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
   if (newSets.empty()) {
      if (!duringDiving && solvedOptimally) {
         SafeDualWeights::weight_type weightSum = dualWeights.fullCost();
         RoundingMode mode = fegetround();
         fesetround(FE_DOWNWARD);
         double safeLowerBound = weightSum;
         safeLowerBound /= dualWeights.scale();
         fesetround(mode);
         assert(m_lpSolver.objective() >= safeLowerBound);
         node.updateLowerBound(safeLowerBound);
      }
      return PricingResult::NO_COLUMNS;
   }

   if(!duringDiving && solvedOptimally){
      RoundingMode mode = fegetround();
      fesetround(FE_UPWARD);
      double pricingProblemUB =
          static_cast<double>(bestPricingValue) / dualWeights.getOne();
      assert(pricingProblemUB >= 1.0 - t_solData.settings().roundingTolerance());
      fesetround(FE_DOWNWARD);
      double derivedLowerBound = m_lpSolver.objective() / pricingProblemUB;
      fesetround(mode);
      assert(m_lpSolver.objective() >= derivedLowerBound);

      node.updateLowerBound(derivedLowerBound);
   }
   //Ensure we keep the solution data up to date if we improve the global lower bound (e.g. we are in the root node).
   //For now, the case where we are not in the root node but improve the global bound is handled by b&b code itself,
   //but here we need to notify the pricing system somehow
   if(node.depth() == 0 && node.lowerBound() > t_solData.lowerBound()){
      writeNodeStatistics(node,t_solData);
   }
   addColumns(newSets, t_solData);
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
void ColorNodeWorker::roundingHeuristic(BBNode &node, SolutionData &t_solData) {

   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   double cutOffBound = static_cast<double>(t_solData.upperBound())-1.0;
   if(m_lpSolver.objective() >= cutOffBound + t_solData.settings().roundingTolerance()){
      return;
   }
   const auto& focusGraph = m_focusGraph;
   const auto& preprocessedToFocus = m_mapToPreprocessed.oldToNewIDs;
   const auto& preprocessedGraph = t_solData.preprocessedGraph();
   const auto& variables = t_solData.variables();
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

   std::vector<std::vector<PartialSolVar>::const_iterator> color_indices;
   auto projectedSolCopy = projectedSolution;
   while(!uncoloredNodes.empty()){
      auto best = std::max_element(projectedSolution.begin(),projectedSolution.end(),[&](const PartialSolVar& a, const PartialSolVar& b){
         return a.value*a.projectedSet.size() < b.value*b.projectedSet.size();
      });
      uncoloredNodes.inplaceDifference(best->projectedSet);
      color_indices.emplace_back(best);
      if(color_indices.size() > LPsol.size()){
         break;
      }
      DenseSet removeSet = best->projectedSet;
      for(auto& projectedVar : projectedSolution){
         projectedVar.projectedSet.inplaceDifference(removeSet);
      }
   }
   if(uncoloredNodes.empty() && color_indices.size() < t_solData.upperBound()){
      //project solution back to original graph; although the sets form a coloring for this graph after branching, we also want to solution
      //to be valid in the root node.
      DenseSet originalUncolored(m_completeFocusGraph.numNodes(),true);
      for(const auto& var : color_indices){
         originalUncolored.inplaceDifference(variables[var->index].set());
      }
      if(originalUncolored.empty()){
         //no need to repair
         std::vector<std::size_t> indices;
         for(const auto& elem : color_indices){
            indices.emplace_back(elem->index);
         }
         //we synchronize data so it looks nice when we display
         writeNodeStatistics(node,t_solData);
         t_solData.addSolution(indices);
      }else{
         // Need to repair the solution.
         // In this case we need to add columns which are also valid in the current node
         // TODO: somehow try minimizing the # of added columns
         SetColoring setColoring;
         for(const auto& var : color_indices){

            long index = std::distance(projectedSolution.cbegin(),var);
            setColoring.addColor(projectedSolCopy[index].projectedSet);
         }
         NodeColoring coloring(m_focusGraph.numNodes(),setColoring);
         NodeColoring newColoring = extendColoring(coloring,m_mapToPreprocessed,m_completeFocusGraph);

         SetColoring fixedColoring(newColoring);
         std::size_t before = t_solData.variables().size();
         std::vector<DenseSet> colors = fixedColoring.colors();
         for(auto& color : colors){
            maximizeStableSet(color,m_completeFocusGraph);
         }
         auto indices = repairColoring(colors,t_solData);

         writeNodeStatistics(node,t_solData);
         t_solData.addSolution(indices);
         std::size_t after = t_solData.variables().size();

         std::cout<<fixedColoring.numColors()<<"-coloring repaired, " <<after -before <<" new stable sets\n";
         //how to do this; repair one at a time?
      };
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
std::vector<std::size_t>
ColorNodeWorker::repairColoring(const std::vector<DenseSet> &coloring,
                                SolutionData &t_solData) {
   assert(t_solData.variables().size() == m_lpSolver.numCols());

   std::size_t oldIndex = t_solData.variables().size();
   std::vector<std::size_t> colorIndices;
   std::vector<DenseSet> toAddSets;
   // First add the columns to our internal storage
   for (const auto &set : coloring) {
      std::size_t index = t_solData.findOrAddStableSet(set);
      colorIndices.push_back(index);
      if(index >= oldIndex){
         toAddSets.push_back(set);
      }
   }
   addColumnsToLP(toAddSets,t_solData);
   solveLP();
   return colorIndices;
}
void ColorNodeWorker::addColumnsToLP(const std::vector<DenseSet> &sets,
                                     SolutionData &t_solData) {
   // Add the columns to the LP
   // Adding all columns at once actually turns out to be slower, as we typically only add a few columns during pricing.
   // Thus, we add one column at a time.
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
      m_lpSolver.addColumn(
          colRows, 1.0, 0.0,
          soplex::infinity); // TODO: upper bound option
   }
   assert(t_solData.variables().size() == m_lpSolver.numCols());
}
void ColorNodeWorker::divingHeuristic(BBNode &t_node, SolutionData &t_solData) {
   const int frequency = t_solData.settings().divingFrequency();
   const int depth = t_node.depth();
   if( frequency == -1 || (frequency == 0 && depth != 0) ||
       (frequency > 0 && depth % frequency != 0) ){
      return;
   }
   //do not dive if we cannot find an incumbent anyways (TODO make setting for finding equivalent solutions for crossover)
   double cutoffLimit = t_solData.upperBound()-1.0 + t_solData.settings().roundingTolerance();
   if(m_lpSolver.objective() > cutoffLimit){
      return;
   }

//   m_lpSolver.setIntegralityPolishing(true);
   //std::cout<<"starting diving!\n";
   const int pricing_frequency = t_solData.settings().divingPricingFrequency();
   bool doColumnGeneration = (pricing_frequency == 0 && depth == 0) ||
                             (pricing_frequency > 0 && depth % pricing_frequency == 0);
   if(!doColumnGeneration){
      m_lpSolver.setObjectiveUpperLimit(cutoffLimit);
   }

   std::vector<ColIdx> fixedUpperBoundVars;
   double alpha = 1.0; //parameter for weighing importance of LP solution
   double beta = 2.0; //parameter for weighing importance of columns
   DenseSet coveredNodes(m_focusGraph.numNodes());
   DenseSet varSet(m_focusGraph.numNodes());
   const auto& vars = t_solData.variables();

   bool solution_found = false;
   while(true){
      //std::cout<<"Diving objective: "<<m_lpSolver.objective()<<"\n";
      auto sol = m_lpSolver.getPrimalSolution();
      auto best_it = sol.end();
      double best_score = - std::numeric_limits<double>::infinity();
      for(auto it = sol.begin(); it != sol.end(); ++it){
         if(t_solData.settings().isFeasZero(it->value)
             || t_solData.settings().isFeasOne(it->value)){
            continue;
         }
         varSet.clear();
         m_mapToPreprocessed.oldToNewIDs.transform(vars[it->column].set(),varSet);
         //extra tolerance term term to ensure that values close to zero favor the sets which cover more nodes
         double score = (std::pow(it->value,alpha)+ t_solData.settings().roundingTolerance()) *
                        std::pow(varSet.difference(coveredNodes).size(),beta);
         if(score > best_score){
            best_score = score;
            best_it = it;
         }
      }
      if(best_it == sol.end()){
         //no fractional points: we think we found an integral solution
         solution_found = true;
         break;
      }

      //Fix the variable from best_it to 1
      varSet.clear();
      m_mapToPreprocessed.oldToNewIDs.transform(vars[best_it->column].set(),varSet);
      coveredNodes.inplaceUnion(varSet);

      m_lpSolver.changeBounds(best_it->column,1.0,soplex::infinity); //TODO: test setting upper bounds?
      fixedUpperBoundVars.push_back(best_it->column);
      //resolve the LP
      bool status = solveLP();
      if(!status){ //exits when the LP solver has an error or the objective is cut off
         break;
      }

      if(doColumnGeneration){
         pricingLoop(t_node,t_solData,true);
         if(m_lpSolver.objective() >= cutoffLimit){
            break;
         }
      }

   }
   if(solution_found){
      std::vector<std::size_t> indices;
      auto solution = m_lpSolver.getPrimalSolution();
      for(const auto& elem : solution){
         assert(t_solData.settings().isFeasOne(elem.value) ||
                t_solData.settings().isFeasZero(elem.value));
         if(t_solData.settings().isFeasOne(elem.value)){
            indices.push_back(elem.column);
         }
      }

      DenseSet originalUncolored(m_completeFocusGraph.numNodes(),true);
      const auto& variables = t_solData.variables();
      for(const auto& index : indices){
         originalUncolored.inplaceDifference(variables[index].set());
      }
      if(originalUncolored.empty()) {
         // no need to repair
         writeNodeStatistics(t_node,t_solData);
         t_solData.addSolution(indices);
      }else{
         std::cout<<"Could not repair diving solution : "<<m_lpSolver.objective()<<"\n";
      }
   }


   //reset the LP to where it was; remove cutoff bound and reset upper bounds which were fixed to zero
   for(const auto& index : fixedUpperBoundVars){
      m_lpSolver.changeBounds(index,0.0,soplex::infinity);
   }
   if(!doColumnGeneration){
      m_lpSolver.setObjectiveUpperLimit(soplex::infinity);
   }
//   m_lpSolver.setIntegralityPolishing(false);
}
void ColorNodeWorker::processNextNode(SolutionData &t_solData) {
   //TODO: somehow refactor to have cleaner interface for updating the global solution data and this needing to be done less often.
   BBNode bb_node = t_solData.popNextNode();

   processNode(bb_node,t_solData);
   // TODO: check for time limit in processNode and return true/false based on hitting time limits or not

   if (bb_node.status() == BBNodeStatus::BRANCHED) {
      t_solData.createChildren(bb_node,*this);
   }
   //Re-check the lower bounds from the trees
   t_solData.updateTreeBounds();
}
void ColorNodeWorker::resetNodeStatistics() {
   m_numLPIterations = 0;
   m_numPricingIterations = 0;
}
void ColorNodeWorker::writeNodeStatistics(BBNode& t_node, SolutionData &t_data) {
   t_data.addLPIterations(m_numLPIterations);
   t_data.addPricingIterations(m_numPricingIterations);
   m_numLPIterations = 0;
   m_numPricingIterations = 0;

   if(t_node.depth() == 0){
      t_data.updateFractionalLowerBound(t_node.fractionalLowerBound());
      t_data.updateLowerBound(t_node.lowerBound());
   }

}
bool ColorNodeWorker::solveLP() {
   //Before every solve, we update the integrality information
   m_lpSolver.markAllColumnsIntegral();

   bool status = m_lpSolver.solve();
   m_numLPIterations += m_lpSolver.numIterations();
   return status;
}

} // namespace pcog
