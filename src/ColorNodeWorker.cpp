//
// Created by rolf on 25-2-23.
//

#include "pcog/ColorNodeWorker.hpp"
#include "pcog/Branching.hpp"
#include "pcog/BranchingSelection.hpp"
#include "pcog/ColorSolver.hpp"
#include "pcog/TabuColoring.hpp"
#include "pcog/mwss/AugmentingSearch.hpp"
#include "pcog/mwss/CombinatorialStableSet.hpp"
#include "pcog/SolutionData.hpp"
#include <thread>

namespace pcog {
void ColorNodeWorker::processNode(BBNode &t_node, SolutionData &t_solData, std::atomic_bool& stop) {
   // Node preprocessing to compute focus graph and complete focus graph
   setupNode(t_node,t_solData,stop);
   if(m_cancelCurrentNode || t_solData.checkTimeLimitHit(stop)){
      cleanUpOnCancel(t_node,t_solData,stop);
      return;
   }
   // Solve lp
   LPSolverStatus status = solveLP(t_solData);
   if(!m_cancelCurrentNode && status == LPSolverStatus::ABORT_TIMEOUT){
      t_solData.stopComputation(stop);
   }else if(status == LPSolverStatus::ERROR){
      std::cout<<"LP Error!\n";
      t_solData.stopComputation(stop);
   }
   if(m_cancelCurrentNode ){
      cleanUpOnCancel(t_node,t_solData,stop);
      return;
   }
   // Farkas pricing: price in new columns until lp solution is feasible
   while (status == LPSolverStatus::INFEASIBLE) {
      farkasPricing(t_node);
      status = solveLP(t_solData);
      if(!m_cancelCurrentNode && status == LPSolverStatus::ABORT_TIMEOUT){
         t_solData.stopComputation(stop);
      }else if(status == LPSolverStatus::ERROR){
         std::cout<<"LP Error!\n";
         t_solData.stopComputation(stop);
      }
      if(m_cancelCurrentNode){
         cleanUpOnCancel(t_node,t_solData,stop);
         return;
      }
   }

   // For graph coloring with ryan-foster branching the subproblem cannot be
   // infeasible, so we can ensure optimality of the lp here.
   assert(status != LPSolverStatus::INFEASIBLE);

   // Price-lp loop
   //TODO: test diving in a separate LP solver
   pricingLoop(t_node, t_solData,false,stop);
   if(m_cancelCurrentNode){ //time limit is checked in pricing loop, which would cancel this node if it was hit
      cleanUpOnCancel(t_node,t_solData,stop);
      return;
   }
   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   auto primalSol = m_lpSolver.getPrimalSolution();
   auto dualSol = m_lpSolver.getDualSolution();

   divingHeuristic(t_node,t_solData,stop);

   if(m_cancelCurrentNode || t_solData.checkTimeLimitHit(stop)){
      cleanUpOnCancel(t_node,t_solData,stop);
      return;
   }

   computeBranchingVertices(t_node, t_solData,primalSol,dualSol);
   if(m_cancelCurrentNode || t_solData.checkTimeLimitHit(stop)){
      //If the node is cancelled after computing branching vertices, this is okay by us
      //In this case, we simply save them
      cleanUpOnCancel(t_node,t_solData,stop);
      return;
   }


   synchronizeStastistics(t_node,t_solData,stop);
}

void ColorNodeWorker::farkasPricing(BBNode &t_node) {
   assert(m_lpSolver.status() == LPSolverStatus::INFEASIBLE);
   const auto &stable_sets = m_localData.variables();
   DenseSet colored_nodes(m_completeFocusGraph.numNodes());
   // color any nodes whom's constraints are not enforced in this subproblem
   // (e.g. have a SAME constraint)

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
                stable_sets[bound.column].set()); // TODO: integer conversions
         }
      }
   }
   DenseSet uncolored_nodes = colored_nodes;
   uncolored_nodes.complement();
   assert(uncolored_nodes.any());

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
   mwss_algorithm.setInterrupt(&m_cancelCurrentNode);

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
      assert(m_localData.isNewSet(computedStableSet));
      assert(t_node.verifyStableSet(computedStableSet));
      pricedSets.push_back(computedStableSet);
   }
   // Add new columns and resolve LP
   addColumns(pricedSets);
}
void ColorNodeWorker::pricingLoop(BBNode &node, SolutionData &t_solData, bool duringDiving,
                                  std::atomic_bool& stop) {
   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   std::size_t numPricingRoundsLimit = std::numeric_limits<std::size_t>::max(); //TODO: limit pricing during diving? Or not, maybe?
   std::size_t numRounds = 0;
   while (numRounds < numPricingRoundsLimit && !m_cancelCurrentNode) {
      // Run rounding heuristics / integrality check on the current LP solution
      roundingHeuristic(node,t_solData,stop);
      if(m_cancelCurrentNode || t_solData.checkTimeLimitHit(stop)){ //check if heuristics found an optimal solution
         break;
      }
      // Solve pricing problem. If new columns can be found, add them and
      // resolve the LP. Otherwise, break away from loop as we have optimally
      // solved this node, and the lower bounds are updated.
      auto result = priceColumn(node, t_solData,duringDiving,stop);
      if(result != PricingResult::FOUND_COLUMN){
         break;
      }
      if (m_cancelCurrentNode || t_solData.checkTimeLimitHit(stop)) {
         return;
      }
      // Resolve lp with newly added columns
      LPSolverStatus status = solveLP(t_solData);
      if(!m_cancelCurrentNode && status == LPSolverStatus::ABORT_TIMEOUT){
         t_solData.stopComputation(stop);
      }else if(status == LPSolverStatus::ERROR){
         std::cout<<"LP Error!\n";
         t_solData.stopComputation(stop);
      }
      if(m_cancelCurrentNode ){
         break;
      }

      ++numRounds;
   }
   if(!duringDiving){
      node.setBasis(toSmallBasis(m_lpSolver.getLPBasis()));
   }
}
void ColorNodeWorker::computeBranchingVertices(BBNode &node,
                                               SolutionData &t_solData,
                                               const RowVector& t_primalSol,
                                               const RowVector& t_dualSol) {
   //TODO: divide function into smaller functions
   if(node.lowerBound() >= t_solData.upperBoundUnscaled()){
      node.setStatus(BBNodeStatus::CUT_OFF);
      return;
   }

   //Get the candidate branching edges
   auto scoredEdges = getAllBranchingEdgesViolatedInBoth(
       m_focusGraph, t_primalSol, m_localData.variables(),
       m_mapToPreprocessed.oldToNewIDs);

   bool alwaysCheckSmallDifference = false;
   if(alwaysCheckSmallDifference){
      //Score them according to some function
      scoreBranchingCandidates(
          scoredEdges, BranchingStrategy::SMALL_DIFFERENCE, m_focusGraph,
          t_primalSol,t_dualSol, m_nodeToLPRow,
          m_localData.variables(),m_mapToPreprocessed.newToOldIDs,m_completeFocusGraph.numNodes());
      std::shuffle(scoredEdges.begin(),scoredEdges.end(),m_random_engine);
      std::stable_sort(scoredEdges.begin(),scoredEdges.end(),[](const ScoredEdge& a, const ScoredEdge& b){return a.score > b.score;});

      const NodeMap& focusToPreprocessed = m_mapToPreprocessed.newToOldIDs;
      auto best_it = selectBestPair(scoredEdges,t_solData.settings().candidateSelectionStrategy(),
                                    t_primalSol,focusToPreprocessed,m_localData.variables());
      assert(best_it != scoredEdges.end());
      assert(best_it->node1 != INVALID_NODE && best_it->node2 != INVALID_NODE);
      //TODO: how to sort/choose between node1 /node2?
      //Do some tests on which node to eliminate when branching (e.g. simple heuristic like degree or so should do)

      Node firstNode = focusToPreprocessed[best_it->node1];
      Node secondNode = focusToPreprocessed[best_it->node2];

      DenseSet differenceOne = m_focusGraph.neighbourhood(best_it->node1).difference(m_focusGraph.neighbourhood(best_it->node2));
      DenseSet differenceTwo = m_focusGraph.neighbourhood(best_it->node2).difference(m_focusGraph.neighbourhood(best_it->node1));

      std::size_t firstSize = differenceOne.size();
      std::size_t secondSize = differenceTwo.size();
      bool firstSmaller = firstSize <= secondSize;
      std::size_t minSize = firstSmaller ? firstSize : secondSize;
      assert(minSize != 0); //node should have been preprocessed in this case.
      if(minSize == 1){
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

         node.setBranchingNodes(firstNode,secondNode);
         node.setBranchingData(data);
         node.setStatus(BBNodeStatus::BRANCHED);
         return;
      }
   }

   //Score them according to some function
   std::shuffle(scoredEdges.begin(),scoredEdges.end(),m_random_engine);
   scoreBranchingCandidates(
       scoredEdges, t_solData.settings().branchingStrategy(), m_focusGraph,
       t_primalSol,t_dualSol, m_nodeToLPRow,
       m_localData.variables(),m_mapToPreprocessed.newToOldIDs,m_completeFocusGraph.numNodes());
   std::stable_sort(scoredEdges.begin(),scoredEdges.end(),[](const ScoredEdge& a, const ScoredEdge& b){return a.score > b.score;});

   //Select pair, checking if they are violated in both branches
   const NodeMap& focusToPreprocessed = m_mapToPreprocessed.newToOldIDs;
   auto best_it = selectBestPair(scoredEdges,t_solData.settings().candidateSelectionStrategy(),
                  t_primalSol,focusToPreprocessed,m_localData.variables());
   assert(best_it != scoredEdges.end());
   assert(best_it->node1 != INVALID_NODE && best_it->node2 != INVALID_NODE);
   //TODO: how to sort/choose between node1 /node2?
   //Do some tests on which node to eliminate when branching (e.g. simple heuristic like degree or so should do)

   constexpr bool neighbourHoodDisjunction = false;
   Node firstNode = focusToPreprocessed[best_it->node1];
   Node secondNode = focusToPreprocessed[best_it->node2];
   node.setBranchingNodes(firstNode,secondNode);

   if(neighbourHoodDisjunction) {
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
   }else{
      node.setBranchingData(
          {{BranchData(firstNode, secondNode, BranchType::SAME)},
           {BranchData(firstNode, secondNode, BranchType::DIFFER)}});
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
void ColorNodeWorker::addColumns(const std::vector<DenseSet> &sets) {
   assert(m_localData.variables().size() == m_lpSolver.numCols());

   // First add the columns to our internal storage
   for (const auto &set : sets) {
      assert(m_localData.isNewSet(set));
      assert(set.capacity() == m_completeFocusGraph.numNodes());
      m_localData.addStableSet(set);
   }
   addColumnsToLP(sets);
}
PricingResult ColorNodeWorker::priceColumn(BBNode &node,
                                           SolutionData &t_solData,
                                           bool duringDiving,
                                           std::atomic_bool& stop) {
   m_localData.addPricingIteration();
   // Get dual weights and make them numerically safe
   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   assert(m_lpSolver.numRows() == m_focusGraph.numNodes());
   auto dual_values = m_lpSolver.getDualSolution();
   std::vector<double> dualSolution(m_focusGraph.numNodes(), 0.0);
   for (const auto &row : dual_values) {
      dualSolution[m_LPRowToNode[row.column]] = row.value;
   }

   // round up weights to prevent negative dual_values
   for (auto &dual : dualSolution) {
      dual = std::max(dual, 0.0);
   }
   SafeDualWeights dualWeights(dualSolution);

   // Find column(s) with negative reduced cost in focusGraph
   bool solvedOptimally = false;
   assert(m_pricedVariables.empty());

   //First attempt to do so heuristically
   bool doHeuristics = true;
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
   if(m_pricedVariables.empty() && (!duringDiving || !doHeuristics)){
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
      auto timeLeft = t_solData.timeLeft();
      if(timeLeft.has_value()){
         m_mwssSolver->setTimeLimit(timeLeft.value());
      }
      m_mwssSolver->setInterrupt(&m_cancelCurrentNode);
      if(duringDiving){
         m_mwssSolver->setNodeLimit(10*m_focusGraph.numNodes()); //TODO: fix
      }
      m_mwssSolver->setUserData(this);
      m_mwssSolver->setCallback(lambda);

      m_mwssSolver->run();
      solvedOptimally = !m_mwssSolver->stoppedBeforeOptimality();
   }
   // Add column(s) to the LP and resolve

   // Add variables which were priced
   SafeDualWeights::weight_type bestPricingValue = 0;
   const auto &map = m_mapToPreprocessed.oldToNewIDs;
   std::vector<DenseSet> newSets;
   for (auto &set : m_pricedVariables) {
      if (m_localData.isNewSet(set)) {
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
      bool improvedLB = false;
      if (!duringDiving && solvedOptimally) {
         SafeDualWeights::weight_type weightSum = dualWeights.fullCost();
         RoundingMode mode = fegetround();
         fesetround(FE_DOWNWARD);
         double safeLowerBound = weightSum;
         safeLowerBound /= dualWeights.scale();
         if(!m_mapToPreprocessed.fixed_sets.empty()){
            double fixedSetsValue = m_mapToPreprocessed.fixed_sets.size(); //TODO: does downward rounding apply here?
            safeLowerBound += fixedSetsValue;
         }
         fesetround(mode);
         assert(m_lpSolver.objective() >= safeLowerBound - m_mapToPreprocessed.fixed_sets.size());
         improvedLB = node.updateLowerBound(safeLowerBound);
      }

      if(node.depth() == 0 && improvedLB){
         synchronizeStastistics(node,t_solData,stop);
      }
      return PricingResult::NO_COLUMNS;
   }
   bool improvedLB = false;
   if(!duringDiving && solvedOptimally){
      RoundingMode mode = fegetround();
      fesetround(FE_UPWARD);
      double pricingProblemUB =
          static_cast<double>(bestPricingValue) / dualWeights.getOne();
      assert(pricingProblemUB >= 1.0 - t_solData.settings().roundingTolerance());
      fesetround(FE_DOWNWARD);
      double derivedLowerBound = m_lpSolver.objective() / pricingProblemUB;
      if(!m_mapToPreprocessed.fixed_sets.empty()){
         double fixedSetsValue = m_mapToPreprocessed.fixed_sets.size(); //TODO: does downward rounding apply here?
         derivedLowerBound += fixedSetsValue;
      }
      fesetround(mode);
      assert(m_lpSolver.objective() >= derivedLowerBound-m_mapToPreprocessed.fixed_sets.size());

      improvedLB = node.updateLowerBound(derivedLowerBound);
   }
   //Ensure we keep the solution data up to date if we improve the global lower bound (e.g. we are in the root node).
   //For now, the case where we are not in the root node but improve the global bound is handled by b&b code itself,
   //but here we need to notify the pricing system somehow
   if(node.depth() == 0 && improvedLB ){
      synchronizeStastistics(node,t_solData,stop);
   }
   addColumns(newSets);
   return PricingResult::FOUND_COLUMN; // TODO: pick up on abort and no found columns here
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

   bool is_new_set = m_localData.isNewSet(original_space_set);
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
void ColorNodeWorker::roundingHeuristic(BBNode &node, SolutionData &t_solData,std::atomic_bool& stop) {

   assert(m_lpSolver.status() == LPSolverStatus::OPTIMAL);
   double cutOffBound = static_cast<double>(t_solData.upperBoundUnscaled())-1.0 - m_mapToPreprocessed.fixed_sets.size();
   if(m_lpSolver.objective() >= cutOffBound + t_solData.settings().roundingTolerance()){
      return;
   }
   const auto& focusGraph = m_focusGraph;
   const auto& preprocessedToFocus = m_mapToPreprocessed.oldToNewIDs;
   const auto& preprocessedGraph = t_solData.preprocessedGraph();
   const auto& variables = m_localData.variables();
   auto LPsol = m_lpSolver.getPrimalSolution();
   assert(!LPsol.empty());

   struct PartialSolVar{
      explicit PartialSolVar(std::size_t n) : index{INVALID_NODE},value{0.0},projectedSet(n){};
      std::size_t index;
      double value;
      DenseSet projectedSet;
   };
   std::vector<PartialSolVar> projectedSolution;
   projectedSolution.reserve(LPsol.size());
   for(const auto& entry :  LPsol){
      PartialSolVar var(focusGraph.numNodes());
      preprocessedToFocus.transform(variables[entry.column].set(),var.projectedSet);
      if(!focusGraph.setIsStable(var.projectedSet)){
         //this happens because variables fixed to zero still enter the LP solution in soplex ('trickle flow')
         assert(fabs(entry.value) < 1e-6); //we filter these out
         continue ;
      }
      var.value = entry.value;
      var.index = entry.column;
      projectedSolution.push_back(var);
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
   bool addSuboptimalSolutions = false;
   if(uncoloredNodes.empty() &&
       (color_indices.size() + m_mapToPreprocessed.fixed_sets.size() < t_solData.upperBoundUnscaled() || addSuboptimalSolutions)){
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
         addSolution(indices,node,t_solData,true,stop);
      }else{
         // Need to repair the solution.
         // In this case we need to add columns which are also valid in the current node
         // TODO: somehow try minimizing the # of added columns
         SetColoring setColoring;
         for(const auto& var : color_indices){
            long index = std::distance(projectedSolution.cbegin(),var);
            setColoring.addColor(projectedSolCopy[index].projectedSet);
         }
         assert(setColoring.isValid(m_focusGraph));
         NodeColoring coloring(m_focusGraph.numNodes(),setColoring);
         NodeColoring newColoring = extendColoring(coloring,m_mapToPreprocessed,m_completeFocusGraph);

         SetColoring fixedColoring(newColoring);
         assert(fixedColoring.isValid(m_completeFocusGraph));
         std::vector<DenseSet> colors = fixedColoring.colors();
         for(auto& color : colors){
            maximizeStableSet(color,m_completeFocusGraph);
         }
         auto indices = repairColoring(colors,t_solData);

         addSolution(indices,node,t_solData,true,stop); //Currently we have already added the 'repaired' columns to the LP; we probably want to only do this if the solution is 'worth saving'
      };
   }

}

NodeMap ColorNodeWorker::focusToPreprocessed() const {
   return m_mapToPreprocessed.newToOldIDs;
}
LPBasis ColorNodeWorker::fixLPBasis(const SmallBasis &basis,
                                 const NodeMap &previous_nodemap) {
   LPBasis largeBasis(m_lpSolver.numRows(),m_lpSolver.numCols()); //Fills status with ON_LOWER initially
   //fix rows of the basis
   assert(m_lpSolver.numRows() == m_mapToPreprocessed.newToOldIDs.size());

   for(const auto& index : basis.basicRows){
      std::size_t preprocessedIndex = previous_nodemap[index];
      std::size_t newIndex = m_mapToPreprocessed.oldToNewIDs[preprocessedIndex];
      if(newIndex != INVALID_NODE){
         largeBasis.rowStatus[newIndex] = soplex::SPxSolverBase<double>::BASIC;
      }
   }

   // Set all fixed columns; This does not seem to actually change the # of iterations,
   // but it stays closer to soplex' basis output so the chance we are doing
   // something wrong is smaller
   auto upperBounds = m_lpSolver.columnUpperBounds();
   for (std::size_t i = 0; i < upperBounds.size(); ++i) {
      if (upperBounds[i].value == 0.0) {
         largeBasis.colStatus[i] = soplex::SPxSolverBase<double>::FIXED;
      }
   }
   for(const auto& index : basis.basicCols){
      largeBasis.colStatus[index] = soplex::SPxSolverBase<double>::BASIC; //TODO: experiment with fixing columns which were zero'd out?
   }

   assert(largeBasis.rowStatus.size() == m_lpSolver.numRows() && largeBasis.colStatus.size() == m_lpSolver.numCols());
   return largeBasis;
}
std::vector<std::size_t>
ColorNodeWorker::repairColoring(const std::vector<DenseSet> &coloring,
                                const SolutionData& t_solData) {
   assert(m_localData.variables().size() == m_lpSolver.numCols());

   std::size_t oldIndex = m_localData.variables().size();
   std::vector<std::size_t> colorIndices;
   std::vector<DenseSet> toAddSets;
   // First add the columns to our internal storage
   for (const auto &set : coloring) {
      std::size_t index = m_localData.findOrAddStableSet(set);
      colorIndices.push_back(index);
      if(index >= oldIndex){
         toAddSets.push_back(set);
      }
   }
   addColumnsToLP(toAddSets);
   //Here we resolve the LP so that the new columns are also taken into account;
   //otherwise, we lose consistency/ cannot query the LP status anymore
   solveLP(t_solData); //TODO: fix resolving error handling etc. here
   return colorIndices;
}
void ColorNodeWorker::addColumnsToLP(const std::vector<DenseSet> &sets) {
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
      double upperBound = m_completeFocusGraph.setIsStable(set) ? soplex::infinity : 0.0;
      m_lpSolver.addColumn(colRows, 1.0, 0.0, upperBound);
   }
   assert(m_localData.variables().size() == m_lpSolver.numCols());
}
void ColorNodeWorker::divingHeuristic(BBNode &t_node, SolutionData &t_solData, std::atomic_bool& stop) {
   const int frequency = t_solData.settings().divingFrequency();
   const int depth = t_node.depth();
   if( frequency == -1 || (frequency == 0 && depth != 0) ||
       (frequency > 0 && depth % frequency != 0) ){
      return;
   }
   //do not dive if we cannot find an incumbent anyways (TODO make setting for finding equivalent solutions for crossover)
   double cutoffLimit = t_solData.upperBoundUnscaled()-1.0 + t_solData.settings().roundingTolerance()
       - m_mapToPreprocessed.fixed_sets.size();
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
   const auto& vars = m_localData.variables();

   bool solution_found = false;
   while(!m_cancelCurrentNode){
//      std::cout<<"Diving objective: "<<m_lpSolver.objective()<<"\n";
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
      LPSolverStatus status = solveLP(t_solData);
      if(!m_cancelCurrentNode && status == LPSolverStatus::ABORT_TIMEOUT){
         t_solData.stopComputation(stop);
      }else if(status == LPSolverStatus::ERROR){
         std::cout<<"LP Error!\n";
         t_solData.stopComputation(stop);
      }else if (status == LPSolverStatus::ABORT_VALUE){
         break;
      }

      if(m_cancelCurrentNode ){ //exits when the LP solver has an error or the objective is cut off
         break;
      }

      if(doColumnGeneration){
         pricingLoop(t_node,t_solData,true,stop);
         if(m_cancelCurrentNode ){
            break;
         }
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
      const auto& variables = m_localData.variables();
      for(const auto& index : indices){
         originalUncolored.inplaceDifference(variables[index].set());
      }
      if(originalUncolored.empty()) {
         // no need to repair
         addSolution(indices,t_node,t_solData,true,stop);
      }else{
         //TODO: fix
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
//   m_lpSolver.setIntegralityPolishing(false); //TODO: check if we need to set/unset this or if this matters at all

}

void ColorNodeWorker::synchronizeStastistics(BBNode& t_node, SolutionData &t_data,std::atomic_bool& stop) {
   t_data.synchronizeLocalDataStatistics(m_localData);
   if(t_node.depth() == 0){
      LowerBoundInfo info(t_node.fractionalLowerBound(), t_node.lowerBound());
      t_data.syncLocalLowerBound(info,m_worker_id,stop);
   }

}
LPSolverStatus ColorNodeWorker::solveLP(const SolutionData &t_data) {
   //Before every solve, we update the integrality information
   m_lpSolver.markAllColumnsIntegral();
   auto timeLeft = t_data.timeLeft();

   m_lpSolver.setTimeLimit(timeLeft);
   LPSolverStatus status = m_lpSolver.solve(&m_cancelCurrentNode);
   m_localData.addLPIterations(m_lpSolver.numIterations());
   return status;
}
void ColorNodeWorker::setupLPFromScratch(BBNode &bb_node) {
   m_lpSolver.clear();

   m_lpSolver.setIntegralityPolishing(true);
   m_lpSolver.setObjectiveSense(ObjectiveSense::MINIMIZE);
   for (Node vertex = 0; vertex < m_mapToPreprocessed.newToOldIDs.size();
        ++vertex) {
      m_lpSolver.addRow({}, 1.0, std::nullopt);
   }
   m_nodeToLPRow = NodeMap::identity(m_focusGraph.numNodes());
   m_LPRowToNode = NodeMap::identity(m_focusGraph.numNodes());

   assert(m_lpSolver.numRows() == m_focusGraph.numNodes());
   std::vector<std::vector<ColElem>> entries;
   std::vector<double> objectives;
   std::vector<double> lbs;
   std::vector<double> ubs;
   for (const auto &variable : m_localData.variables()) {
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
      if(!m_completeFocusGraph.setIsStable(set)){
         upperBound = 0.0;
      }

      entries.push_back(columnEntries);
      objectives.push_back(1.0);
      lbs.push_back(0.0);
      ubs.push_back(upperBound);
   }
   m_lpSolver.addColumns(entries,objectives,lbs,ubs);

   //set up LP Basis:
   if(bb_node.depth() != 0){
      LPBasis basis = fixLPBasis(bb_node.basis(),bb_node.previousNodeMap());
      m_lpSolver.setBasis(basis);
   }
}
void ColorNodeWorker::setupNode(BBNode &node, SolutionData &solver, std::atomic_bool& stop) {
   bool isChildNode = node.parent() == m_focusNode;
   if(node.depth() == 0){
      // For the root node; simply copy
      m_completeFocusGraph = solver.preprocessedGraph();
      m_focusGraph = solver.preprocessedGraph();
      m_mapToPreprocessed =
          PreprocessedMap({},{}, NodeMap::identity(m_focusGraph.numNodes()),
                          NodeMap::identity(m_focusGraph.numNodes()));
      solver.syncLocalVarsWithGlobal(m_localData);
      setupLPFromScratch(node);
   }
   else if(isChildNode){ //TODO: fix with multithreading initialization
      //We can more efficiently initialize child nodes, as we already have
      //the graphs and the LP in memory which are 'almost' correct,
      //which means we can skip some work

      auto completeGraph =
          constructBranchedFullGraphFromChild(m_completeFocusGraph,
                                              node.branchDecisions(),
                                              node.getNumAddedBranchingDecisions());
      m_completeFocusGraph = completeGraph;
      auto result = preprocessedGraphFromChild(m_focusGraph,
                                               node.branchDecisions(),
                                               node.getNumAddedBranchingDecisions(),
                                               m_mapToPreprocessed.oldToNewIDs,
                                               node.lowerBound()-m_mapToPreprocessed.fixed_sets.size());

      m_focusGraph = result.graph;
      m_mapToPreprocessed.extend(result.map);

      if(m_mapToPreprocessed.newToOldIDs.size() == 0){
         std::cout<<"Empty after preprocessing!\n";
      }
      //TODO: check if preprocessing found new optimal solution, and terminate if preprocessing removed the graph
      solver.writeLocalVarsToGlobal(m_localData,stop);
      setupChildLP(node,solver,result.map);
   }
   else{
      // When not a child node, we do not bother with a scheme to attempt to 'rollback' changes;
      // this is likely to be more pain than it is worth

      // Perform node preprocessing
      auto [completeGraph, result] =
          constructBranchedGraphs(solver.preprocessedGraph(),
                                  node.branchDecisions(), node.lowerBound());
      m_completeFocusGraph = completeGraph;
      m_focusGraph = result.graph;
      m_mapToPreprocessed = result.map;
      if(m_mapToPreprocessed.newToOldIDs.size() == 0){
         std::cout<<"Empty after preprocessing!\n";
      }
      //TODO: check if preprocessing found new optimal solution, and terminate if preprocessing removed the graph
      solver.syncLocalVarsWithGlobal(m_localData);
      setupLPFromScratch(node);
   }

   m_focusNode = node.id();

#ifndef NDEBUG
   {
      //Sanity check if the LP and the complete focus graph were constructed
      //correctly

       auto bounds = m_lpSolver.columnUpperBounds();

       const auto& vars = m_localData.variables();
       assert(bounds.size() == vars.size());
       for(std::size_t i = 0; i < bounds.size(); ++i){
         if(bounds[i].value == 0.0){
            assert(!m_completeFocusGraph.setIsStable(vars[i].set()));
         }else{
            assert(m_completeFocusGraph.setIsStable(vars[i].set()));
         }
       }
   };
#endif
}
void ColorNodeWorker::setupChildLP(BBNode &bb_node, const SolutionData &t_solData,
                                  const PreprocessedMap& t_childMap) {
   if(!t_childMap.removed_nodes.empty()){
      assert(m_focusGraph.numNodes() < m_lpSolver.numRows());

      std::vector<int> permutation(m_lpSolver.numRows(),0);
      for(const auto& node : t_childMap.removed_nodes){
         permutation[m_nodeToLPRow[node.removedNode()]] = -1;
      }
      m_lpSolver.removeRows(permutation);
      // permutation vector now gives old LP row to new LP row mapping e.g.
      // perm[old_row ] returns the new row (note; terms of the OLD node indices)

      NodeMap newNodeToLP = NodeMap::identity(m_focusGraph.numNodes());
      for(std::size_t i = 0; i < m_focusGraph.numNodes(); ++i){
         Node oldNode = t_childMap.newToOldIDs[i];
         std::size_t oldRow = m_nodeToLPRow[oldNode];
         std::size_t newRow = permutation[oldRow] ;
         newNodeToLP[i] = newRow;
      }

      m_nodeToLPRow = newNodeToLP;
      m_LPRowToNode = NodeMap::inverse(m_nodeToLPRow,m_lpSolver.numRows());
      //Save permutation of rows
   }
   auto decisions = bb_node.branchDecisions();
   std::size_t numDecisions = decisions.size();
   std::size_t numAdded = bb_node.getNumAddedBranchingDecisions();

   auto bounds = m_lpSolver.columnUpperBounds();
   const auto& variables = m_localData.variables();
   assert(bounds.size() == variables.size()); //TODO: can fail sometimes (fails on brock200_4.clq; somehow more variables in the LP)
   //Fix variables to zero for columns which do not fit the current subproblem
   for(std::size_t j = 0; j < bounds.size(); ++j){
      if(bounds[j].value == 0.0) continue;
      const auto &set = variables[j].set();
      if(!m_completeFocusGraph.setIsStable(set)){
         m_lpSolver.changeBounds(j,0.0,0.0);
      }
   }
   if(bb_node.depth() != 0){
      SmallBasis sb = bb_node.basis();
      auto mapping = m_localData.getGlobalToLocalMapping();
      for(auto& index : sb.basicCols){
         std::size_t newIndex = mapping.mapGlobalToLocal(index);
         assert(m_localData.mapToGlobalIndex(newIndex) == index);
         assert(newIndex < variables.size());
         index = newIndex;
      }
      LPBasis basis = fixLPBasis(sb,bb_node.previousNodeMap());
      m_lpSolver.setBasis(basis);
   }
}


void ColorNodeWorker::improvementHeuristic(const std::vector<std::size_t>& solVarIndices,
                                           BBNode &t_node,
                                           SolutionData &t_solData,
                                           std::atomic_bool& stop) {

   const auto& preprocessedGraph = t_solData.preprocessedGraph();
   const auto& variables = m_localData.variables();
   SetColoring initialColoring;
   for(const auto& index : solVarIndices){
      initialColoring.addColor(variables[index].set());
   }

   std::size_t numSearchColors = initialColoring.numColors()-1;
   std::size_t lowerBound = 0;
   NodeColoring searchColoring(preprocessedGraph.numNodes(),initialColoring);
   NodeColoring bestColoring = searchColoring;

   std::size_t lb = t_solData.lowerBound();
   std::size_t ub = t_solData.upperBound();
   TabuColoring tabuAlgorithm(preprocessedGraph);
   while( lb <= numSearchColors) {
      if(bestColoring.numColors() <= ub){
         tabuAlgorithm.setMaxIterations(100'000);
      }else{
         tabuAlgorithm.setMaxIterations(5'000);
      }
      std::size_t removeColor = numSearchColors; //We remove the 'highest' color, but different strategies are possible
      searchColoring.setNumColors(numSearchColors);
      for(std::size_t i = 0; i < preprocessedGraph.numNodes(); ++i) {
         if(searchColoring[i] == removeColor) {
            searchColoring[i]--; // We set it to be the next color; again different strategies are possible.
         }
      }
      auto foundColoring = tabuAlgorithm.run(searchColoring);
      if(foundColoring.has_value()) {
         searchColoring = foundColoring.value();
         bestColoring = foundColoring.value();

      }else {
         break;
      }
      --numSearchColors;
   }
   //Add best coloring if it was improved
   if(bestColoring.numColors() < initialColoring.numColors()){
      //std::cout<<"Improved "<<initialColoring.numColors()<<"-coloring into a "<< bestColoring.numColors()<<"-coloring\n";
      SetColoring bestSetColoring(bestColoring);
      for(auto& color : bestSetColoring.colors()) {
         maximizeStableSet(color,preprocessedGraph);
      }

      std::vector<std::size_t> indices = repairColoring(bestSetColoring.colors(),t_solData);
      addSolution(indices,t_node,t_solData,false,stop);
   }else{
      addSolution(solVarIndices,t_node,t_solData,false,stop); //add original solution otherwise
   }

}
void ColorNodeWorker::addSolution(const std::vector<std::size_t>& indices,
                                  BBNode& t_node,
                                  SolutionData& t_solData,
                                  bool doImprovements,
                                  std::atomic_bool& stop) {
   if(doImprovements){
      improvementHeuristic(indices,t_node,t_solData,stop);
   }else{
      synchronizeStastistics(t_node,t_solData,stop); //Maybe only sync statistics if an improving solution was found
      m_localData.addSolution(indices);
      t_solData.writeLocalVarsToGlobal(m_localData,stop);
   }
}
std::size_t
ColorNodeWorker::pickChildNode(NodeChildSelectionStrategy strategy,
                               const std::vector<std::size_t> &children) {
   assert(!children.empty());
   switch (strategy) {
   case NodeChildSelectionStrategy::PREFER_SAME:
      return children.front(); //TODO: A bit ugly; it assumes that the worker also puts the SAME and DIFFER branches in this order. Probably should check here
   case NodeChildSelectionStrategy::PREFER_DIFFER:
      return children.back();
   case NodeChildSelectionStrategy::RANDOMLY:{
      std::uniform_int_distribution<std::size_t> dist(0,children.size()-1);
      return children[dist(m_random_engine)];
   }
   }
}
void ColorNodeWorker::runLoop(SolutionData &t_soldata, std::atomic_bool& stop) {
   while(!stop){
      //we pop the node with best bound when the thread has been idling.
      //This can only happen if the thread did no work or if the last node was cut off
      if(t_soldata.settings().nodeLimit() != NO_NODE_LIMIT &&
          t_soldata.numProcessedNodes() >= t_soldata.settings().nodeLimit()){
         break;
      }
      std::optional<BBNode> node = t_soldata.popNextNode(*this);
      bool processed = false;
      while(!m_cancelCurrentNode && !stop && node.has_value()){
         processed = true;
         processNode(node.value(),t_soldata,stop);
         //Although we might have stopped here already, we want to save generated columns and statistics still!
         t_soldata.synchronizeLocalDataStatistics(m_localData);
         t_soldata.writeLocalVarsToGlobal(m_localData,stop); //Need to do this before branching, (otherwise, we have a race condition)
         bool nodeLimitHit = t_soldata.settings().nodeLimit() != NO_NODE_LIMIT &&
                             t_soldata.numProcessedNodes() >= t_soldata.settings().nodeLimit();
         if(stop || m_cancelCurrentNode || nodeLimitHit){
            bool improved = t_soldata.doRecomputeLowerBound(stop);
            break;
         }
         node = t_soldata.branchAndPopNode(node.value(),*this);
         bool improved = t_soldata.doRecomputeLowerBound(stop);
         if(t_soldata.numProcessedNodes() % t_soldata.settings().nodeDisplayFrequency() == 0 || t_soldata.numProcessedNodes() < 2){
            t_soldata.display(std::cout);
         }

      }
      if(stop){
         break;
      }
      if(!processed){
         std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }
      //TODO: can probably wait using condition variables or something smarter which notifies threads that new B&B nodes have arrived
   }
}
std::size_t ColorNodeWorker::successiveChildrenProcessed() const {
   return m_successiveChildNodesProcessed;
}
std::size_t ColorNodeWorker::mapToGlobalIndex(std::size_t localIndex) const {
   return m_localData.mapToGlobalIndex(localIndex);
}
void ColorNodeWorker::cleanUpOnCancel(BBNode &t_node, SolutionData& t_solData,std::atomic_bool& stop) {
   if(t_node.lowerBound() >= t_solData.upperBoundUnscaled()){
      t_node.setStatus(BBNodeStatus::CUT_OFF);
   }else{
      //The node might have been cancelled after branching vertices were already computed
      //In case we have not solved the problem, we still relay this information
      if(!t_node.branchingData().empty()){
         t_node.setStatus(BBNodeStatus::BRANCHED);
      }else{
         t_node.setStatus(BBNodeStatus::INITIALIZED);
      }
   }
   synchronizeStastistics(t_node,t_solData,stop);
   //TODO: check if we also need to write back variables or if this is done everywhere already
}

} // namespace pcog
