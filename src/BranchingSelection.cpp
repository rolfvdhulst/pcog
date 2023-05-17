//
// Created by rolf on 4-3-23.
//

#include "pcog/BranchingSelection.hpp"
namespace pcog {
std::size_t countSameAddedTriangles(Node node1, Node node2,
                                    const DenseGraph &graph) {
   assert(!graph.isEdge(node1, node2));
   std::size_t numTriangles = 0;

   // making a small picture helps: when contracting the node, only triangles
   // can be formed by contraction if two edges which were not previously
   // connected to one node are now connected to one node.
   DenseSet neighboursNotInFirst =
       graph.neighbourhood(node2).difference(graph.neighbourhood(node1));
   DenseSet neighboursNotInSecond =
       graph.neighbourhood(node1).difference(graph.neighbourhood(node2));

   // One can count these
   for (const auto &node : neighboursNotInFirst) {
      numTriangles +=
          graph.neighbourhood(node).intersection(neighboursNotInSecond).size();
   }
   return numTriangles;
}

std::vector<ScoredEdge> getAllBranchingEdges(const DenseGraph &graph) {
   std::size_t nodes = graph.numNodes();
   std::vector<ScoredEdge> candidates;
   for (Node i = 0; i < nodes; ++i) {
      for (Node j = 0; j < i; ++j) {
         if (graph.isEdge(i, j)) { // cannot branch on nodes which have edges
            continue;
         }
         candidates.emplace_back(i, j);
      }
   }
   return candidates;
}

void scoreBranchingCandidates(std::vector<ScoredEdge> &candidates,
                              BranchingStrategy strategy,
                              const DenseGraph &graph,
                              LPSolver& t_lpSolver,
                              const NodeMap& t_nodeToLPRow,
                              const std::vector<StableSetVariable>& variables,
                              const NodeMap& mapToPreprocessed,
                              std::size_t numPreprocessedNodes
                              ) {
   switch (strategy) {
   case BranchingStrategy::INTERSECTION_SIZE:
      for (auto &candidate : candidates) {
         candidate.score = static_cast<double>(
             graph.neighbourhood(candidate.node1)
                 .intersection(graph.neighbourhood(candidate.node2))
                 .size());
      }
      return;
   case BranchingStrategy::UNION_SIZE:
      for (auto &candidate : candidates) {
         candidate.score = static_cast<double>(
             graph.neighbourhood(candidate.node1)
                 .setUnion(graph.neighbourhood(candidate.node2))
                 .size());
      }
      return;
   case BranchingStrategy::INTERSECTION_UNION_SIZE:
   {
      std::vector<std::size_t> degrees(graph.numNodes());
      auto dual = t_lpSolver.getDualSolution();

      for(std::size_t i = 0; i < degrees.size(); ++i){
         degrees[i] = graph.neighbourhood(i).size();
      }
      for(std::size_t i = 0; i < degrees.size(); ++i){
         if(fabs(dual[i].value) < 1e-8){
            for(const auto& node : graph.neighbourhood(i)){
               if(degrees[node] > 0){
                  degrees[node] -= 1;
               }
            }
            degrees[i] = 0;
         }
      }


      for (auto &candidate : candidates) {

         std::size_t size_1 =  degrees[candidate.node1];
         std::size_t size_2 =  degrees[candidate.node2];
         candidate.score = static_cast<double>(size_1 + size_2);
         // IS + union size is the same as just counting
         // every element once for both individually
         // elements in the intersection
      }
      return;
   }
   case BranchingStrategy::SYMMETRIC_DIFFERENCE_SIZE:
      for (auto &candidate : candidates) {
         candidate.score = static_cast<double>(
             graph.neighbourhood(candidate.node1)
                 .symmetricDifference(graph.neighbourhood(candidate.node2))
                 .size());
      }
      return;
   case BranchingStrategy::TRIANGLES_ADDED:
      for (auto &candidate : candidates) {
         std::size_t num_triangles_added_differ =
             graph.neighbourhood(candidate.node1)
                 .intersection(graph.neighbourhood(candidate.node2))
                 .size();
         std::size_t num_triangles_added_same =
             countSameAddedTriangles(candidate.node1, candidate.node2, graph);
         candidate.score = static_cast<double>(num_triangles_added_differ +
                                               num_triangles_added_same);
      }
      return;
   case BranchingStrategy::TRIANGLES_ADDED_SCALED:
      scoreScaledTriangles(candidates, graph);
      return;
   case BranchingStrategy::HELDS_RULE:
   {
      auto lpSol = t_lpSolver.getPrimalSolution();
      scoreHeldsRule(candidates,lpSol,variables,mapToPreprocessed);
      return;
   }
   case BranchingStrategy::RANDOMLY:
      return;
   case BranchingStrategy::DUAL_MAX:{
      auto dualValues = t_lpSolver.getDualSolution();
      scoreDualMaximization(candidates,dualValues,t_nodeToLPRow);
      return;
   }

   case BranchingStrategy::DUAL_MIN:{
      auto dualValues = t_lpSolver.getDualSolution();
      scoreDualMinimization(candidates,dualValues,t_nodeToLPRow);
      return;
   }

   case BranchingStrategy::FRACTIONAL:{
      auto lpSol = t_lpSolver.getPrimalSolution();
      scoreFractional(candidates,lpSol,variables,mapToPreprocessed,numPreprocessedNodes);
      return;
   }
   case BranchingStrategy::REMOVAL_SIZE:{
      auto lpSol = t_lpSolver.getPrimalSolution();
      scoreRemovalSize(candidates,lpSol,variables,mapToPreprocessed);
      return;
   }

   case BranchingStrategy::MIN_REMOVAL_SIZE:{
      auto lpSol = t_lpSolver.getPrimalSolution();
      scoreMinRemovalSize(candidates,lpSol,variables,mapToPreprocessed);
      return;
   }
   case BranchingStrategy::SMALL_DIFFERENCE:{
      for(auto& candidate : candidates){
         std::size_t size_1 = graph.neighbourhood(candidate.node1).difference(graph.neighbourhood(candidate.node2)).size();
         std::size_t size_2 = graph.neighbourhood(candidate.node2).difference(graph.neighbourhood(candidate.node1)).size();
         std::size_t min_diff = std::min(size_1,size_2);
         std::size_t sum = graph.neighbourhood(candidate.node1).size() + graph.neighbourhood(candidate.node2).size();
         candidate.score = -double(min_diff) + double(sum)/(2*double(graph.numNodes()));
      }

   }
   }
}
void scoreScaledTriangles(std::vector<ScoredEdge> &candidates,
                          const DenseGraph &graph) {
   // first, compute the number of triangles in both branches for each candidate
   // we also compute the maximum number of triangles
   std::vector<std::pair<std::size_t, std::size_t>> triangles_pairs;
   triangles_pairs.reserve(candidates.size());

   std::size_t max_differ = 1;
   std::size_t max_same = 1;
   for (auto &candidate : candidates) {
      std::size_t num_triangles_added_differ =
          graph.neighbourhood(candidate.node1)
              .intersection(graph.neighbourhood(candidate.node2))
              .size();
      std::size_t num_triangles_added_same =
          countSameAddedTriangles(candidate.node1, candidate.node2, graph);
      triangles_pairs.emplace_back(num_triangles_added_differ,
                                   num_triangles_added_same);
      if (max_differ < num_triangles_added_differ) {
         max_differ = num_triangles_added_differ;
      }
      if (max_same < num_triangles_added_same) {
         max_same = num_triangles_added_same;
      }
   }

   // the final score for same is scaled so that it is 'more fair' to DIFFER
   for (std::size_t i = 0; i < candidates.size(); ++i) {
      double scaled_same_score =
          static_cast<double>(triangles_pairs[i].second*max_differ) / static_cast<double>(max_same);
      candidates[i].score = static_cast<double>(triangles_pairs[i].first) + scaled_same_score;
   }
}

void scoreHeldsRule(std::vector<ScoredEdge> &candidates,
                    const RowVector & lpSolution,
                    const std::vector<StableSetVariable>& variables,
                    const NodeMap& mapToPreprocessed) {
   for (auto &candidate : candidates) {
      Node oldNode1 = mapToPreprocessed[candidate.node1];
      Node oldNode2 = mapToPreprocessed[candidate.node2];
      double numerator = 0.0;
      double denominator = 0.0;
      for (const auto &var : lpSolution) {
         const auto &variable = variables[var.column];
         bool contains_1 = variable.set().contains(oldNode1);
         bool contains_2 = variable.set().contains(oldNode2);
         if (contains_1 && contains_2) {
            numerator += var.value;
         }
         if (contains_1) {
            denominator += var.value;
         }
         if (contains_2) {
            denominator += var.value;
         }
      }
      constexpr double alpha = 0.55;
      // negate the score (we want to be as close as possible to the alpha fraction)
      candidate.score =-fabs(numerator / (0.5 * denominator) - alpha);
   }
}

void scoreDualMaximization(std::vector<ScoredEdge> &t_candidates,
                           const RowVector& t_dualValues,
                           const NodeMap& t_nodeToLPRow) {
   for (auto &candidate : t_candidates) {
      double dual_sum = t_dualValues[t_nodeToLPRow[candidate.node1]].value + t_dualValues[t_nodeToLPRow[candidate.node2]].value;
      candidate.score = dual_sum;
   }
}
void scoreDualMinimization(std::vector<ScoredEdge> &t_candidates,
                           const RowVector& t_dualValues,
                           const NodeMap& t_nodeToLPRow) {
   for (auto &candidate : t_candidates) {
      double dual_sum = t_dualValues[t_nodeToLPRow[candidate.node1]].value + t_dualValues[t_nodeToLPRow[candidate.node2]].value;
      candidate.score = -dual_sum;
   }
}
void scoreFractional(std::vector<ScoredEdge> &candidates,
                     const RowVector & lpSolution,
                     const std::vector<StableSetVariable>& variables,
                     const NodeMap& mapToPreprocessed,
                     std::size_t numPreprocessedGraphNodes) {

   std::vector<double> fractionalities(numPreprocessedGraphNodes, 0.0);
   for (const auto &var : lpSolution) {
      const auto &variable = variables[var.column];
//      assert(var.value <= 1.0 && var.value >= 0.0); //TODO: margins OR rounding
      double fractionality = std::min(1.0 - var.value,var.value);
      for (const auto &node : variable.set()) {
         fractionalities[node] = std::max(fractionalities[node], fractionality);
      }
   }
   for (auto &candidate : candidates) {
      Node oldNode1 = mapToPreprocessed[candidate.node1];
      Node oldNode2 = mapToPreprocessed[candidate.node2];
      candidate.score = fractionalities[oldNode1] + fractionalities[oldNode2];
   }
}
void scoreRemovalSize(std::vector<ScoredEdge> &candidates,
                      const RowVector & lpSolution,
                      const std::vector<StableSetVariable>& variables,
                      const NodeMap& mapToPreprocessed) {

   for (auto &candidate : candidates) {
      Node oldNode1 = mapToPreprocessed[candidate.node1];
      Node oldNode2 = mapToPreprocessed[candidate.node2];
      double sum = 0.0;
      for (const auto &var : lpSolution) {
         const auto &variable = variables[var.column];
         bool contains_1 = variable.set().contains(oldNode1);
         bool contains_2 = variable.set().contains(oldNode2);
         if (contains_1 || contains_2) {
            sum += var.value; // DIFFER is broken if both are true, SAME if
                              // exactly one is true
         }
      }
      candidate.score = sum;
   }
}

void scoreMinRemovalSize(std::vector<ScoredEdge> &candidates,
                         const RowVector & lpSolution,
                         const std::vector<StableSetVariable>& variables,
                         const NodeMap& mapToPreprocessed) {

   for (auto &candidate : candidates) {
      Node oldNode1 = mapToPreprocessed[candidate.node1];
      Node oldNode2 = mapToPreprocessed[candidate.node2];
      double sum = 0.0;
      for (const auto &var : lpSolution) {
         const auto &variable = variables[var.column];
         bool contains_1 = variable.set().contains(oldNode1);
         bool contains_2 = variable.set().contains(oldNode2);
         if (contains_1 || contains_2) {
            sum += var.value; // DIFFER is broken if both are true, SAME if
                              // exactly one is true
         }
      }
      //No violation; we want to disallow pairs with no violation
      if(fabs(sum) <= 1e-8){ //TODO: move tolerances to a central options spot
         sum = 1e9; //penalize significantly
      }
      candidate.score = -sum;
   }
}

std::vector<ScoredEdge> getAllBranchingEdgesViolatedInBoth(
    const DenseGraph &focusGraph, const RowVector &lp_sol,
    const std::vector<StableSetVariable> &variables,
    const NodeMap &mapToFocussed) {
   std::vector<DenseSet> nodeSols(focusGraph.numNodes(),
                                  DenseSet(lp_sol.size()));
   for (std::size_t i = 0; i < lp_sol.size(); ++i) {
      if (lp_sol[i].value != 0.0) {
         const auto &stable_set = variables[lp_sol[i].column];
         for (const auto &node : stable_set.set()) {
            Node mappedNode = mapToFocussed[node];
            if (mappedNode != INVALID_NODE) {
               nodeSols[mappedNode].add(i);
            }
         }
      }
   }
   std::vector<ScoredEdge> edges;
   for (Node i = 1; i < focusGraph.numNodes(); ++i) {
      for (Node j = 0; j < i; ++j) {
         std::size_t first_size = nodeSols[i].size();
         std::size_t second_size = nodeSols[j].size();
         std::size_t intersection_size =
             nodeSols[i].intersection(nodeSols[j]).size();
         // intersection needs to be nonempty for DIFFER to be violated.
         // We need some different set to be in at least one of the two in order
         // to have a stable set which is violated in SAME
         if (!focusGraph.isEdge(i, j) && intersection_size > 0 &&
             intersection_size < std::min(first_size, second_size)) {
            edges.emplace_back(i, j);
         }
      }
   }
   return edges;
}

std::vector<ScoredEdge>::const_iterator
selectBestPair(const std::vector<ScoredEdge> &t_sortedPairs,
               CandidateSelectionStrategy t_selectionStrategy, const RowVector &solution,
               const NodeMap &focusToPreprocessed,
               const std::vector<StableSetVariable> &stable_sets) {
   auto best_pair = t_sortedPairs.end();

   switch (t_selectionStrategy) {
   case CandidateSelectionStrategy::VIOLATED_IN_BOTH: {
      for (auto pair_it = t_sortedPairs.begin(); pair_it != t_sortedPairs.end();
           ++pair_it) {
         for (std::size_t i = 1; i < solution.size(); ++i) { //TODO: limit LP solution to only fractional variables
            for (std::size_t j = 0; j < i; ++j) {
               const auto &set_1 =
                   stable_sets[solution[i].column]; // TODO: fix conversions
                                                    // somehow?
               const auto &set_2 = stable_sets[solution[j].column];

               Node node1 = focusToPreprocessed[pair_it->node1];
               Node node2 = focusToPreprocessed[pair_it->node2];
               assert(node1 != INVALID_NODE && node2 != INVALID_NODE);

               bool set1_contains_node1 = set_1.set().contains(node1);
               bool set1_contains_node2 = set_1.set().contains(node2);
               bool set2_contains_node1 = set_2.set().contains(node1);
               bool set2_contains_node2 = set_2.set().contains(node2);
               int sum = set1_contains_node1 + set1_contains_node2 +
                         set2_contains_node1 + set2_contains_node2;
               if (sum == 3) {
                  return pair_it;
               }
            }
         }
      }
   } break;
   case CandidateSelectionStrategy::VIOLATED_IN_SAME: // TODO fix
   case CandidateSelectionStrategy::VIOLATED_IN_DIFFER:
   case CandidateSelectionStrategy::VIOLATED_IN_ONE:
      break;
   case CandidateSelectionStrategy::FIRST:
      break;
   }
   if (best_pair == t_sortedPairs.end()) {
      std::cout << "Could not find unviolated LP pair, picking the first one "
                   "arbitrarily!\n";
      best_pair = t_sortedPairs.begin();
   }
   return best_pair;
}

} // namespace pcog