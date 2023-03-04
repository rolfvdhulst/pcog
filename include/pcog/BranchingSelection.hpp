//
// Created by rolf on 4-3-23.
//

#ifndef PCOG_INCLUDE_PCOG_BRANCHINGSELECTION_HPP
#define PCOG_INCLUDE_PCOG_BRANCHINGSELECTION_HPP

#include "ColorSolver.hpp"
#include "DenseGraph.hpp"

#include <iostream>

namespace pcog {
enum class BranchingStrategy : int {
   INTERSECTION_SIZE = 0,
   UNION_SIZE = 1,
   INTERSECTION_UNION_SIZE = 2,
   SYMMETRIC_DIFFERENCE_SIZE = 3,
   TRIANGLES_ADDED = 4,
   TRIANGLES_ADDED_SCALED = 5,
   HELDS_RULE = 6,
   RANDOMLY = 7,
   DUAL_MAX = 8,
   DUAL_MIN = 9,
   FRACTIONAL = 10,
   REMOVAL_SIZE = 11,
   MIN_REMOVAL_SIZE = 12,
   MAXIMUM_STRATEGY = 12 // TODO MAKE SURE TO UPDATE THIS CORRECTLY WHEN
                         // UPDATING THE STRATEGIES!
};

enum class SelectionStrategy {
   VIOLATED_IN_BOTH = 0,
   VIOLATED_IN_SAME = 1,
   VIOLATED_IN_DIFFER = 2,
   VIOLATED_IN_ONE = 3,
   FIRST = 4
};

struct ScoredEdge {
   ScoredEdge() : node1{INVALID_NODE}, node2{INVALID_NODE}, score{0.0} {};
   ScoredEdge(Node node1, Node node2)
       : node1{node1}, node2{node2}, score{0.0} {};
   ScoredEdge(Node node1, Node node2, double score)
       : node1{node1}, node2{node2}, score{score} {};
   Node node1;
   Node node2;
   double score;
};

std::vector<ScoredEdge>::const_iterator
selectBestPair(const std::vector<ScoredEdge> &t_sortedPairs,
               SelectionStrategy t_selectionStrategy, const RowVector &lp_sol,
               const NodeMap &focusToPreprocessed,
               const std::vector<StableSetVariable> &stable_sets);
std::size_t countSameAddedTriangles(Node node1, Node node2,
                                    const DenseGraph &graph);
std::vector<ScoredEdge> getAllBranchingEdges(const DenseGraph &graph);

std::vector<ScoredEdge> getAllBranchingEdgesViolatedInBoth(
    const DenseGraph &graph, const RowVector &lp_sol,
    const std::vector<StableSetVariable> &variables, const NodeMap &toFocussed);

void scoreBranchingCandidates(std::vector<ScoredEdge> &edges,
                              BranchingStrategy strategy,
                              const DenseGraph &graph);
void scoreScaledTriangles(std::vector<ScoredEdge> &edges,
                          const DenseGraph &graph);
//void scoreHeldsRule(std::vector<ScoredEdge> &edges, const DenseGraph &graph,
//                    SCIPProblem *problem, SCIP *scip);
//void scoreDualMaximization(std::vector<ScoredEdge> &edges,
//                           const DenseGraph &graph, SCIPProblem *problem,
//                           SCIP *scip);
//void scoreDualMinimization(std::vector<ScoredEdge> &edges,
//                           const DenseGraph &graph, SCIPProblem *problem,
//                           SCIP *scip);
//void scoreRemovalSize(std::vector<ScoredEdge> &candidates,
//                      const DenseGraph &graph, SCIPProblem *probdata,
//                      SCIP *scip);
//void scoreMinRemovalSize(std::vector<ScoredEdge> &candidates,
//                         const DenseGraph &graph, SCIPProblem *probdata,
//                         SCIP *scip);
//
//void scoreFractional(std::vector<ScoredEdge> &candidates,
//                     const DenseGraph &graph, SCIPProblem *probdata,
//                     SCIP *scip);

} // namespace pcog
#endif // PCOG_INCLUDE_PCOG_BRANCHINGSELECTION_HPP
