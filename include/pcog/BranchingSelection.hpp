//
// Created by rolf on 4-3-23.
//

#ifndef PCOG_INCLUDE_PCOG_BRANCHINGSELECTION_HPP
#define PCOG_INCLUDE_PCOG_BRANCHINGSELECTION_HPP

#include "ColorSolver.hpp"
#include "utilities/DenseGraph.hpp"
#include "Settings.hpp"
#include <iostream>

namespace pcog {


struct ScoredEdge {
   ScoredEdge() : node1{INVALID_NODE}, node2{INVALID_NODE}, score{0.0} {};
   ScoredEdge(Node t_node1, Node t_node2)
       : node1{t_node1}, node2{t_node2}, score{0.0} {};
   ScoredEdge(Node t_node1, Node t_node2, double t_score)
       : node1{t_node1}, node2{t_node2}, score{t_score} {};
   Node node1;
   Node node2;
   double score;
};

std::vector<ScoredEdge>::const_iterator
selectBestPair(const std::vector<ScoredEdge> &t_sortedPairs,
               CandidateSelectionStrategy t_selectionStrategy, const RowVector &lp_sol,
               const NodeMap &focusToPreprocessed,
               const std::vector<StableSetVariable> &stable_sets);
std::size_t countSameAddedTriangles(Node node1, Node node2,
                                    const DenseGraph &graph);
std::vector<ScoredEdge> getAllBranchingEdges(const DenseGraph &graph);

std::vector<ScoredEdge> getAllBranchingEdgesViolatedInBoth(
    const DenseGraph &graph, const RowVector &lp_sol,
    const std::vector<StableSetVariable> &variables, const NodeMap &toFocussed);

void scoreBranchingCandidates(std::vector<ScoredEdge> &candidates,
                              BranchingStrategy strategy,
                              const DenseGraph &graph,
                              const RowVector& t_primalSol,
                              const RowVector& t_dualSol,
                              const NodeMap& t_nodeToLPRow,
                              const std::vector<StableSetVariable>& variables,
                              const NodeMap& mapToPreprocessed,
                              std::size_t numPreprocessedNodes);

void scoreScaledTriangles(std::vector<ScoredEdge> &edges,
                          const DenseGraph &graph);
void scoreHeldsRule(std::vector<ScoredEdge> &candidates,
                    const RowVector & lpSolution,
                    const std::vector<StableSetVariable>& variables,
                    const NodeMap& mapToPreprocessed);

void scoreDualMaximization(std::vector<ScoredEdge> &t_candidates,
                           const RowVector & t_dualValues,
                           const NodeMap & t_nodeToLPRow);
void scoreDualMinimization(std::vector<ScoredEdge> &t_candidates,
                           const RowVector & t_dualValues,
                           const NodeMap& t_nodeToLPRow
                           );

void scoreRemovalSize(std::vector<ScoredEdge> &candidates,
                      const RowVector & lpSolution,
                      const std::vector<StableSetVariable>& variables,
                      const NodeMap& mapToPreprocessed);

void scoreMinRemovalSize(std::vector<ScoredEdge> &candidates,
                         const RowVector & lpSolution,
                         const std::vector<StableSetVariable>& variables,
                         const NodeMap& mapToPreprocessed);

void scoreFractional(std::vector<ScoredEdge> &candidates,
                     const RowVector & lpSolution,
                     const std::vector<StableSetVariable>& variables,
                     const NodeMap& mapToPreprocessed,
                     std::size_t numPreprocessedGraphNodes);

} // namespace pcog
#endif // PCOG_INCLUDE_PCOG_BRANCHINGSELECTION_HPP
