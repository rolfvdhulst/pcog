//
// Created by rolf on 29-11-22.
//

#ifndef PCOG_SRC_BRANCHING_HPP
#define PCOG_SRC_BRANCHING_HPP

#include "pcog/utilities/DenseGraph.hpp"
#include "Preprocessing.hpp"
namespace pcog {
///Represents the last type of branching applied to the node
enum class BranchType { DIFFER, SAME, ROOT };
///Struct which contains the two branching nodes and the type of branch
struct BranchData {
   BranchData() = default;
   BranchData(Node first, Node second, BranchType type)
       : first{first}, second{second}, type{type} {};
   Node first = INVALID_NODE;
   Node second = INVALID_NODE;
   BranchType type = BranchType::ROOT;
};

std::pair<DenseGraph, PreprocessingResult>
constructBranchedGraphs(const DenseGraph &graph,
                        const std::vector<BranchData> &branch_information,
                        std::size_t lower_bound);
DenseGraph
constructBranchedFullGraph(const DenseGraph &graph,
                           const std::vector<BranchData> &branch_information);

DenseGraph
constructBranchedFullGraphFromChild(const DenseGraph& t_fullChildGraph,
                                    const std::vector<BranchData>& branch_information,
                                    std::size_t numNewBranchingConstraints);

PreprocessingResult preprocessedGraphFromChild(const DenseGraph& t_preprocessedChildGraph,
                                               const std::vector<BranchData>& branch_information,
                                               std::size_t numNewBranchingConstraints,
                                               const NodeMap& preprocessedToChildMap,
                                               std::size_t coloring_bound);

///
/// Method which constructs the current graph given a set of branching
/// decisions. Also performs preprocessing again to make pricing as efficient as
/// possible
///@param original Root node graph
///@param branch_information the branching decision taken (in order from first
/// to last)
///@param lower bound on the number of colors needed
/// Returns the processed graph and information on how each node was
/// preprocessed
PreprocessingResult
preprocessBranchedGraph(const DenseGraph &branched_full_graph,
                        const std::vector<BranchData> &branch_information,
                        std::size_t coloring_bound);
};     // namespace pcog
#endif // PCOG_SRC_BRANCHING_HPP
