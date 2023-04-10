//
// Created by rolf on 27-11-22.
//

#ifndef PCOG_SRC_PREPROCESSING_HPP
#define PCOG_SRC_PREPROCESSING_HPP

#include "pcog/utilities/DenseGraph.hpp"
#include "pcog/utilities/NodeMap.hpp"
#include "pcog/utilities/Coloring.hpp"

namespace pcog {

enum class PreprocessedReason { LOW_DEGREE, DOMINATED_NODE };

/// A class which holds a node which was removed during preprocessing
class PreprocessedNode {
 public:
   /// Construct a Preprocessed node
   /// \param t_node Node which was preprocessed
   /// \param t_reason Reason for its removal
   /// \param t_dominatedByNode If the node was dominated, this field holds the
   /// node which dominated it. Otherwise, it holds an INVALID_NODE
   PreprocessedNode(Node t_node, PreprocessedReason t_reason,
                    Node t_dominatedByNode = INVALID_NODE);

   /// \return The node which was removed
   [[nodiscard]] Node removedNode() const;
   /// \return The reason why this node was preprocessed
   [[nodiscard]] PreprocessedReason removedReason() const;
   /// \return The node which dominates this node. In case this node was not
   /// dominated, this should be an INVALID_NODE.
   [[nodiscard]] Node dominatingNode() const;

 private:
   Node m_node;
   PreprocessedReason m_reason;
   Node m_dominated_by;
};

/// This class holds the preprocessed nodes and a bijection between the new and
/// old nodes
struct PreprocessedMap {
   PreprocessedMap() = default;
   PreprocessedMap(std::vector<PreprocessedNode> removed_nodes,
                   NodeMap newToOld, NodeMap oldToNew);
   void clear();
   std::vector<PreprocessedNode> removed_nodes;
   NodeMap newToOldIDs;
   NodeMap oldToNewIDs;
};

/// This class holds the preprocessed graph and a mapping with the old vertices
struct PreprocessingResult {
   PreprocessingResult(DenseGraph graph, PreprocessedMap map);
   DenseGraph graph;
   PreprocessedMap map;
};

/// Preprocess the original root node graph
/// \param graph The graph to preprocess
/// \return The preprocessing result with graph and map.
PreprocessingResult preprocessOriginalGraph(const DenseGraph &graph);

/// Preprocess low degree vertices from the graph, given a clique which is not
/// to be removed. This function modifies the present nodes, and only considers
/// the graph induced by present_nodes.
/// \param graph The graph to modify. Only the graph induced by present_nodes is
/// considered.
/// \param present_nodes The nodes which induce the subgraph. Modified if any nodes are removed.
/// \param clique The clique to use for preprocessing
/// \return A vector with the removed nodes
std::vector<PreprocessedNode>
removeLowDegreeVerticesClique(const DenseGraph &graph, DenseSet &present_nodes,
                              const DenseSet &clique);

/// Preprocess dominated vertices from the graph, given a clique which is not
/// to be removed. This function modifies the present nodes, and only considers
/// the graph induced by present_nodes.
/// \param graph The graph to modify. Only the graph induced by present_nodes is
/// considered.
/// \param present_nodes The nodes which induce the subgraph. Modified if any nodes are removed.
/// \param checkForDominatedVertices The vertices which should be checked if they are dominated by another vertex
/// \return A vector with the removed nodes
std::vector<PreprocessedNode>
removeDominatedVerticesClique(const DenseGraph &graph, DenseSet &present_nodes,
                              const DenseSet &checkForDominatedVertices);

/// Preprocess low degree vertices from the graph.  This function modifies the
/// present nodes, and only considers the graph induced by present_nodes.
/// \param graph The graph to modify. Only the graph induced by present_nodes is
/// considered.
/// \param present_nodes The nodes which induce the subgraph. Modified if any nodes are removed.
/// \param clique The clique to use for preprocessing
/// \return A vector with the removed nodes
std::vector<PreprocessedNode>
removeLowDegreeVertices(const DenseGraph &graph, DenseSet &present_nodes,
                        std::size_t coloring_bound);

/// Preprocess strictly dominated vertices from the graph.  This function modifies the
/// present nodes, and only considers the graph induced by present_nodes.
/// \param graph The graph to modify. Only the graph induced by present_nodes is
/// considered.
/// \param present_nodes The nodes which induce the subgraph. Modified if any nodes are removed.
/// \param clique The clique to use for preprocessing
/// \return A vector with the removed nodes
std::vector<PreprocessedNode>
removeStrictlyDominatedVertices(const DenseGraph &graph,
                                DenseSet &present_nodes);

NodeColoring extendColoring(const NodeColoring& coloring,
                            const PreprocessedMap &map,
                            const DenseGraph &originalGraph);

} // namespace pcog
#endif // PCOG_SRC_PREPROCESSING_HPP
