//
// Created by rolf on 6-12-22.
//

#include "pcog/BBNode.hpp"
namespace pcog {

void BBTree::clear() {
   m_node_data.clear();
   m_open_nodes.clear();
   m_selection_strategy =
       NodeSelectionStrategy::DFS; // TODO: reset selection strategy?
}
void BBTree::createRootNode() {
   assert(m_node_data.empty());
   assert(m_open_nodes.empty());

   m_node_data.emplace_back(0, INVALID_BB_NODE, 0, 0.0, 0,
                            std::vector<BranchData>{BranchData(
                                INVALID_NODE, INVALID_NODE, BranchType::ROOT)});
   m_open_nodes.emplace_back(0);
}
bool BBTree::hasOpenNodes() const { return !m_open_nodes.empty(); }
BBNode &BBTree::popNextNode() {
   node_id node = m_open_nodes.back();
   m_open_nodes.pop_back();
   return m_node_data[node];
}
void BBTree::createChildren(node_id t_node) {
   BBNode &parent_data = m_node_data[t_node];

   {
      node_id left_id = m_node_data.size();
      // insert node according to selection rule
      std::vector<BranchData> branchingPath = parent_data.branchDecisions();
      branchingPath.emplace_back(parent_data.firstBranchingNode(),
                                 parent_data.secondBranchNode(),
                                 BranchType::DIFFER);
      BBNode left_node(left_id, t_node, parent_data.depth(),
                       parent_data.fractionalLowerBound(),
                       parent_data.lowerBound(), branchingPath);
      m_node_data.push_back(left_node);
      m_open_nodes.push_back(left_id); // TODO: insert according to node
                                       // ordering
   }
   {
      node_id right_id = m_node_data.size();
      // insert node according to selection rule
      std::vector<BranchData> branchingPath = parent_data.branchDecisions();
      branchingPath.emplace_back(parent_data.firstBranchingNode(),
                                 parent_data.secondBranchNode(),
                                 BranchType::SAME);
      BBNode right_node(right_id, t_node, parent_data.depth(),
                        parent_data.fractionalLowerBound(),
                        parent_data.lowerBound(), branchingPath);
      m_node_data.push_back(right_node);
      m_open_nodes.push_back(
          right_id); // TODO: insert according to node ordering
   }
}
} // namespace pcog