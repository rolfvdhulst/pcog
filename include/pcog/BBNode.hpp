//
// Created by rolf on 6-12-22.
//

#ifndef PCOG_SRC_BBNODE_HPP
#define PCOG_SRC_BBNODE_HPP
#include "Branching.hpp"
#include <queue>
#include <utility>

namespace pcog {
enum class BBNodeStatus {
   INITIALIZED,
   PROCESSING,
   BRANCHED, // Node was solved, but
   CUT_OFF,  // Node was cut off e.g. the lower bound was greater than the upper
             // bound
};

using node_id =
    std::size_t; // TODO: rename to something not confusable with graph nodes
static constexpr node_id INVALID_BB_NODE = std::numeric_limits<node_id>::max();

class BBNode {
 public:
   BBNode(node_id id, node_id parent_id, node_id depth,
          double fractionalLowerBound, std::size_t lowerBound,
          std::vector<BranchData> branchInfo)
       : m_id{id}, m_parent_id{parent_id}, m_depth{depth},
         m_branchingDecisions{std::move(branchInfo)}, m_status{BBNodeStatus::INITIALIZED},
         m_fractionalLowerBound{fractionalLowerBound}, m_lowerBound{lowerBound},
         m_firstBranchNode{INVALID_NODE},m_secondBranchNode{INVALID_NODE}
         {};
   [[nodiscard]] node_id id() const { return m_id; }
   [[nodiscard]] double fractionalLowerBound() const { return m_fractionalLowerBound; }
   [[nodiscard]] std::size_t lowerBound() const {return m_lowerBound;}
   [[nodiscard]] node_id parent() const { return m_parent_id; }
   [[nodiscard]] const std::vector<node_id> &children() const {
      return m_childrenID;
   }
   [[nodiscard]] std::size_t depth() const { return m_depth; }
   [[nodiscard]] BBNodeStatus status() const { return m_status; }

   [[nodiscard]] std::vector<BranchData> branchDecisions() const {return m_branchingDecisions;}
   [[nodiscard]] Node firstBranchingNode() const { return m_firstBranchNode; }
   [[nodiscard]] Node secondBranchNode() const { return m_secondBranchNode; }

   void updateLowerBound(double lowerBound) {
      m_fractionalLowerBound = std::max(m_fractionalLowerBound,lowerBound);
      m_lowerBound = std::max(m_lowerBound,static_cast<std::size_t>(std::ceil(lowerBound)));//TODO: fix conversion?
   }
   void setStatus(BBNodeStatus t_status) {m_status = t_status;}
   /// Check if a stable set meets the branching decisions taken for this node. Used only for debugging.
   [[nodiscard]] bool verifyStableSet(const DenseSet& set) const{
      for(const auto& data : m_branchingDecisions){
         if((data.type == BranchType::SAME && (set.contains(data.first) != set.contains(data.second))) ||
             (data.type == BranchType::DIFFER && set.contains(data.first) && set.contains(data.second))){
            return false;
         }
      }
      return true;
   }
 private:
   node_id m_id;
   node_id m_parent_id;
   std::vector<node_id> m_childrenID;
   std::size_t m_depth;
   std::vector<BranchData> m_branchingDecisions;
   BBNodeStatus m_status;
   double m_fractionalLowerBound;
   std::size_t m_lowerBound;
   Node m_firstBranchNode;
   Node m_secondBranchNode;
   // TODO: store basis and other information
};

enum class NodeSelectionStrategy { DFS, BFS, BEST_BOUND };
class BBTree {
 public:
   /// Clear all the information stored in this object.
   void clear();
   void createRootNode();

   [[nodiscard]] bool hasOpenNodes() const;

   BBNode &popNextNode();

   void createChildren(node_id t_node);

 private:
   // TODO: for now we store the entire b&b tree, but by increasing the
   // per-memory, we can get away with only storing the open nodes Evaluate if
   // the memory savings are worth it.
   std::vector<BBNode> m_node_data;
   std::vector<node_id> m_open_nodes; // sorted by selection priority,
   NodeSelectionStrategy m_selection_strategy;
};

} // namespace pcog
#endif // PCOG_SRC_BBNODE_HPP
