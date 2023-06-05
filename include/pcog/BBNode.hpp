//
// Created by rolf on 6-12-22.
//

#ifndef PCOG_SRC_BBNODE_HPP
#define PCOG_SRC_BBNODE_HPP
#include "Branching.hpp"
#include "LPSolver.hpp"
#include "utilities/RbTree.hpp"
#include <queue>
#include <utility>

namespace pcog {
enum class BBNodeStatus {
   INITIALIZED,
   PROCESSING,
   BRANCHED, // Node was solved, and branching was performed
   CUT_OFF,  // Node was cut off e.g. the lower bound was greater than the upper
};

using node_id = std::size_t;
static constexpr node_id INVALID_BB_NODE = std::numeric_limits<node_id>::max();

class ColorNodeWorker;

class BBNode {
 public:
   BBNode(node_id id, node_id parent_id, node_id depth,
          double fractionalLowerBound, std::size_t lowerBound,
          std::vector<BranchData> branchInfo, std::size_t numAddedBranches,
          SmallBasis t_basis, NodeMap t_map)
       : m_id{id}, m_parent_id{parent_id}, m_depth{depth},
         m_branchingDecisions{std::move(branchInfo)},
         m_numAddedBranchingDecisions{numAddedBranches},
         m_status{BBNodeStatus::INITIALIZED},
         m_fractionalLowerBound{fractionalLowerBound}, m_lowerBound{lowerBound},
         m_firstBranchNode{INVALID_NODE}, m_secondBranchNode{INVALID_NODE},
         m_previousNodeMap(std::move(t_map)),
         m_initialBasis(std::move(t_basis))
         {};
   [[nodiscard]] node_id id() const { return m_id; }
   [[nodiscard]] double fractionalLowerBound() const {
      return m_fractionalLowerBound;
   }
   [[nodiscard]] std::size_t lowerBound() const { return m_lowerBound; }
   [[nodiscard]] node_id parent() const { return m_parent_id; }
   [[nodiscard]] const std::vector<node_id> &children() const {
      return m_childrenID;
   }
   [[nodiscard]] std::size_t depth() const { return m_depth; }
   [[nodiscard]] BBNodeStatus status() const { return m_status; }

   [[nodiscard]] std::vector<BranchData> branchDecisions() const {
      return m_branchingDecisions;
   }
   [[nodiscard]] Node firstBranchingNode() const { return m_firstBranchNode; }
   [[nodiscard]] Node secondBranchNode() const { return m_secondBranchNode; }

   bool updateLowerBound(double lowerBound) {
      m_fractionalLowerBound = std::max(m_fractionalLowerBound, lowerBound);
      auto lb = static_cast<std::size_t>(std::ceil(lowerBound));// TODO: fix conversion
      bool improved = lb > m_lowerBound;
      m_lowerBound = improved ? lb : m_lowerBound;

      return improved;
   }
   void setBranchingNodes(Node a, Node b) {
      m_firstBranchNode = a;
      m_secondBranchNode = b;
   }
   void setStatus(BBNodeStatus t_status) { m_status = t_status; }
   /// Check if a stable set meets the branching decisions taken for this node.
   /// Used only for debugging.
   //TODO: sets with multiple SAME constraints can still be problematic,
   //really need to check if set is stable in current local graph instead (Remove this function)
   [[nodiscard]] bool verifyStableSet(const DenseSet &set) const {
      for (const auto &data : m_branchingDecisions) {
         if ((data.type == BranchType::SAME &&
              (set.contains(data.first) != set.contains(data.second))) ||
             (data.type == BranchType::DIFFER && set.contains(data.first) &&
              set.contains(data.second))) {
            return false;
         }
      }
      return true;
   }
   [[nodiscard]] SmallBasis basis() const;

   void setBasis(SmallBasis basis);
   [[nodiscard]] NodeMap previousNodeMap() const {return m_previousNodeMap;}


   RbTreeLinks<int64_t> &getLowerLinks() { return m_lowerLinks; }
   [[nodiscard]] const RbTreeLinks<int64_t> &getLowerLinks() const { return m_lowerLinks; }

   [[nodiscard]] const std::vector<std::vector<BranchData>> &
   branchingData() const {
      return m_branchingData;
   }
   void setBranchingData(std::vector<std::vector<BranchData>> branchData) {
      m_branchingData = std::move(branchData);
   }

   std::size_t getNumAddedBranchingDecisions() const;

 private:
   node_id m_id;
   node_id m_parent_id;
   std::vector<node_id> m_childrenID;
   std::size_t m_depth;

   std::vector<BranchData> m_branchingDecisions;
   std::size_t m_numAddedBranchingDecisions;

   BBNodeStatus m_status;
   double m_fractionalLowerBound;
   std::size_t m_lowerBound;

   Node m_firstBranchNode;
   Node m_secondBranchNode;
   std::vector<std::vector<BranchData>> m_branchingData;

   NodeMap m_previousNodeMap; // needed to interpret previous basis
   SmallBasis m_initialBasis;

   RbTreeLinks<int64_t> m_lowerLinks;
};

enum class NodeSelectionStrategy { DFS, BFS, BEST_BOUND };
class BBTree {
 public:
   class LowerRbTree;
   /// Clear all the information stored in this object.
   void clear();
   void createRootNode(std::size_t numRootGraphNodes);

   [[nodiscard]] bool hasOpenNodes() const;

   const BBNode& peekNode(node_id) const;
   BBNode &&popNextNode();
   BBNode &&popNodeWithID(node_id id);

   /**
    * Remove branch-and-bound nodes whoms
    * lower bound is greater or equal to the newly found upper bound
    * @param numColors
    */
   void pruneUpperBound(std::size_t numColors);
   std::vector<node_id> createChildren(const BBNode &t_parentNode,
                       ColorNodeWorker &t_nodeWorker);

   [[nodiscard]] std::size_t numOpenNodes() const;
   [[nodiscard]] std::size_t numTotalNodes() const;
   [[nodiscard]] std::size_t numProcessedNodes() const;
   [[nodiscard]] double fractionalLowerBound() const;
   [[nodiscard]] std::size_t lowerBound() const;

 private:
   node_id createNode(const BBNode &t_parentNode, ColorNodeWorker &worker,
                   std::vector<BranchData> t_addedBranchingDecisions);

   /// Newly links the adding it to the distributed ordered set structures.
   /// This ensures we can efficiently select the next open node and keep track
   /// of the lower bound
   void addToTrees(node_id t_nodeId);
   void removeFromTrees(node_id t_nodeId);

   std::vector<BBNode> m_nodes;
   std::priority_queue<node_id, std::vector<node_id>, std::greater<>>
       m_freeSlots;

   NodeSelectionStrategy m_selection_strategy;

   int64_t lowerRoot = RbTreeLinks<int64_t>::noLink();
   int64_t lowerMin = RbTreeLinks<int64_t>::noLink();

   std::size_t m_totalNodes = 0;
};

} // namespace pcog
#endif // PCOG_SRC_BBNODE_HPP
