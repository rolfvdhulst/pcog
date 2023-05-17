//
// Created by rolf on 6-12-22.
//

#include <utility>

#include "pcog/BBNode.hpp"
#include "pcog/ColorNodeWorker.hpp"
namespace pcog {

template<>
struct RbTreeTraits<BBTree::LowerRbTree> {
   using KeyType = std::tuple<double,int64_t>;
   using LinkType = int64_t;
};
class BBTree::LowerRbTree : public RbTree<LowerRbTree> {
   BBTree* m_tree;
 public:
   explicit LowerRbTree(BBTree* t_tree) : RbTree<LowerRbTree>(t_tree->lowerRoot,
                                                              t_tree->lowerMin),
       m_tree{t_tree}{}

   RbTreeLinks<int64_t>& getRbTreeLinks(int64_t t_node){
      return m_tree->m_nodes[t_node].getLowerLinks();
   }
   [[nodiscard]] const RbTreeLinks<int64_t>& getRbTreeLinks(int64_t t_node) const {
      return m_tree->m_nodes[t_node].getLowerLinks();
   }
   [[nodiscard]] std::tuple<double,int64_t> getKey(int64_t t_node) const {
      return std::make_tuple(m_tree->m_nodes[t_node].fractionalLowerBound(),t_node);
   }
};

void BBTree::clear() {
   while(!m_freeSlots.empty()){ //TODO: can clear more efficiently
      m_freeSlots.pop();
   }
   m_nodes.clear();
   m_selection_strategy = NodeSelectionStrategy::DFS; // TODO: reset selection strategy?
}
void BBTree::createRootNode(std::size_t numRootGraphNodes) {
   assert(m_nodes.empty());
   assert(m_freeSlots.empty());


   m_nodes.emplace_back(0, INVALID_BB_NODE, 0, 0.0, 0,
                            std::vector<BranchData>{BranchData(
                                INVALID_NODE, INVALID_NODE, BranchType::ROOT)},1
                        ,SmallBasis{},NodeMap::identity(numRootGraphNodes));
   m_totalNodes = 1;
   addToTrees(0);
}
bool BBTree::hasOpenNodes() const { return numOpenNodes() != 0; }

BBNode&& BBTree::popNextNode() {
   int64_t node = lowerMin;
   removeFromTrees(node);
   return std::move(m_nodes[node]);
}

BBNode&& BBTree::popNodeWithID(node_id id){

   removeFromTrees(id);
   return std::move(m_nodes[id]);
}
void BBTree::pruneUpperBound(std::size_t numColors) {
   //TODO: pruned nodes are counted as 'processed' nodes in b&b loop visualization. Fix

   //Iterate from largest to smallest lowerBound, pruning the nodes with too high lowerBounds away
   LowerRbTree rbTree(this);
   int64_t node = rbTree.upper_bound();

   while(node != RbTreeLinks<int64_t>::noLink()){
      if(m_nodes[node].lowerBound() < numColors) break;
      int64_t next = rbTree.predecessor(node); //First get the predecessor, before invalidating the node
      removeFromTrees(node);
      node = next;
   }
}
std::vector<node_id> BBTree::createChildren(const BBNode& t_parentNode, ColorNodeWorker &t_nodeWorker) {
   std::vector<node_id> children;
   for(const auto& data : t_parentNode.branchingData()){
      node_id node = createNode(t_parentNode,t_nodeWorker,data);
      children.push_back(node);
   }
   return children;
}
std::size_t BBTree::numOpenNodes() const {
   return m_nodes.size()-m_freeSlots.size();
}
std::size_t BBTree::numTotalNodes() const {
   return m_totalNodes;
}
std::size_t BBTree::numProcessedNodes() const {
   return numTotalNodes()-numOpenNodes();
}
std::size_t BBTree::lowerBound() const {
   std::size_t lb = lowerMin == RbTreeLinks<int64_t>::noLink() ? std::numeric_limits<std::size_t>::max()
                                                               : m_nodes[lowerMin].lowerBound();
   return lb;
}
double BBTree::fractionalLowerBound() const {
   double lb = lowerMin == RbTreeLinks<int64_t>::noLink() ? std::numeric_limits<double>::infinity() : m_nodes[lowerMin].fractionalLowerBound();
   return lb;
}
node_id BBTree::createNode(const BBNode &t_parentNode, ColorNodeWorker &t_nodeWorker,
                        std::vector<BranchData> t_addedBranchinDecisions) {
   node_id pos;
   std::vector<BranchData> branchingPath = t_parentNode.branchDecisions();
   branchingPath.insert(branchingPath.end(),t_addedBranchinDecisions.begin(),t_addedBranchinDecisions.end());
   if(m_freeSlots.empty()){
      pos = m_nodes.size();
      m_nodes.emplace_back(m_totalNodes, t_parentNode.id(), t_parentNode.depth()+1,
                           t_parentNode.fractionalLowerBound(),
                           t_parentNode.lowerBound(), branchingPath, t_addedBranchinDecisions.size(),
                           t_parentNode.basis(),t_nodeWorker.focusToPreprocessed());
   }else{
      pos = m_freeSlots.top();
      m_freeSlots.pop();
      m_nodes[pos] = BBNode(m_totalNodes, t_parentNode.id(), t_parentNode.depth()+1,
                            t_parentNode.fractionalLowerBound(),
                            t_parentNode.lowerBound(), branchingPath,t_addedBranchinDecisions.size(),
                            t_parentNode.basis(),t_nodeWorker.focusToPreprocessed());
   }
   ++m_totalNodes;
   addToTrees(pos);
   return pos;
}
void BBTree::addToTrees(node_id t_nodeId) {
   LowerRbTree rbTree(this);
   rbTree.insert(static_cast<int64_t>(t_nodeId));
}
void BBTree::removeFromTrees(node_id t_nodeId) {
   LowerRbTree rbTree(this);
   rbTree.erase(static_cast<int64_t>(t_nodeId));
   m_freeSlots.push(t_nodeId);
}
const BBNode &BBTree::peekNode(node_id t_nodeId) const {
   assert(t_nodeId < m_nodes.size());
   //TODO: assert that node is not 'freed' up
   return m_nodes[t_nodeId];
}
SmallBasis BBNode::basis() const { return m_initialBasis; }
std::size_t BBNode::getNumAddedBranchingDecisions() const { return m_numAddedBranchingDecisions; }
void BBNode::setBasis(SmallBasis basis) {
   m_initialBasis = std::move(basis);
}
} // namespace pcog