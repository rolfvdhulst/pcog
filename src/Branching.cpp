//
// Created by rolf on 29-11-22.
//

#include "pcog/Branching.hpp"
#include <numeric>
namespace pcog {

std::pair<DenseGraph, PreprocessingResult>
constructBranchedGraphs(const DenseGraph &graph,
                        const std::vector<BranchData> &branch_information,
                        std::size_t lower_bound) {
   DenseGraph fullGraph = constructBranchedFullGraph(graph, branch_information);
   PreprocessingResult result =
       preprocessBranchedGraph(fullGraph, branch_information, lower_bound);
   return {fullGraph, result};
}
DenseGraph
constructBranchedFullGraph(const DenseGraph &graph,
                           const std::vector<BranchData> &branch_information) {
   DenseGraph branched_full_graph(graph);
   std::size_t num_nodes = graph.numNodes();

   std::vector<Node> representative(num_nodes);
   std::iota(representative.begin(), representative.end(), 0);
   std::vector<std::vector<Node>> same_classes(num_nodes, {0});
   for (std::size_t i = 0; i < num_nodes; i++) {
      same_classes[i][0] = i;
   }
   for (const auto &data : branch_information) {
      assert(data.type == BranchType::ROOT ||
             !graph.isEdge(data.first, data.second));
      assert(data.type == BranchType::ROOT ||
             !graph.isEdge(data.second, data.first));
      assert(data.type == BranchType::ROOT ||
             representative[data.first] == data.first);
      assert(data.type == BranchType::ROOT ||
             representative[data.second] == data.second);

#ifndef NDEBUG
      // Check if all same classes are still correct
      for (std::size_t i = 0; i < num_nodes; ++i) {
         if (!same_classes[i].empty()) {
            for (const auto &node : same_classes[i]) {
               assert(branched_full_graph.neighbourhood(node) ==
                      branched_full_graph.neighbourhood(same_classes[i][0]));
            }
         }
      }
#endif
      if (data.type == BranchType::SAME) {
         // update all nodes to have the same neighbourhood.
         // If a same class was previously taken, all nodes in it already have the same neighbourhood

#ifndef NDEBUG
         DenseSet neighbourUnion =
             branched_full_graph.neighbourhood(data.first)
                 .setUnion(branched_full_graph.neighbourhood(data.second));
#endif
         // TODO: below code for adding edges may be somewhat inefficient? Still could be better than setting too many individual bits, however.
         const DenseSet &first_neighbourhood =
             branched_full_graph.neighbourhood(
                 data.first); // TODO: branched_full_graph or original graph?
         const DenseSet &second_neighbourhood =
             branched_full_graph.neighbourhood(data.second);

         DenseSet addToFirst =
             second_neighbourhood.difference(first_neighbourhood);
         DenseSet addToSecond =
             first_neighbourhood.difference(second_neighbourhood);

         // here we do not need to worry about other same classes because all nodes within the same class are both (dis)connected to the two branching nodes
         for (const Node &node : same_classes[data.first]) {
            for (const auto &other_node : addToFirst) {
               branched_full_graph.addEdge(node, other_node);
            }
         }
         for (const Node &node : same_classes[data.second]) {
            for (const auto &other_node : addToSecond) {
               branched_full_graph.addEdge(node, other_node);
            }
         }

         // update the representative and the same class (data.second is 'deleted' later as well!).
         //  We need to do this for all nodes in the second node's same class
         for (const Node &node : same_classes[data.second]) {
            representative[node] = data.first;
         }
         // merge the two same classes into one same class (represented now all by the first node
         same_classes[data.first].insert(same_classes[data.first].end(),
                                         same_classes[data.second].begin(),
                                         same_classes[data.second].end());

#ifndef NDEBUG
         // check if all neighbourhoods are indeed equal after the case
         for (const auto &node : same_classes[data.first]) {
            assert(branched_full_graph.neighbourhood(node) == neighbourUnion);
         }
         for (const auto &node : same_classes[data.second]) {
            assert(branched_full_graph.neighbourhood(node) == neighbourUnion);
         }
#endif
         same_classes[data.second].clear();
      } else if (data.type == BranchType::DIFFER) {
         Node repFirst = representative[data.first];
         Node repSecond = representative[data.second];
         for (const Node &first_node : same_classes[repFirst]) {
            for (const Node &second_node : same_classes[repSecond]) {
               branched_full_graph.addEdge(first_node, second_node);
            }
         }
      }
   }

   assert(branched_full_graph.numSelfLoops() == 0);

   return branched_full_graph;
}
PreprocessingResult
preprocessBranchedGraph(const DenseGraph &t_graph,
                        const std::vector<BranchData> &t_branch_information,
                        std::size_t lower_bound) {
   std::vector<PreprocessedNode> removed_nodes;

   DenseSet present_nodes(t_graph.numNodes(), true);
   assert(present_nodes.full());
   // first; remove all dominated neighbourhoods implied by 'same' nodes
   for (const auto &data : t_branch_information) {
      if (data.type == BranchType::SAME) {
         removed_nodes.emplace_back(
             data.second, PreprocessedReason::DOMINATED_NODE, data.first);
         present_nodes.remove(data.second);
      }
   }
   assert(present_nodes.size() == (present_nodes.capacity()-removed_nodes.size()));

   DenseSet checkForDominatedVertices = present_nodes;
   std::ptrdiff_t lastRoundRemovedNodesSize = 0;

   bool firstRound = true;
   while (true) {
      while (true) {
         auto low_degree_removed = removeLowDegreeVertices(t_graph,present_nodes,lower_bound);
         if (low_degree_removed.empty()) {
            break;
         }
         removed_nodes.insert(removed_nodes.end(), low_degree_removed.begin(),
                              low_degree_removed.end());
      }
      //We only check nodes whoms neighbourhood has changed for if they are dominated
      //In the first iteration we check all present nodes
      if(!firstRound){
         checkForDominatedVertices.clear();
         for(auto node_it = removed_nodes.begin() + lastRoundRemovedNodesSize; node_it != removed_nodes.end(); ++node_it){
            checkForDominatedVertices.inplaceUnion(t_graph.neighbourhood(node_it->removedNode()));
         }
      }else{
         firstRound = false;
      }
      checkForDominatedVertices.inplaceIntersection(present_nodes);
      lastRoundRemovedNodesSize = static_cast<ptrdiff_t>(removed_nodes.size());

      auto dominated_removed =
          removeDominatedVerticesClique(t_graph, present_nodes, checkForDominatedVertices);
      if (dominated_removed.empty()) {
         break;
      }
      removed_nodes.insert(removed_nodes.end(), dominated_removed.begin(),
                           dominated_removed.end());
   }

   auto [newGraph, newToOld] = t_graph.nodeInducedSubgraph(present_nodes);
   NodeMap oldToNew = NodeMap::inverse(newToOld, t_graph.numNodes());

   return {newGraph, PreprocessedMap(removed_nodes, newToOld, oldToNew)};
}
}