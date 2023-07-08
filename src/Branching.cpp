//
// Created by rolf on 29-11-22.
//

#include "pcog/Branching.hpp"
#include <numeric>

#include <iostream> //TODO: remove again
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
         // TODO: below code for adding edges could be more efficient
         const DenseSet &first_neighbourhood =
             branched_full_graph.neighbourhood(data.first);
         const DenseSet &second_neighbourhood =
             branched_full_graph.neighbourhood(data.second);

         DenseSet addToFirst =
             second_neighbourhood.difference(first_neighbourhood);
         DenseSet addToSecond =
             first_neighbourhood.difference(second_neighbourhood);

         // here we do not need to worry about other same classes because all nodes within the same class are both (dis)connected to the two branching nodes
         //TODO: benchmark
         for(const Node& node : same_classes[data.first]){
            branched_full_graph.addEdges(node,addToFirst);
         }
         for(const Node& node : same_classes[data.second]){
            branched_full_graph.addEdges(node,addToSecond);
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

DenseGraph
constructBranchedFullGraphFromChild(const DenseGraph& t_childGraph,
                                    const std::vector<BranchData>& branch_information,
                                    std::size_t numAddedDecisions){
   std::size_t startChildIndex = branch_information.size()-numAddedDecisions;
   DenseGraph branched_full_graph(t_childGraph);
   std::size_t num_nodes = t_childGraph.numNodes();

   std::vector<Node> representative(num_nodes);
   std::iota(representative.begin(), representative.end(), 0);
   std::vector<std::vector<Node>> same_classes(num_nodes, {0});
   for (std::size_t i = 0; i < num_nodes; i++) {
      same_classes[i][0] = i;
   }
   for(std::size_t i = 0; i < startChildIndex; i++){
      const auto& data = branch_information[i];
      if (data.type == BranchType::SAME) {

         // update the representative and the same class (data.second is 'deleted' later as well!).
         //  We need to do this for all nodes in the second node's same class
         for (const Node &node : same_classes[data.second]) {
            representative[node] = data.first;
         }
         // merge the two same classes into one same class (represented now all by the first node
         same_classes[data.first].insert(same_classes[data.first].end(),
                                         same_classes[data.second].begin(),
                                         same_classes[data.second].end());
         same_classes[data.second].clear();
         //TODO: check that all edges have the same neighbourhood
      }
#ifndef NDEBUG
      if (data.type == BranchType::DIFFER) {
         Node repFirst = representative[data.first];
         Node repSecond = representative[data.second];
         for (const Node &first_node : same_classes[repFirst]) {
            for (const Node &second_node : same_classes[repSecond]) {
               assert(branched_full_graph.isEdge(first_node,second_node));
            }
         }
      }
#endif

   }
   //Now actually add edges rather than checking + constructing same classes
   for(std::size_t i = startChildIndex; i < branch_information.size(); ++i){
      const auto& data = branch_information[i];
      assert(data.type == BranchType::ROOT ||
             !t_childGraph.isEdge(data.first, data.second));
      assert(data.type == BranchType::ROOT ||
             !t_childGraph.isEdge(data.second, data.first));
      assert(data.type == BranchType::ROOT ||
             representative[data.first] == data.first);
      assert(data.type == BranchType::ROOT ||
             representative[data.second] == data.second);

#ifndef NDEBUG
      // Check if all same classes are still correct
      for (std::size_t j = 0; j < num_nodes; ++j) {
         if (!same_classes[j].empty()) {
            for (const auto &node : same_classes[j]) {
               assert(branched_full_graph.neighbourhood(node) ==
                      branched_full_graph.neighbourhood(same_classes[j][0]));
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
         // TODO: below code for adding edges could be more efficient
         const DenseSet &first_neighbourhood =
             branched_full_graph.neighbourhood(data.first);
         const DenseSet &second_neighbourhood =
             branched_full_graph.neighbourhood(data.second);

         DenseSet addToFirst =
             second_neighbourhood.difference(first_neighbourhood);
         DenseSet addToSecond =
             first_neighbourhood.difference(second_neighbourhood);


         // here we do not need to worry about other same classes because all nodes within the same class are both (dis)connected to the two branching nodes
         //TODO: benchmark
         for(const Node& node : same_classes[data.first]){
            branched_full_graph.addEdges(node,addToFirst);
         }
         for(const Node& node : same_classes[data.second]){
            branched_full_graph.addEdges(node,addToSecond);
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
   std::vector<DenseSet> fixed_sets;

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

   std::size_t lastSuccessfulMethod = 0;
   std::size_t round = 0;
   //TODO: Fix; below loop does not preprocess everything.
   while (true) {
      {
         bool lowDegreeSuccesful = false;
         while (true) {
            auto low_degree_removed = removeLowDegreeVertices(t_graph,present_nodes,lower_bound);
            if (low_degree_removed.empty()) {
               break;
            }
            lowDegreeSuccesful = true;
            removed_nodes.insert(removed_nodes.end(),
                                 low_degree_removed.begin(),
                                 low_degree_removed.end());
         }
         if (lowDegreeSuccesful) {
            lastSuccessfulMethod = 0;
         } else if (lastSuccessfulMethod == 1) {
            break;
         }
      }
      {
         auto roundFixedSets = fixStableNeighbourhoods(t_graph,present_nodes);
         for(const auto& set : roundFixedSets){
            for(Node node : set){
               removed_nodes.emplace_back(node,PreprocessedReason::STABLE_NEIGHBOURHOOD,fixed_sets.size());
            }
            fixed_sets.push_back(set);
         }
         assert(lower_bound >= roundFixedSets.size());
         lower_bound -= roundFixedSets.size();
         if(!roundFixedSets.empty()){
            lastSuccessfulMethod = 1;
         }else if(removed_nodes.size() == t_graph.numNodes() ||
                    lastSuccessfulMethod == 2){
            break;
         }
      }

      //We only check nodes whoms neighbourhood has changed if they are removed or not.
      //In the first iteration we check all present nodes
      {
         if (round != 0) {
            checkForDominatedVertices.clear();
            for (auto node_it =
                     removed_nodes.begin() + lastRoundRemovedNodesSize;
                 node_it != removed_nodes.end(); ++node_it) {
               checkForDominatedVertices.inplaceUnion(
                   t_graph.neighbourhood(node_it->removedNode()));
            }
         }
         checkForDominatedVertices.inplaceIntersection(present_nodes);

         lastRoundRemovedNodesSize =
             static_cast<ptrdiff_t>(removed_nodes.size());

         bool dominatedSuccesful = false;
         std::size_t dominated_round = 0;
         while(true){
            auto dominated_removed = removeDominatedVerticesClique(
                t_graph, present_nodes, checkForDominatedVertices);
            if(dominated_removed.empty()){
               break;
            }
            dominatedSuccesful = true;

            removed_nodes.insert(removed_nodes.end(), dominated_removed.begin(),
                                 dominated_removed.end());
            ++dominated_round;
            checkForDominatedVertices.clear();
            for(const auto& node : dominated_removed){
               checkForDominatedVertices.inplaceUnion(t_graph.neighbourhood(node.removedNode()));
            }
            checkForDominatedVertices.inplaceIntersection(present_nodes);
         }

         if (dominatedSuccesful) {
            lastSuccessfulMethod = 2;
         } else if (lastSuccessfulMethod == 0) {
            break;
         }
      }
      ++round;
   }
//   bool firstRound = true;
//   while (true) {
//      while (true) {
//         auto low_degree_removed = removeLowDegreeVertices(t_graph,present_nodes,lower_bound);
//         if (low_degree_removed.empty()) {
//            break;
//         }
//         removed_nodes.insert(removed_nodes.end(), low_degree_removed.begin(),
//                              low_degree_removed.end());
//      }
//      //We only check nodes whoms neighbourhood has changed for if they are dominated
//      //In the first iteration we check all present nodes
//      if(!firstRound){
//         checkForDominatedVertices.clear();
//         for(auto node_it = removed_nodes.begin() + lastRoundRemovedNodesSize; node_it != removed_nodes.end(); ++node_it){
//            checkForDominatedVertices.inplaceUnion(t_graph.neighbourhood(node_it->removedNode()));
//         }
//      }else{
//         firstRound = false;
//      }
//      checkForDominatedVertices.inplaceIntersection(present_nodes);
//      lastRoundRemovedNodesSize = static_cast<ptrdiff_t>(removed_nodes.size());
//
//      auto dominated_removed =
//          removeDominatedVerticesClique(t_graph, present_nodes, checkForDominatedVertices);
//      if (dominated_removed.empty()) {
//         break;
//      }
//      removed_nodes.insert(removed_nodes.end(), dominated_removed.begin(),
//                           dominated_removed.end());
//   }
#ifndef NDEBUG
   //ensure that there are no more dominated nodes in the graph
   for(Node node : present_nodes){
      for(Node otherNode : present_nodes){
         if(otherNode == node){
            continue;
         }
         DenseSet nb1 = t_graph.neighbourhood(node).intersection(present_nodes);
         DenseSet nb2 = t_graph.neighbourhood(otherNode).intersection(present_nodes);
         assert(!nb1.isProperSubsetOf(nb2));
      }
   }
#endif
   auto [newGraph, newToOld] = t_graph.nodeInducedSubgraph(present_nodes);
   NodeMap oldToNew = NodeMap::inverse(newToOld, t_graph.numNodes());

   return {newGraph, PreprocessedMap(removed_nodes,fixed_sets, newToOld, oldToNew)};
}


PreprocessingResult preprocessedGraphFromChild(const DenseGraph& t_preprocessedChildGraph,
                                               const std::vector<BranchData>& t_branch_information,
                                               std::size_t numNewBranchingConstraints,
                                               const NodeMap& preprocessedToChildMap,
                                               std::size_t lower_bound){

   DenseGraph updatedGraph(t_preprocessedChildGraph);
   //All output is in terms of the current node index: after this function,
   //you will typically want to combine this information with the current
   std::vector<PreprocessedNode> removed_nodes;
   DenseSet present_nodes(updatedGraph.numNodes(), true);
   // first; remove all dominated neighbourhoods implied by 'same' nodes
   std::size_t startIndex = t_branch_information.size()- numNewBranchingConstraints;

   //TODO: should keep track of SAME classes here when adding more than one node
   assert(t_branch_information.size() - startIndex <= 1); //For now, to ensure everything works okay. Later remove after above is fixed

   for(std::size_t i = startIndex; i < t_branch_information.size(); ++i){
      const auto& data = t_branch_information[i];
      Node childFirst = preprocessedToChildMap[data.first];
      Node childSecond = preprocessedToChildMap[data.second];
      assert(childFirst != INVALID_NODE && childSecond != INVALID_NODE);
      if (data.type == BranchType::SAME) {
         DenseSet addToSecond = updatedGraph.neighbourhood(childFirst).difference(updatedGraph.neighbourhood(childSecond));
         DenseSet addToFirst = updatedGraph.neighbourhood(childSecond).difference(updatedGraph.neighbourhood(childFirst));
         updatedGraph.addEdges(childFirst,addToFirst);
         updatedGraph.addEdges(childFirst,addToSecond);

         removed_nodes.emplace_back(childSecond, PreprocessedReason::DOMINATED_NODE, childFirst);
         present_nodes.remove(childSecond);
      }else if (data.type == BranchType::DIFFER){
         updatedGraph.addEdge(childFirst,childSecond);
      }
   }
   assert(present_nodes.size() == (present_nodes.capacity()-removed_nodes.size()));

   DenseSet checkForDominatedVertices = present_nodes;
   std::ptrdiff_t lastRoundRemovedNodesSize = 0;

   std::vector<DenseSet> fixed_sets;

   std::size_t lastSuccessfulMethod = 0;
   std::size_t round = 0;
   //TODO: Fix; below loop does not preprocess everything.
   while (true) {
      {
         bool lowDegreeSuccesful = false;
         while (true) {
            auto low_degree_removed = removeLowDegreeVertices(updatedGraph,present_nodes,lower_bound);
            if (low_degree_removed.empty()) {
               break;
            }
            lowDegreeSuccesful = true;
            removed_nodes.insert(removed_nodes.end(),
                                 low_degree_removed.begin(),
                                 low_degree_removed.end());
         }
         if (lowDegreeSuccesful) {
            lastSuccessfulMethod = 0;
         } else if (lastSuccessfulMethod == 1) {
            break;
         }
      }
      {
         auto roundFixedSets = fixStableNeighbourhoods(updatedGraph,present_nodes);
         for(const auto& set : roundFixedSets){
            for(Node node : set){
               removed_nodes.emplace_back(node,PreprocessedReason::STABLE_NEIGHBOURHOOD,fixed_sets.size());
            }
            fixed_sets.push_back(set);
         }
         assert(lower_bound >= roundFixedSets.size());
         lower_bound -= roundFixedSets.size();
         if(!roundFixedSets.empty()){
            lastSuccessfulMethod = 1;
         }else if(removed_nodes.size() == updatedGraph.numNodes() ||
                    lastSuccessfulMethod == 2){
            break;
         }
      }

      //We only check nodes whoms neighbourhood has changed if they are removed or not.
      //In the first iteration we check all present nodes
      {
         if (round != 0) {
            checkForDominatedVertices.clear();
            for (auto node_it =
                     removed_nodes.begin() + lastRoundRemovedNodesSize;
                 node_it != removed_nodes.end(); ++node_it) {
               checkForDominatedVertices.inplaceUnion(
                   updatedGraph.neighbourhood(node_it->removedNode()));
            }
         }
         checkForDominatedVertices.inplaceIntersection(present_nodes);

         lastRoundRemovedNodesSize =
             static_cast<ptrdiff_t>(removed_nodes.size());

         bool dominatedSuccesful = false;
         std::size_t dominated_round = 0;
         while(true){
            auto dominated_removed = removeDominatedVerticesClique(
                updatedGraph, present_nodes, checkForDominatedVertices);
            if(dominated_removed.empty()){
               break;
            }
            dominatedSuccesful = true;

            removed_nodes.insert(removed_nodes.end(), dominated_removed.begin(),
                                 dominated_removed.end());
            ++dominated_round;
            checkForDominatedVertices.clear();
            for(const auto& node : dominated_removed){
               checkForDominatedVertices.inplaceUnion(updatedGraph.neighbourhood(node.removedNode()));
            }
            checkForDominatedVertices.inplaceIntersection(present_nodes);
         }

         if (dominatedSuccesful) {
            lastSuccessfulMethod = 2;
         } else if (lastSuccessfulMethod == 0) {
            break;
         }
      }
      ++round;
   }

//   bool firstRound = true;
//   while (true) {
//      while (true) {
//         auto low_degree_removed = removeLowDegreeVertices(updatedGraph,present_nodes,lower_bound);
//         if (low_degree_removed.empty()) {
//            break;
//         }
//         removed_nodes.insert(removed_nodes.end(), low_degree_removed.begin(),
//                              low_degree_removed.end());
//      }
//      //We only check nodes whoms neighbourhood has changed for if they are dominated
//      //In the first iteration we check all present nodes
//      if(!firstRound){
//         checkForDominatedVertices.clear();
//         for(auto node_it = removed_nodes.begin() + lastRoundRemovedNodesSize; node_it != removed_nodes.end(); ++node_it){
//            checkForDominatedVertices.inplaceUnion(updatedGraph.neighbourhood(node_it->removedNode()));
//         }
//      }else{
//         firstRound = false;
//      }
//      checkForDominatedVertices.inplaceIntersection(present_nodes);
//      lastRoundRemovedNodesSize = static_cast<ptrdiff_t>(removed_nodes.size());
//
//      auto dominated_removed =
//          removeDominatedVerticesClique(updatedGraph, present_nodes, checkForDominatedVertices);
//      if (dominated_removed.empty()) {
//         break;
//      }
//      removed_nodes.insert(removed_nodes.end(), dominated_removed.begin(),
//                           dominated_removed.end());
//   }

#ifndef NDEBUG
   //ensure that there are no more dominated nodes in the graph
   for(Node node : present_nodes){
      for(Node otherNode : present_nodes){
         if(otherNode == node){
            continue;
         }
         DenseSet nb1 = updatedGraph.neighbourhood(node).intersection(present_nodes);
         DenseSet nb2 = updatedGraph.neighbourhood(otherNode).intersection(present_nodes);
         assert(!nb1.isProperSubsetOf(nb2));
      }
   }
#endif
   auto [newGraph, newToOld] = updatedGraph.nodeInducedSubgraph(present_nodes);
   NodeMap oldToNew = NodeMap::inverse(newToOld, updatedGraph.numNodes());

   return {newGraph, PreprocessedMap(removed_nodes,fixed_sets, newToOld, oldToNew)};
}

}