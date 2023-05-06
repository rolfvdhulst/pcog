//
// Created by rolf on 27-11-22.
//

#include "pcog/Preprocessing.hpp"
#include "pcog/mwss/CombinatorialStableSet.hpp"
#include <utility>

namespace pcog {
PreprocessingResult::PreprocessingResult(DenseGraph graph, PreprocessedMap map)
    : graph{std::move(graph)}, map{std::move(map)} {}

DenseSet findInitialClique(const DenseGraph& t_graph){
   DenseSet data;

   auto mss_callback = [](const DenseSet& current_nodes, MaxStableSetCombinatorial::weight_type total_weight, void * user_data, bool first_solution, bool& stop_solving, bool& accepted_solution){
      auto * data = static_cast<DenseSet *>(user_data);
      *data = current_nodes;
      accepted_solution = true;
      stop_solving = false;
   };
   DenseGraph complementGraph = t_graph;
   complementGraph.complement();

   MaxStableSetCombinatorial mss(UniformWeightFunction(),complementGraph);
   mss.setUserData(&data);
   mss.setCallback(mss_callback);
   mss.setInfiniteUpperBound();
   mss.setNodeLimit(t_graph.numNodes()*10);

   mss.run();

   std::cout<<mss.numBranchAndBoundNodes()<<" b&b nodes for finding clique  ns/node: "<<double(mss.timeTaken().count()) / mss.numBranchAndBoundNodes() <<std::endl;
   std::cout<<"total time: "<<double(mss.timeTaken().count())/1e9<<" seconds"<<std::endl;
   std::cout<<"found clique size: "<<data.size()<<std::endl;

   return data;
}
PreprocessingResult preprocessOriginalGraph(const DenseGraph &t_graph) {
   std::vector<PreprocessedNode> removed_nodes;
   std::vector<DenseSet> fixed_sets;

   DenseSet present_nodes(t_graph.numNodes(), true);
   assert(present_nodes.full());

   DenseSet clique = findInitialClique(t_graph);

   DenseSet checkForDominatedVertices = present_nodes;
   std::ptrdiff_t lastRoundRemovedNodesSize = 0;

   std::size_t lastSuccessfulMethod = 0;
   std::size_t round = 0;
   while (true) {
      auto roundFixedSets = fixStableNeighbourhoods(t_graph,present_nodes);
      for(const auto& set : roundFixedSets){
         for(Node node : set){
            removed_nodes.emplace_back(node,PreprocessedReason::STABLE_NEIGHBOURHOOD,fixed_sets.size());
         }
         fixed_sets.push_back(set);
      }
      if(!roundFixedSets.empty()){
         lastSuccessfulMethod = 0;
      }else if(removed_nodes.size() == t_graph.numNodes() ||
                 lastSuccessfulMethod == 1){
         break;
      }
      {
         bool lowDegreeSuccesful = false;
         while (true) {
            auto low_degree_removed =
                removeLowDegreeVerticesClique(t_graph, present_nodes, clique);
            if (low_degree_removed.empty()) {
               break;
            }
            lowDegreeSuccesful = true;
            removed_nodes.insert(removed_nodes.end(), low_degree_removed.begin(),
                                 low_degree_removed.end());
         }
         if(lowDegreeSuccesful){
            lastSuccessfulMethod = 1;
         }else if(lastSuccessfulMethod == 2){
            break;
         }
      }

      //We only check nodes whoms neighbourhood has changed if they are removed or not.
      //In the first iteration we check all present nodes
      if(round != 0){
         checkForDominatedVertices.clear();
         for(auto node_it = removed_nodes.begin() + lastRoundRemovedNodesSize; node_it != removed_nodes.end(); ++node_it){
            checkForDominatedVertices.inplaceUnion(t_graph.neighbourhood(node_it->removedNode()));
         }
      }
      checkForDominatedVertices.inplaceIntersection(present_nodes);
      checkForDominatedVertices.inplaceDifference(clique);
      lastRoundRemovedNodesSize = static_cast<ptrdiff_t>(removed_nodes.size());

      auto dominated_removed =
          removeDominatedVerticesClique(t_graph, present_nodes, checkForDominatedVertices);
      if(!dominated_removed.empty()){
         lastSuccessfulMethod = 2;
      }else if(lastSuccessfulMethod == 0){
         break;
      }
      removed_nodes.insert(removed_nodes.end(), dominated_removed.begin(),
                           dominated_removed.end());

      ++round;
   }

   auto [newGraph, newToOld] = t_graph.nodeInducedSubgraph(present_nodes);
   NodeMap oldToNew = NodeMap::inverse(newToOld, t_graph.numNodes());

   std::cout<<"Preprocessing removed "<< removed_nodes.size()<<" nodes, fixed "<<fixed_sets.size()<<" sets\n";
   return {newGraph, PreprocessedMap(removed_nodes, fixed_sets, newToOld, oldToNew)};
}

std::vector<PreprocessedNode>
removeLowDegreeVerticesClique(const DenseGraph &graph, DenseSet &present_nodes,
                              const DenseSet &clique) {
   std::size_t lower_bound = clique.size();
   DenseSet iterate_nodes = present_nodes;
   iterate_nodes.inplaceDifference(clique);

   DenseSet removed_nodes(graph.numNodes());
   assert(removed_nodes.empty());

   for (const auto &node : iterate_nodes) {
      if (graph.neighbourhood(node).intersection(present_nodes).size() <
          lower_bound) {
         removed_nodes.add(node);
      }
   }

   present_nodes.inplaceDifference(removed_nodes);
   std::vector<PreprocessedNode> removed;
   for (const auto &node : removed_nodes) {
      removed.emplace_back(
          PreprocessedNode(node, PreprocessedReason::LOW_DEGREE));
   }
   return removed;
}

std::vector<PreprocessedNode>
removeDominatedVerticesClique(const DenseGraph &graph, DenseSet &present_nodes,
                              const DenseSet &checkForDominatedVertices) {

   std::vector<PreprocessedNode> removed_nodes;
   for(const auto& node_a : checkForDominatedVertices){
      const DenseSet masked_neighbourhood_a = graph.neighbourhood(node_a).intersection(present_nodes);
      Node neighbourhood_node = masked_neighbourhood_a.first();
      //TODO: it might be worth to heuristically select a node with a small neighbourhood as initial point to speed up the termination of below loop,
      //particularly for larger graphs
      DenseSet possible_dominating_nodes = graph.neighbourhood(neighbourhood_node).intersection(present_nodes);
      possible_dominating_nodes.remove(node_a);

      while(neighbourhood_node != INVALID_NODE){
         if(possible_dominating_nodes.empty()){
            break;
         }
         possible_dominating_nodes.inplaceIntersection(graph.neighbourhood(neighbourhood_node));
         neighbourhood_node = masked_neighbourhood_a.find_next(neighbourhood_node);
      }
      if(!possible_dominating_nodes.empty()){
         removed_nodes.emplace_back(node_a, PreprocessedReason::DOMINATED_NODE, possible_dominating_nodes.first());
         present_nodes.remove(node_a);
      }
   }
   return removed_nodes;

}

PreprocessedNode::PreprocessedNode(Node t_node, PreprocessedReason t_reason,
                                   Node t_dominatedByNode)
    : m_node{t_node}, m_reason{t_reason}, m_dominated_by{t_dominatedByNode} {}

Node PreprocessedNode::removedNode() const { return m_node; }

PreprocessedReason PreprocessedNode::removedReason() const { return m_reason; }

Node PreprocessedNode::dominatingNode() const { return m_dominated_by; }

PreprocessedMap::PreprocessedMap(std::vector<PreprocessedNode> removed_nodes,
                                 std::vector<DenseSet> t_fixedSets,
                                 NodeMap newToOld, NodeMap oldToNew)
    : removed_nodes{std::move(removed_nodes)},
      fixed_sets{std::move(t_fixedSets)},
      newToOldIDs{std::move(newToOld)},
      oldToNewIDs{std::move(oldToNew)} {}

void PreprocessedMap::clear() {
   removed_nodes.clear();
   newToOldIDs.clear();
   oldToNewIDs.clear();
}
void PreprocessedMap::extend(const PreprocessedMap &t_map) {
   //change the id's of the map which extends this
    std::vector<PreprocessedNode> newly_removed_nodes = t_map.removed_nodes;
    for(auto& node : newly_removed_nodes){
      std::size_t oldRemoved = newToOldIDs[node.removedNode()];
      std::size_t oldDominated = node.dominatingNode() == INVALID_NODE ? INVALID_NODE : newToOldIDs[node.dominatingNode()];
      node = PreprocessedNode(oldRemoved,node.removedReason(),oldDominated);
    }
    removed_nodes.insert(removed_nodes.end(),
                         newly_removed_nodes.begin(),newly_removed_nodes.end());
    //update the old->new and new->old mappings
    NodeMap newToOld = NodeMap::identity(t_map.newToOldIDs.size());
    for(std::size_t i = 0; i < t_map.newToOldIDs.size(); ++i){
      newToOld[i] = newToOldIDs[t_map.newToOldIDs[i]];
    }
    newToOldIDs = newToOld;

    for(std::size_t i = 0; i < oldToNewIDs.size(); ++i){
      std::size_t childID = oldToNewIDs[i];
      if(childID != INVALID_NODE){
         oldToNewIDs[i] = t_map.oldToNewIDs[oldToNewIDs[i]];
      }
    }
}

std::vector<PreprocessedNode> removeLowDegreeVertices(const DenseGraph &graph,
                                                      DenseSet &present_nodes,
                                                      std::size_t lower_bound) {
   if(lower_bound == 0){
      return {};
   }
   DenseSet removed_nodes(graph.numNodes());
   assert(removed_nodes.empty());

   for (const auto &node : present_nodes) {
      // Note the -1, otherwise we could remove clique nodes from the clique
      // If we have a subgraph which we know for certain attains this bound,
      // we can safely remove nodes outside of it
      if (graph.neighbourhood(node).intersection(present_nodes).size() < lower_bound - 1) {
         removed_nodes.add(node);
      }
   }

   present_nodes.inplaceDifference(removed_nodes);
   std::vector<PreprocessedNode> removed;
   for (const auto &node : removed_nodes) {
      removed.emplace_back(node, PreprocessedReason::LOW_DEGREE);
   }
   return removed;
}

std::vector<PreprocessedNode>
removeStrictlyDominatedVertices(const DenseGraph &graph,
                                DenseSet &present_nodes) {

   std::vector<PreprocessedNode> removed_nodes;
   for (const auto &node_a : present_nodes) {
      for (const auto &node_b : present_nodes) {
         if (node_a == node_b) {
            continue;
         }
         const DenseSet masked_neighbourhood_a =
             graph.neighbourhood(node_a).intersection(present_nodes);
         const DenseSet masked_neighbourhood_b =
             graph.neighbourhood(node_b).intersection(present_nodes);
         if (masked_neighbourhood_a.isProperSubsetOf(masked_neighbourhood_b)) {
            removed_nodes.emplace_back(
                node_a, PreprocessedReason::DOMINATED_NODE, node_b);
            break;
         }
      }
   }
   for (const auto &processed_node : removed_nodes) {
      present_nodes.remove(processed_node.removedNode());
   }
   return removed_nodes;
}

NodeColoring extendColoring(const NodeColoring &coloring,
                            const PreprocessedMap &map,
                            const DenseGraph &originalGraph) {
   NodeColoring newColoring(map.oldToNewIDs.size());
   newColoring.setNumColors(coloring.numColors());

   for (std::size_t node = 0; node < coloring.numNodes(); ++node) {
      newColoring[map.newToOldIDs[node]] = coloring[node];
   }
   if (!map.removed_nodes.empty()) {
      DenseSet unused_colors(coloring.numColors());
      for (auto it = map.removed_nodes.rbegin(); it != map.removed_nodes.rend();
           it++) {
         const PreprocessedNode &preprocessedNode = *it;
         if (preprocessedNode.removedReason() ==
             PreprocessedReason::LOW_DEGREE) {
            assert(newColoring[preprocessedNode.removedNode()] >=
                   coloring.numColors());
            unused_colors.setAll();
            Node removed_node = preprocessedNode.removedNode();
            for (const auto &node : originalGraph.neighbourhood(removed_node)) {
               if (newColoring[node] <
                   coloring.numColors()) { // the removed node may have some
                                           // other removed
                  // nodes in its neighbourhood which were
                  // removed earlier
                  unused_colors.remove(newColoring[node]);
               }
            }
            std::size_t assigned_color = unused_colors.first();
            newColoring[removed_node] = assigned_color;
            assert(newColoring[preprocessedNode.removedNode()] <
                   coloring.numColors());
         } else if (preprocessedNode.removedReason() ==
                    PreprocessedReason::DOMINATED_NODE) {
            // in case the node was dominated, we can extend it by using the
            // same color
            assert(newColoring[preprocessedNode.removedNode()] >=
                   coloring.numColors());
            assert(newColoring[preprocessedNode.dominatingNode()] <
                   coloring.numColors());
            newColoring[preprocessedNode.removedNode()] =
                newColoring[preprocessedNode.dominatingNode()];
            assert(newColoring[preprocessedNode.removedNode()] <
                   coloring.numColors());
         }
      }
   }

   assert(newColoring.hasNoInvalidNodes());
   return newColoring;
}

std::vector<DenseSet>
fixStableNeighbourhoods(const DenseGraph &graph, DenseSet &present_nodes){
   std::vector<DenseSet> fixedSets;
   bool found;
   do {
      found = false;
      for (Node node : present_nodes) {
         DenseSet set = graph.neighbourhood(node);
         set.complement();
         set.inplaceIntersection(present_nodes);

         bool setIsStable = true;
         DenseSet disallowed_nodes(set.capacity());
         for(const auto& setNode : set){
            if(disallowed_nodes.contains(setNode)){
               setIsStable = false;
               break;
            }
            disallowed_nodes.add(setNode);
            disallowed_nodes.inplaceUnion(graph.neighbourhood(setNode)); //no need to constrain to intersection here
         }
         if(setIsStable){
            found = true;
            fixedSets.push_back(set);
            present_nodes.inplaceDifference(set);
            break;
         }
      }
   } while(found);

   return fixedSets;
}
} // namespace pcog