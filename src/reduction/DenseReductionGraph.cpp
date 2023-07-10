//
// Created by rolf on 4-7-23.
//

#include "pcog/reduction/DenseReductionGraph.hpp"
namespace pcog{
DenseReductionGraph::DenseReductionGraph(const DenseGraph &graph)
: adjacencyMatrix(graph.numNodes()), presentNodes(graph.numNodes(),true){
   degrees.resize(graph.numNodes());
   for(std::size_t i = 0; i < graph.numNodes(); ++i){
      adjacencyMatrix[i] = graph.neighbourhood(i);
      degrees[i] = graph.neighbourhood(i).size();
   }
}
void DenseReductionGraph::removeGraphNode(Node node) {
   assert(presentNodes.contains(node));
   presentNodes.remove(node);
   for(Node neighbour : adjacencyMatrix[node]){
      --degrees[neighbour];
      assert(degrees[neighbour] < presentNodes.size());
      adjacencyMatrix[neighbour].remove(node);
   }
   degrees[node] = 0;
   adjacencyMatrix[node].clear();
}
const DenseSet &DenseReductionGraph::neighbourhood(Node node) const {
   return adjacencyMatrix[node];
}
std::size_t DenseReductionGraph::nodeDegree(Node node) const {
   return degrees[node];
}
bool DenseReductionGraph::setIsStable(const DenseSet &set) const {
   assert(set.intersection(presentNodes).size() == set.size());
   Node node = set.first();
   DenseSet disallowed_nodes(set.capacity());
   while (node != INVALID_NODE) {
      if (disallowed_nodes.contains(node)) {
         return false;
      }
      const auto &neighbours = neighbourhood(node);
      disallowed_nodes.add(node);
      disallowed_nodes.inplaceUnion(neighbours);
      node = set.find_next(node);
   }
   return true;
}
bool DenseReductionGraph::setIsStableMaximal(const DenseSet &set) const {
   DenseSet covered_nodes(set.capacity());
   for (const auto &node : set) {
      covered_nodes.inplaceUnion(neighbourhood(node));
      covered_nodes.add(node);
   }
   bool result = covered_nodes.full();
   return result;
}
const DenseSet &DenseReductionGraph::nodes() const { return presentNodes; }
bool DenseReductionGraph::containsNode(Node node) const {
   return presentNodes.contains(node);
}
bool DenseReductionGraph::hasLowerBound() const {
   return cliqueNodes.capacity() != 0;
}
std::size_t DenseReductionGraph::lowerBound() const {
   return lowerbound;
}
const DenseSet &DenseReductionGraph::lowerBoundNodes() const {
   return cliqueNodes;
}
void DenseReductionGraph::setLowerBoundClique(const DenseSet &clique,
                                              std::size_t numCliqueNodes ) {
   assert(clique.size() == numCliqueNodes);
   assert(clique.intersection(presentNodes).size() == numCliqueNodes);
   cliqueNodes = clique;
   lowerbound = numCliqueNodes;
}
void DenseReductionGraph::removeNode(Node node) {
   removeGraphNode(node);
   if(hasLowerBound() && lowerBoundNodes().contains(node)){
      --lowerbound;
   }
}
void DenseReductionGraph::removeStableSet(const DenseSet &set) {
   for(const auto& node : set){
      removeGraphNode(node);
   }
   if(hasLowerBound() && !lowerBoundNodes().intersection(set).empty()){
      --lowerbound;
   }
}

}