//
// Created by rolf on 26-11-22.
//

#include "pcog/utilities/DenseGraph.hpp"

namespace pcog {
DenseGraph::DenseGraph(degree_type t_numNodes)
    : m_adjacencyMatrix(t_numNodes, DenseSet(t_numNodes)) {}

void DenseGraph::addEdge(Node a, Node b) {
   m_adjacencyMatrix[a].add(b);
   m_adjacencyMatrix[b].add(a);
}

bool DenseGraph::isEdge(Node a, Node b) const {
   return m_adjacencyMatrix[a].contains(b);
}

degree_type DenseGraph::numNodes() const { return m_adjacencyMatrix.size(); }

const DenseSet &DenseGraph::neighbourhood(Node t_node) const {
   return m_adjacencyMatrix[t_node];
}

void DenseGraph::removeEdge(Node a, Node b) {
   m_adjacencyMatrix[a].remove(b);
   m_adjacencyMatrix[b].remove(a);
}

bool DenseGraph::nodeDominates(Node a, Node b) const {
   return m_adjacencyMatrix[b].isSubsetOf(m_adjacencyMatrix[a]);
}

bool DenseGraph::nodeStrictlyDominates(Node a, Node b) const {
   return m_adjacencyMatrix[b].isProperSubsetOf(m_adjacencyMatrix[a]);
}

std::pair<DenseGraph, NodeMap>
DenseGraph::nodeInducedSubgraph(const DenseSet &t_set) const {
   assert(t_set.capacity() == numNodes());
   std::size_t new_size = t_set.size();
   DenseGraph graph(new_size);

   auto matrix_it = graph.m_adjacencyMatrix.begin();
   for (const Node &node : t_set) {
      std::size_t index = 0;
      const auto &old_neighboorhood = m_adjacencyMatrix[node];
      for (const Node &neighbourNode : t_set) {
         matrix_it->set(index, old_neighboorhood.contains(neighbourNode));
         index++;
      }
      matrix_it++;
   }

   NodeMap newToOld(new_size);
   Node oldNode = t_set.first();
   for (std::size_t i = 0; i < new_size; i++) {
      newToOld[i] = oldNode;
      oldNode = t_set.find_next(oldNode);
   }

   return {graph, newToOld};
}

degree_type DenseGraph::numEdges() const {
   degree_type sum = 0;
   for (const auto &list : m_adjacencyMatrix) {
      sum += list.size();
   }
   return sum / 2;
}

degree_type DenseGraph::numSelfLoops() const {
   degree_type sum = 0;
   for (Node node = 0; node < m_adjacencyMatrix.size(); node++) {
      if (m_adjacencyMatrix[node].contains(node)) {
         sum++;
      }
   }
   return sum;
}

// void DenseGraph::setNeighbourhood(Node node, const DenseSet& set) {
//    assert(node < adjacency_matrix.size());
//    assert(set.capacity() == adjacency_matrix.size() && set.capacity() ==
//    adjacency_matrix[node].capacity()); adjacency_matrix[node] = set;
// }

bool DenseGraph::setIsStable(const DenseSet &t_set) const {
   Node node = t_set.first();
   DenseSet disallowed_nodes(t_set.capacity());
   while (node != INVALID_NODE) {
      if (disallowed_nodes.contains(node)) {
         return false;
      }
      const auto &neighbours = neighbourhood(node);
      disallowed_nodes.add(node);
      disallowed_nodes.inplaceUnion(neighbours);
      node = t_set.find_next(node);
   }
   return true;
}

void DenseGraph::complement() {
   for (auto &list : m_adjacencyMatrix) {
      list.complement();
   }

   for (Node i = 0; i < m_adjacencyMatrix.size(); ++i) {
      m_adjacencyMatrix[i].remove(i);
   }
}

bool DenseGraph::setIsStableMaximal(const DenseSet &t_set) const {
   DenseSet covered_nodes(t_set.capacity());
   for (const auto &node : t_set) {
      covered_nodes.inplaceUnion(neighbourhood(node));
      covered_nodes.add(node);
   }
   bool result = covered_nodes.full();
   return result;
}

bool DenseGraph::isConsistent() const {
   for (const auto &list : m_adjacencyMatrix) {
      if (list.capacity() != m_adjacencyMatrix.size()) {
         return false;
      }
   }
   for (Node a = 0; a < m_adjacencyMatrix.size(); ++a) {
      for (Node b = 0; b <= a; ++b) {
         if (m_adjacencyMatrix[a].contains(b) !=
             m_adjacencyMatrix[b].contains(a)) {
            return false;
         }
      }
   }
   return true;
}

void DenseGraph::clearEdges() {
   for (auto &list : m_adjacencyMatrix) {
      list.clear();
   }
}

Node DenseGraph::addNodeWithNeighbourhood(DenseSet t_neighbourhood) {
   assert(t_neighbourhood.capacity() == m_adjacencyMatrix.size());
   Node node = m_adjacencyMatrix.size();
   for (Node i = 0; i < t_neighbourhood.capacity(); ++i) {
      m_adjacencyMatrix[i].extend(t_neighbourhood.contains(i));
   }
   t_neighbourhood.extend(false);
   m_adjacencyMatrix.push_back(t_neighbourhood);
   return node;
}
degree_type DenseGraph::nodeDegree(Node t_node) const {
   return neighbourhood(t_node).size();
}
void DenseGraph::clear() {
   m_adjacencyMatrix.clear();
}
void DenseGraph::addEdges(Node node, const DenseSet &other_nodes) {
   m_adjacencyMatrix[node].inplaceUnion(other_nodes);
   for(const auto& other_node : other_nodes){
      m_adjacencyMatrix[other_node].add(node);
   }
}
double DenseGraph::density() const {
   std::size_t nodes = numNodes();
   if(nodes == 0 ) return 0.0;
   double density = static_cast<double>(2*numEdges()) /
                    static_cast<double>(nodes*(nodes-1));
   return density;
}
} // namespace pcog