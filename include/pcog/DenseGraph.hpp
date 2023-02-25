//
// Created by rolf on 26-11-22.
//

#ifndef PCOG_SRC_DENSEGRAPH_HPP
#define PCOG_SRC_DENSEGRAPH_HPP

#include "DenseSet.hpp"
#include "NodeMap.hpp"

namespace pcog {
/// Models an undirected graph. Loops are allowed, but discouraged.
/// Optimized for dense graphs, in particular graphs with more than 1% edge
/// density
class DenseGraph {
 public:
   /// Construct an empty DenseGraph
   DenseGraph() = default;
   /// Construct an empty dense graph with the given number of nodes
   /// \param t_numNodes The number of nodes the graph is initialized in.
   explicit DenseGraph(degree_type t_numNodes);

   //   void setNeighbourhood(Node node, const DenseSet& set);

   /// Adds an edge between node a and node b. If the edge is there already,
   /// nothing happens.
   /// \param a First edge node
   /// \param b Second edge node
   void addEdge(Node a, Node b);

   /// Removes an edge between node a and node b. If the edge was already
   /// removed, nothing happens
   /// \param a First edge node
   /// \param b Second edge node
   void removeEdge(Node a, Node b);
   /// Removes all edges, but keeps the vertices
   void clearEdges();

   /// Removes all edges and vertices, making the graph empty
   void clear();

   /// Complements the graph, removing all edges which were present and adding
   /// all edges which were not. Does not add any self loops
   void complement();
   /// Adds a new node given by the neighbourhood. Neighbourhood must be a dense
   /// set with capacity of numNodes()
   Node addNodeWithNeighbourhood(DenseSet t_neighbourhood);

   /// Returns true if (a,b) is an edge in the graph, false otherwise
   /// \param a The first node
   /// \param b The second ndoe
   /// \return true if (a,b) is an edge in the graph, false otherwise
   [[nodiscard]] bool isEdge(Node a, Node b) const;

   /// \return The number of nodes in the graph.
   [[nodiscard]] degree_type numNodes() const;
   /// Counts the number of edges in the graph. Note this does NOT count self
   /// loops correctly! \return The number of edges in the graph
   [[nodiscard]] degree_type numEdges() const;
   /// Gets the dense set of the neighbourhood of the given node
   /// \param t_node The node to get the neighbourhood for
   /// \return The neighbourhood of the node
   [[nodiscard]] const DenseSet &neighbourhood(Node t_node) const;
   /// Returns the degree of the given node, the number of adjacent adges
   /// \param t_node The node to return the degree for
   /// \return The degree of the given node
   [[nodiscard]] degree_type nodeDegree(Node t_node) const;
   /// Counts the number of nodes which has a selfloop (e.g. an (a,a) edge for
   /// node a) \return The number of self loops
   [[nodiscard]] degree_type numSelfLoops() const;

   /// Checks if node a dominates node b
   /// \param a The first node
   /// \param b The second node
   /// \return Returns true if the neighbourhood of b is a subset of the
   /// neighbourhood of a
   [[nodiscard]] bool nodeDominates(Node a, Node b) const;
   /// Checks if node a dominates node b, i.e. N(b) is a proper subset of N(a)
   /// \param a The first node
   /// \param b The second node
   /// \return Returns true if the neighbourhood of b is a proper subset of the
   /// neighbourhood of a
   [[nodiscard]] bool nodeStrictlyDominates(Node a, Node b) const;
   /// Checks if the given set of nodes is stable in the graph
   /// \param t_set The set to check
   /// \return True if the set is stable, false otherwise
   [[nodiscard]] bool setIsStable(const DenseSet &t_set) const;
   /// Checks if the given set of nodes is maximally stable
   /// e.g. if no node v can be added such that v is independent of all nodes in
   /// t_set \param t_set The set to check \return True if the set is maximally
   /// stable.
   [[nodiscard]] bool setIsStableMaximal(const DenseSet &t_set) const;

   /// Computes the induced graph of the current graph based on a subset of
   /// vertices. The induced graph is formed by all vertices of the given subset
   /// and any edges between nodes of this subset.
   /// \param t_set The (sub)set of vertices for which to compute the induced
   /// graph.
   /// \return A pair containing the induced graph and a NodeMap which
   /// maps the nodes from the new graph to the old graphs' nodes.
   [[nodiscard]] std::pair<DenseGraph, NodeMap>
   nodeInducedSubgraph(const DenseSet &t_set) const;

   /// Checks if the graph data structures are internally consistent.
   /// In particular as we store edges twice it checks if both or neither are
   /// present.
   /// \return True if the data structures are in a valid state
   [[nodiscard]] bool isConsistent() const;

 private:
   std::vector<DenseSet> m_adjacencyMatrix;
};
} // namespace pcog
#endif // PCOG_SRC_DENSEGRAPH_HPP
