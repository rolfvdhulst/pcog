//
// Created by rolf on 7-4-23.
//

#ifndef PCOG_SOLUTIONDATA_HPP
#define PCOG_SOLUTIONDATA_HPP

#include <utility>

#include "BBNode.hpp"
#include "pcog/utilities/DenseGraph.hpp"
#include "pcog/utilities/DenseSet.hpp"
#include "Preprocessing.hpp"
#include "Settings.hpp"

namespace pcog {
class StableSetVariable{
 public:
   explicit StableSetVariable(DenseSet set) : m_set{std::move(set)}{};
   [[nodiscard]] const DenseSet& set() const { return m_set;}
 private:
   DenseSet m_set;
};

/// Class which holds the data produced during a solve of the graph.
/// This data can have frequent reads and writes from multiple threads
class SolutionData {
 public:
   explicit SolutionData( Settings& t_settings);

   /// This data is only changed during presolve, not during branch-and-bound
   [[nodiscard]] const DenseGraph& originalGraph() const;
   [[nodiscard]] const DenseGraph& preprocessedGraph() const;
   [[nodiscard]] const PreprocessedMap& preprocessingMap() const;
   [[nodiscard]] const NodeMap& preprocessedToOriginal() const;
   [[nodiscard]] const NodeMap& originalToPreprocessed() const;

   void doPresolve();

   void initializeBBTree();

   [[nodiscard]] const std::vector<StableSetVariable>& variables() const;
   [[nodiscard]] bool isNewSet(const DenseSet& t_set) const;
   std::size_t addStableSet(DenseSet t_set);
   std::size_t findOrAddStableSet(const DenseSet& t_set);

   void addSolution(std::vector<std::size_t> t_stable_set_indices);

   [[nodiscard]] std::size_t upperBoundUnscaled() const;
   [[nodiscard]] std::size_t lowerBoundUnscaled() const;

   [[nodiscard]] std::size_t upperBound() const;
   [[nodiscard]] std::size_t lowerBound() const;
   [[nodiscard]] double fractionalLowerBound() const;
   [[nodiscard]] SetColoring incumbentUnscaled() const;
   [[nodiscard]] NodeColoring incumbent() const;

   //These functions should only be called either during (or before) processing of the root node,
   //or after a node has been solved and we evaluate the b&b tree
   void updateLowerBound(std::size_t t_lb); //Note the lower bounds are wtihout the offset from the presolving of the root node!
   void updateFractionalLowerBound(double t_fractional_lb);
   void updateTreeBounds();

   BBNode&& popNextNode();
   BBNode &&popNodeWithID(node_id id);

   std::vector<node_id> createChildren(const BBNode& t_node, ColorNodeWorker& t_nodeWorker);

   [[nodiscard]] std::size_t numProcessedNodes() const;
   [[nodiscard]] std::size_t numOpenNodes() const;
   [[nodiscard]] bool hasOpenNodes() const;

   void startSolveTime();
   void reset(DenseGraph t_graph);

   [[nodiscard]] std::chrono::duration<double> timeSinceStart() const;
   [[nodiscard]] bool checkTimelimitHit() const;
   [[nodiscard]] bool checkNodeLimitHit() const;

   [[nodiscard]] const Settings& settings() const;
   void displayHeader(std::ostream& t_stream) const;
   void display(std::ostream& t_stream);

   void addPricingIterations(std::size_t count);
   void addLPIterations(std::size_t count);
 private:
   // Original Problem data
   DenseGraph m_originalGraph;
   // Problem after presolving.
   DenseGraph m_preprocessedGraph;
   PreprocessedMap m_preprocessedToOriginal;

   // Shared solution data, not constant. The variables constraints and colorings are in terms of the preprocessed graph.
   std::vector<StableSetVariable> m_variables;
   // upper bound info
   std::vector<std::vector<std::size_t>> m_colorings;
   std::size_t m_incumbent_index;
   std::size_t m_upperBound;
   // lower bound info; //TODO: track or query?
   std::size_t m_lowerBound;
   double m_fractionalLowerBound;
   //TODO: maybe add lower bound certificates in form of clique / mycielski graphs here

   BBTree m_tree;

   std::chrono::high_resolution_clock::time_point m_start_solve_time;

   //printing options
   int m_printheader_counter;

   Settings& m_settings;
   std::size_t m_lpIterations;
   std::size_t m_pricingIterations;
};
}
#endif // PCOG_SOLUTIONDATA_HPP
