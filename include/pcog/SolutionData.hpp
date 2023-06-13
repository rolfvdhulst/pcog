//
// Created by rolf on 7-4-23.
//

#ifndef PCOG_SOLUTIONDATA_HPP
#define PCOG_SOLUTIONDATA_HPP

#include <utility>

#include "BBNode.hpp"
#include "pcog/utilities/DenseGraph.hpp"
#include "pcog/utilities/DenseSet.hpp"
#include "pcog/ColorNodeWorker.hpp"
#include "LocalSolutionData.hpp"
#include "Preprocessing.hpp"
#include "Settings.hpp"

namespace pcog {

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
   void runBranchAndBound();


   [[nodiscard]] const std::vector<StableSetVariable>& variables() const;
   [[nodiscard]] bool isNewSet(const DenseSet& t_set) const;
   std::size_t addStableSet(DenseSet t_set);
   std::size_t findOrAddStableSet(const DenseSet& t_set, std::size_t checkNewFromIndex = 0);

   void addSolution(std::vector<std::size_t> t_stable_set_indices);

   [[nodiscard]] std::size_t upperBoundUnscaled() ;
   [[nodiscard]] std::size_t lowerBoundUnscaled() ;

   [[nodiscard]] std::size_t upperBound() ;
   [[nodiscard]] std::size_t lowerBound() ;
   [[nodiscard]] double fractionalLowerBound();
   [[nodiscard]] SetColoring incumbentUnscaled() const;
   [[nodiscard]] NodeColoring incumbent() const;

   std::optional<BBNode> popNextNode(ColorNodeWorker& t_nodeWorker);
   std::optional<BBNode> branchAndPopNode(BBNode& node, ColorNodeWorker& t_nodeWorker);
   [[nodiscard]] std::size_t numProcessedNodes();
   [[nodiscard]] std::size_t numOpenNodes();

   void startSolveTime();
   void reset(DenseGraph t_graph);

   [[nodiscard]] std::chrono::duration<double> timeSinceStart() const;
   [[nodiscard]] bool checkTimelimitHit() const;
   [[nodiscard]] bool checkNodeLimitHit() const;

   [[nodiscard]] const Settings& settings() const;

   void displayHeader(std::ostream& t_stream) const;
   void display(std::ostream& t_stream);

   void synchronizeLocalDataStatistics(LocalSolutionData& t_localData);
   void writeLocalVarsToGlobal(LocalSolutionData& t_localSolutionData,std::atomic_bool& stop);
   void writeLocalSolutionsToGlobal(LocalSolutionData& t_localSolutionData, std::atomic_bool& stop);
   void syncLocalVarsWithGlobal(LocalSolutionData& t_localSolutionData);
   void syncLocalLowerBound(LowerBoundInfo lbInfo, std::size_t worker_id);

   bool doRecomputeLowerBound(std::atomic_bool& stop);
 private:
   std::vector<node_id> createChildren(const BBNode& t_node, ColorNodeWorker& t_nodeWorker);
   void pruneUpperBound(std::size_t upperBound,std::atomic_bool& stop);
   void initializeBBTree();
   std::optional<BBNode> pickNextNode(BBNode& t_node, ColorNodeWorker &t_nodeWorker, const std::vector<node_id>& children);
   void writeLocalVariablesToGlobal(LocalSolutionData& t_localSolutionData);

   //Recomputes the lower bound based on the tree and worker info.
   //It's your responsibility to lock the mutex before calling these two functions (hence private!)
   //Returns true if the lower bound was improved
   bool recomputeLowerBound();
   bool assignLB(double t_frac_lb, std::size_t t_lb);

   // Original Problem data
   DenseGraph m_originalGraph;
   // Problem after presolving.
   DenseGraph m_preprocessedGraph;
   PreprocessedMap m_preprocessedToOriginal;

   // Shared solution data, not constant. The variables, constraints and colorings are in terms of the preprocessed graph.
   std::mutex m_variable_mutex;
   std::vector<StableSetVariable> m_variables;
   // upper bound info
   std::mutex m_upperBound_mutex;
   std::vector<std::vector<std::size_t>> m_colorings;
   std::size_t m_incumbent_index;
   std::atomic_size_t m_upperBound; //This is atomic because the worker threads access them from time to time but it is unlikely to be updated much.
                                    //Putting it into the Mutex would likely be slower

   std::vector<ColorNodeWorker> m_workers; //Worker threads which execute the branch& price loop for a single node

   // lower bound info
   std::mutex m_lowerBound_mutex; //This mutex is locked whenever the B&B tree or the lower bound changes
   std::vector<std::optional<LowerBoundInfo>> m_processing_node_lower_bounds; //The lower bounds of the nodes which are currently being processed
   std::size_t m_lowerBound;
   double m_fractionalLowerBound;
   BBTree m_tree;
   //TODO: maybe add lower bound certificates in form of clique / mycielski graphs here

   std::chrono::high_resolution_clock::time_point m_start_solve_time;

   //printing options
   std::mutex m_printing_mutex;
   int m_printheader_counter;

   Settings& m_settings;
   std::atomic_size_t m_lpIterations;
   std::atomic_size_t m_pricingIterations;
};
}
#endif // PCOG_SOLUTIONDATA_HPP
