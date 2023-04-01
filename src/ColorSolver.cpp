//
// Created by rolf on 6-12-22.
//

#include <utility>

#include "pcog/ColorSolver.hpp"
namespace pcog {

void ColorSolver::setProblem(std::string t_name, DenseGraph t_graph) {
   m_problemName = std::move(t_name);
   m_originalGraph = std::move(t_graph);
   m_preprocessedGraph.clear();
   m_preprocessedToOriginal.clear();
   m_variables.clear();
   m_colorings.clear();
   m_tree.clear();

   m_upperBound = std::numeric_limits<std::size_t>::max();
   m_incumbent_index = std::numeric_limits<std::size_t>::max();
   m_status = SolverStatus::PROBLEM_INITIALIZED;
}
SolverStatus ColorSolver::solve() {
   // Behavior at start of solve:
   // If NO_PROBLEM => return
   // If PROBLEM_INITIALIZED => presolve the problem and solve normally
   // If SOLVING or PRESOLVING => returns ERROR but does not change solve
   // status, no concurrent solves currently allowed If PRESOLVED => skip
   // presolving and solve normally If SOLVED_OPTIMALLY => return

   // TODO: check path of the following and add some asserts
   // ERROR?
   // SOLVED_SUBOPTIMALLY => continue branch-and-bound

   switch (m_status) {
   case SolverStatus::NO_PROBLEM:
   case SolverStatus::SOLVED_OPTIMALLY:
      return m_status;
   case SolverStatus::SOLVING:
   case SolverStatus::PRESOLVING:
      return SolverStatus::ERROR; // We do not set status to error so that
                                  // concurrent solve could continue
   default: // TODO: ERROR and SOLVED_SUBOPTIMALLY; what to do with them
      break;
   }
   m_start_solve_time = std::chrono::high_resolution_clock::now();
   if (m_status == SolverStatus::PROBLEM_INITIALIZED) {
      cleanStatistics();
      // First presolve the problem
      presolve();
      if (m_status == SolverStatus::ERROR ||
          m_status == SolverStatus::SOLVED_OPTIMALLY) {
         recordStatistics();
         return m_status;
      }
   }
   assert(m_status == SolverStatus::PRESOLVED);
   branchAndBound();
   recordStatistics();
   return m_status;
}
void ColorSolver::setStatus(SolverStatus t_status) {
   m_status = t_status; // TODO: add lock?
}
SolverStatus ColorSolver::presolve() {
   auto time_start = std::chrono::high_resolution_clock::now();
   assert(m_status == SolverStatus::PROBLEM_INITIALIZED);
   setStatus(SolverStatus::PRESOLVING);
   // TODO: find an initial coloring using a simple greedy method.

   // This can help to judge if it's 'worth' to find a clique for the lower
   // bound; If the minimal degree of the graph is larger than the chromatic
   // number (which typically happens for dense graphs), then the low-degree
   // preprocessing rule will not find anything unless some other rule finds a
   // reduction which lowers a node degree of the graph below the chromatic
   // number
   auto result = preprocessOriginalGraph(m_originalGraph);
   m_preprocessedGraph = result.graph;
   m_preprocessedToOriginal = result.map;
   // Also add coloring variables as initial solution

   setStatus(SolverStatus::PRESOLVED);
   auto time_end = std::chrono::high_resolution_clock::now();
   m_presolve_time =
       std::chrono::duration<double>(time_end - time_start).count();
   return m_status;
}
void ColorSolver::branchAndBound() {
   auto time_start = std::chrono::high_resolution_clock::now();
   setStatus(SolverStatus::SOLVING);
   // initialize tree
   m_tree.createRootNode(m_preprocessedGraph.numNodes());
   // Branch-and-bound loop
   while (m_tree.hasOpenNodes()) {
      BBNode &bb_node = m_tree.popNextNode();

      m_worker.processNode(bb_node,
                           *this); // TODO: ugly *this usage, use
                                   // intermediary struct to save data.
      // TODO: check for time limit in processNode
      std::cout << "Processed node, lb: " << bb_node.fractionalLowerBound()
                << " open nodes: " << m_tree.numOpenNodes()
                << " total nodes: " << m_tree.numTotalNodes()
                << " vars: " << m_variables.size() << std::endl;

      if (bb_node.status() == BBNodeStatus::BRANCHED) {
         m_tree.createChildren(bb_node.id(), m_worker);
      }
      ++m_num_processed_nodes;

      if (m_num_processed_nodes == m_settings.nodeLimit() ||
          checkTimelimitHit()) {
         break;
      }
   }
   auto time_end = std::chrono::high_resolution_clock::now();
   m_branch_and_bound_time =
       std::chrono::duration<double>(time_end - time_start).count();
   // TODO: be careful with termination from time limit in processNode before
   // branching when only one node is left
   m_status = m_tree.hasOpenNodes() ? SolverStatus::SOLVED_SUBOPTIMALLY
                                    : SolverStatus::SOLVED_OPTIMALLY;
}
void ColorSolver::recordStatistics() {
   auto end_time = std::chrono::high_resolution_clock::now();
   m_total_solve_time =
       std::chrono::duration<double>(end_time - m_start_solve_time).count();
}
void ColorSolver::cleanStatistics() {
   m_presolve_time = 0.0;
   m_branch_and_bound_time = 0.0;
   m_total_solve_time = 0.0;
   m_num_processed_nodes = 0;
}
void ColorSolver::printStatistics(std::ostream &t_out) const {
   t_out << "Total time: " << m_total_solve_time << "\n";
   t_out << "Presolve time: " << m_presolve_time << "\n";
   t_out << "B&B time: " << m_branch_and_bound_time << "\n";
   t_out << "B&B nodes: " << m_num_processed_nodes << "\n";
}
std::size_t ColorSolver::findOrAddStableSet(const DenseSet &set) {
   for(std::size_t i = 0; i < m_variables.size(); ++i){
      if(m_variables[i].set() == set){
         return i;
      }
   }
   std::size_t index = m_variables.size();
   addStableSet(set);
   return index;
}
void ColorSolver::addSolution(const std::vector<std::size_t> &t_color_indices) {
#ifndef NDEBUG
   {
      // assert that coloring is indeed a valid coloring
      DenseSet coveredNodes(m_preprocessedGraph.numNodes());
      for (const auto &index : t_color_indices) {
         assert(index < m_variables.size());
         coveredNodes.inplaceUnion(m_variables[index].set());
      }
      assert(coveredNodes.full());
   }
#endif

   std::size_t ub = t_color_indices.size();
   if (ub < m_upperBound) {
      // prune redundant nodes from the tree
      m_tree.pruneUpperBound(ub);
      // TODO: Additionally, update the LP cutoff limit for the current open workers
      m_upperBound = ub;
      m_incumbent_index = m_colorings.size();
      std::cout << "New incumbent found with " << ub << "-coloring!\n";
   }
   m_colorings.push_back(t_color_indices);
}
} // namespace pcog