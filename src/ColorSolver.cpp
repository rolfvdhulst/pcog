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
   m_statistics.m_start_solve_time = std::chrono::high_resolution_clock::now();
   if (m_status == SolverStatus::PROBLEM_INITIALIZED) {
      m_statistics.reset();
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
   m_statistics.m_presolve_time =
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

      if (bb_node.status() == BBNodeStatus::BRANCHED) {
         m_tree.createChildren(bb_node.id(), m_worker);
      }
      ++m_statistics.m_num_processed_nodes;

      if (m_statistics.m_num_processed_nodes == m_settings.nodeLimit() ||
          checkTimelimitHit()) {
         break;
      }
      m_statistics.display(std::cout);
   }
   auto time_end = std::chrono::high_resolution_clock::now();
   m_statistics.m_branch_and_bound_time =
       std::chrono::duration<double>(time_end - time_start).count();
   // TODO: be careful with termination from time limit in processNode before
   // branching when only one node is left
   m_status = m_tree.hasOpenNodes() ? SolverStatus::SOLVED_SUBOPTIMALLY
                                    : SolverStatus::SOLVED_OPTIMALLY;
}
void ColorSolver::recordStatistics() {
   m_statistics.m_total_solve_time = m_statistics.timeSinceStart();
}

void ColorSolver::printStatistics(std::ostream &t_out) const {
   t_out << "Total time: " << m_statistics.m_total_solve_time << "\n";
   t_out << "Presolve time: " << m_statistics.m_presolve_time << "\n";
   t_out << "B&B time: " << m_statistics.m_branch_and_bound_time << "\n";
   t_out << "B&B nodes: " << m_statistics.m_num_processed_nodes << "\n";
}
std::size_t ColorSolver::findOrAddStableSet(const DenseSet &set) {

}
void ColorSolver::addSolution(const std::vector<std::size_t> &t_color_indices) {

}
} // namespace pcog