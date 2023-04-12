//
// Created by rolf on 6-12-22.
//

#include <utility>

#include "pcog/ColorSolver.hpp"
namespace pcog {

void ColorSolver::setProblem(std::string t_name, DenseGraph t_graph) {
   m_problemName = std::move(t_name);
   m_solData.reset(std::move(t_graph));
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
   m_solData.startSolveTime();
   if (m_status == SolverStatus::PROBLEM_INITIALIZED) {
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
   m_solData.doPresolve();
   setStatus(SolverStatus::PRESOLVED);
   auto time_end = std::chrono::high_resolution_clock::now();
   m_statistics.m_presolve_time = std::chrono::duration<double>(time_end - time_start).count();
   return m_status;
}
void ColorSolver::branchAndBound() {
   auto time_start = std::chrono::high_resolution_clock::now();
   setStatus(SolverStatus::SOLVING);
   m_solData.initializeBBTree();
   // Branch-and-bound loop
   while (m_solData.hasOpenNodes()) {

      m_worker.processNextNode(m_solData);

      if (m_solData.checkNodeLimitHit() ||
          m_solData.checkTimelimitHit()) {
         break;
      }
      m_solData.display(std::cout);
   }
   auto time_end = std::chrono::high_resolution_clock::now();
   m_statistics.m_branch_and_bound_time =
       std::chrono::duration<double>(time_end - time_start).count();
   // TODO: be careful with termination from time limit in processNode before
   // branching when only one node is left
   m_status = m_solData.lowerBound() == m_solData.upperBound()
                  ? SolverStatus::SOLVED_OPTIMALLY : SolverStatus::SOLVED_SUBOPTIMALLY ;
}
void ColorSolver::recordStatistics() {
   m_statistics.m_total_solve_time = m_solData.timeSinceStart();
}


} // namespace pcog