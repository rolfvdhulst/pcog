//
// Created by rolf on 6-12-22.
//

#ifndef PCOG_SRC_COLORSOLVER_HPP
#define PCOG_SRC_COLORSOLVER_HPP

#include "ColorNodeWorker.hpp"
#include "SolutionData.hpp"
#include "Settings.hpp"
#include "Statistics.hpp"
#include <chrono>

namespace pcog {


enum class SolverStatus{
   NO_PROBLEM,
   PROBLEM_INITIALIZED, /// The solver has read a problem but has not started solving yet
   PRESOLVING, /// The solver is currently preprocessing the problem
   PRESOLVED, /// The solver has read the problem and preprocessed it
   SOLVING, /// The solver is still currently solving
   SOLVED_SUBOPTIMALLY, /// The solution process was terminated early due to some valid condition, such as a time-limit or node-limit
   SOLVED_OPTIMALLY, /// An optimal solution was found during presolving or solving
   ERROR /// An unexpected error occurred during the solution process
};
/// Class which solves the graph coloring problem exactly.
/// This class' purpose is to contain global information such as the graph,
/// a solution pool, the search tree and it provides a simple user interface
class ColorSolver {
 public:
   ColorSolver() = default;

   void setProblem(std::string t_name, DenseGraph t_graph);

   SolverStatus solve();

   [[nodiscard]] const DenseGraph& preprocessedGraph() const {return m_preprocessedGraph;}
   [[nodiscard]] const std::vector<StableSetVariable>& variables() const {return m_variables;}

   Settings& settings() { return m_settings;}
   const Settings& settings() const {return m_settings;}
 private:
   SolverStatus presolve();
   void branchAndBound();
   void setStatus(SolverStatus t_status);
   void recordStatistics();
   bool checkTimelimitHit() const{
      if(m_settings.timeLimit() == NO_TIME_LIMIT){
         return false;
      }
      auto time = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration(m_statistics.m_start_solve_time-time);
      return duration.count() >= m_settings.timeLimit();
   }
   SolverStatus m_status = SolverStatus::NO_PROBLEM;
   ColorNodeWorker m_worker;

   //Shared data, constant
   std::string m_problemName;


   //settings (TODO implement checks)
   Settings m_settings;

   //statistics
   Statistics m_statistics;
};
} // namespace pcog
#endif // PCOG_SRC_COLORSOLVER_HPP
