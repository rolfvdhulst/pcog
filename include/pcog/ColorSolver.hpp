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
   ColorSolver() : m_solData{m_settings}{};

   void setProblem(std::string t_name, DenseGraph t_graph);

   SolverStatus solve();

   Settings& settings() { return m_settings;}
   const Settings& settings() const {return m_settings;}
 private:
   SolverStatus presolve();
   void branchAndBound();
   void setStatus(SolverStatus t_status);
   void recordStatistics();
   void displayEndResult(std::ostream& stream);

   SolverStatus m_status = SolverStatus::NO_PROBLEM;

   //Shared data, constant
   std::string m_problemName;
   //settings (TODO implement gap checks etc.)
   Settings m_settings;

   SolutionData m_solData;

   //statistics
   Statistics m_statistics;
};
} // namespace pcog
#endif // PCOG_SRC_COLORSOLVER_HPP
