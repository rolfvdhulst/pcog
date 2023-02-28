//
// Created by rolf on 6-12-22.
//

#ifndef PCOG_SRC_COLORSOLVER_HPP
#define PCOG_SRC_COLORSOLVER_HPP

#include "BBNode.hpp"
#include "ColorNodeWorker.hpp"
#include "DenseGraph.hpp"
#include "Preprocessing.hpp"

namespace pcog {

class StableSetVariable{
 public:
   StableSetVariable(DenseSet set) : m_set{set}{};
   [[nodiscard]] const DenseSet& set() const { return m_set;}
 private:
   DenseSet m_set;
};

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
   bool isNewSet(const DenseSet& set) const{
      for(const auto& var : m_variables){
         if(var.set() == set){
            return false;
         }
      }
      return true;
   }
   void addStableSet(const DenseSet& set){
      assert(m_preprocessedGraph.setIsStable(set));
      m_variables.emplace_back(set);
   }
 private:
   void setStatus(SolverStatus t_status);
   SolverStatus m_status = SolverStatus::NO_PROBLEM;
   ColorNodeWorker m_worker;

   //Shared data, constant
   std::string m_problemName;
   //Original Problem data
   DenseGraph m_originalGraph;
   //Problem after presolving.
   DenseGraph m_preprocessedGraph;
   PreprocessedMap m_preprocessedToOriginal;

   //Shared solution data, not constant. The variables constraints and colorings are in terms of the preprocessed graph.
   std::vector<StableSetVariable> m_variables;
   std::vector<std::vector<std::size_t>> m_colorings;
   BBTree m_tree;

   //TODO: statistics (LP its, nodes,
};
} // namespace pcog
#endif // PCOG_SRC_COLORSOLVER_HPP
