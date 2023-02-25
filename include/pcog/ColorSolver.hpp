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
 private:
   DenseSet m_set;
};

class NodeConstraint{
 public:
   explicit NodeConstraint(Node t_node) : m_node{t_node}{};
 private:
   Node m_node;
};

enum class SolverStatus{
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
 private:
   void setStatus(SolverStatus t_status);
   SolverStatus m_status;

   std::string m_problemName;
   //Original Problem data
   DenseGraph m_originalGraph;
   //Problem after presolving
   DenseGraph m_preprocessedGraph;
   PreprocessedMap m_preprocessedToOriginal;

   ColorNodeWorker m_worker;
   //Shared solution data. The variables constraints and colorings are in terms of the preprocessed graph.
   std::vector<StableSetVariable> m_variables;
   std::vector<NodeConstraint> m_constraints;
   std::vector<std::vector<std::size_t>> m_colorings;

   BBTree m_tree;

   //TODO: statistics (LP its, nodes,
};
} // namespace pcog
#endif // PCOG_SRC_COLORSOLVER_HPP
