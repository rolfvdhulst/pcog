//
// Created by rolf on 6-12-22.
//

#include <utility>

#include "pcog/ColorSolver.hpp"

namespace pcog{

void ColorSolver::setProblem(std::string t_name, DenseGraph t_graph) {
   m_problemName = std::move(t_name);
   m_originalGraph = std::move(t_graph);
   m_preprocessedGraph.clear();
   m_preprocessedToOriginal.clear();
   m_variables.clear();
   m_constraints.clear();
   m_colorings.clear();
   m_tree.clear();

}
SolverStatus ColorSolver::solve() {
   //TODO: check if problem is initialized and handle cases where we have already solved the problem?

   //First presolve the problem
   {
      setStatus(SolverStatus::PRESOLVING);
      DenseSet clique; //TODO: find clique
      {
         auto result = preprocessOriginalGraph(m_originalGraph,clique);
         m_preprocessedGraph = result.graph;
         m_preprocessedToOriginal = result.map;
      }
      //set up the constraints, one for each node.
      for(Node node = 0; node < m_preprocessedGraph.numNodes(); ++node){
         m_constraints.emplace_back(node);
      }
      setStatus(SolverStatus::PRESOLVED);
   }

   {
      setStatus(SolverStatus::SOLVING);
      //initialize tree
      m_tree.createRootNode();
      //Branch-and-bound loop
      while(m_tree.hasOpenNodes()){
         BBNode& bb_node = m_tree.popNextNode();
         m_worker.processNode(bb_node);
         if (bb_node.status() == BBNodeStatus::BRANCHED){

            m_tree.createChildren(bb_node.id());
         }
      }
      //TODO: set status correctly here
   }

   return m_status;
}
void ColorSolver::setStatus(SolverStatus t_status) {
   m_status = t_status; //TODO: add lock?
}
}