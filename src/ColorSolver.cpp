//
// Created by rolf on 6-12-22.
//

#include <utility>

#include "pcog/ColorSolver.hpp"
#include <chrono>
namespace pcog{

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
}
SolverStatus ColorSolver::solve() {
   //TODO: check if problem is initialized and handle cases where we have already solved the problem?

   auto time_start = std::chrono::high_resolution_clock::now();
   //First presolve the problem
   {
      setStatus(SolverStatus::PRESOLVING);
      {
         auto result = preprocessOriginalGraph(m_originalGraph);
         m_preprocessedGraph = result.graph;
         m_preprocessedToOriginal = result.map;
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

         m_worker.processNode(bb_node,*this); //TODO: ugly *this usage, use intermediary struct to save data.
         std::cout<<"Processed node, lb: "<<bb_node.fractionalLowerBound()<<std::endl;
         if (bb_node.status() == BBNodeStatus::BRANCHED){

            m_tree.createChildren(bb_node.id());
         }
      }
      //TODO: set status correctly here
   }
   auto time_end = std::chrono::high_resolution_clock::now();
   std::cout<<(time_end-time_start).count() /1e9 <<" seconds to optimize!\n";

   return m_status;
}
void ColorSolver::setStatus(SolverStatus t_status) {
   m_status = t_status; //TODO: add lock?
}
}