//
// Created by rolf on 25-2-23.
//

#include "pcog/ColorNodeWorker.hpp"
#include "pcog/ColorSolver.hpp"
#include "pcog/Branching.hpp"

namespace pcog {
void ColorNodeWorker::processNode(BBNode &t_node, const ColorSolver &t_solver) {
   // Node preprocessing to compute focus graph and complete focus graph
   setupGraphs(t_node,t_solver);
   // Set up LP representation
   setupLP(t_node,t_solver);
   m_focusNode = t_node.id();
   // Price-lp loop
   // Branch: Check if node is cut off. If not, determine the branching vertex
   // pair of the graph.
}
void ColorNodeWorker::setupGraphs(BBNode &node, const ColorSolver &solver) {

   auto [completeGraph, result] = constructBranchedGraphs(
       solver.preprocessedGraph(), node.branchDecisions(), node.lowerBound());
   //TODO: check if moving is necesary here or if auto-detected.
   m_completeFocusGraph = completeGraph;
   m_focusGraph = result.graph;
   m_mapToPreprocessed = result.map;
}
void ColorNodeWorker::setupLP(BBNode &bb_node, const ColorSolver& t_solver) {
   //TODO: adjust method below so that the lp is not reinitialized from scratch every time
   m_lpSolver.clear();
   m_lpSolver.setObjectiveSense(ObjectiveSense::MINIMIZE);
   for(Node vertex = 0; vertex < m_mapToPreprocessed.newToOldIDs.size(); ++vertex){
      m_lpSolver.addRow({},1.0,std::nullopt);
   }
   for(const auto& variable : t_solver.variables()){
      std::vector<ColElem> columnEntries;
      for(const auto& node : variable.set()){
         auto index = m_mapToPreprocessed.oldToNewIDs[node]; //TODO: check if there is not a 'double mapping' here
         if(index != INVALID_NODE){
            columnEntries.push_back(ColElem(index,1.0)); //TODO: fix conversion to int
         }
      }

      //TODO: experiment with upper bound set to infinity (e.g. the 'virtual' bounds scip uses)
      //Fix variables to zero from ryan-foster branching
      double upperBound = 1.0;
      const auto& set = variable.set();
      for(const auto& decision : bb_node.branchDecisions()){
         //Check if the stable set is viable in the Subproblem
         if(
             (decision.type == BranchType::DIFFER && set.contains(decision.first) && set.contains(decision.second)) ||
             (decision.type == BranchType::SAME && (set.contains(decision.first) != set.contains(decision.second)))
             ){
            upperBound = 0.0;
            break;
         }
      }
      m_lpSolver.addColumn(columnEntries,1.0,0.0,upperBound);
   }
}
} // namespace pcog
