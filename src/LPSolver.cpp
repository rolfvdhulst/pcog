//
// Created by rolf on 18-11-22.
//

#include "pcog/LPSolver.hpp"
using namespace soplex;
namespace pcog {
void LPSolver::setObjectiveSense(ObjectiveSense t_objective) {
   switch (t_objective) {
   case ObjectiveSense::MINIMIZE:
      m_soplex.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
      return;
   case ObjectiveSense::MAXIMIZE:
      m_soplex.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
      return;
   }
}
void LPSolver::addRow(const std::vector<RowElem>& t_rowElems, std::optional<double> t_lhs, std::optional<double> t_rhs) {
   if(!t_lhs.has_value()){
      t_lhs = soplex::infinity;
   }
   if(!t_rhs.has_value()){
      t_rhs = soplex::infinity;
   }
   DSVector row(static_cast<int>(t_rowElems.size()));
   for(const auto& elem : t_rowElems){
      row.add(elem.column,elem.value);
   }

   m_soplex.addRowReal(LPRow(t_lhs.value(),row,t_rhs.value()));
}
void LPSolver::addColumn(const std::vector<ColElem>& t_colElements, double t_objective, double t_lowerBound, double t_upperBound){
   DSVector col(static_cast<int>(t_colElements.size()));
   for(const auto& elem : t_colElements){
      col.add(elem.row,elem.value);
   }
   m_soplex.addColReal(LPCol(t_objective,col,t_upperBound,t_lowerBound));
}
bool LPSolver::solve() {

   auto status = m_soplex.optimize();

   return status == SPxSolver::Status::OPTIMAL;
}
RowVector LPSolver::getDualSolution() {
   assert(m_soplex.status() == SPxSolver::OPTIMAL);
   DVector spxVector(m_soplex.numRows()); //TODO: how do we know how much space to allocate?
   m_soplex.getDual(spxVector);

   RowVector rowVector; //TODO reserve/allocate space
   for (int i = 0; i < spxVector.dim(); ++i) {
      rowVector.emplace_back(i,spxVector[i]);
   }
   return rowVector;
}

ColVector LPSolver::getPrimalSolution() {
   assert(m_soplex.status() == SPxSolver::OPTIMAL);
   DVector spxVector(m_soplex.numCols()); //TODO: how do we know how much space to allocate?
   m_soplex.getPrimal(spxVector);

   ColVector colVector; //TODO reserve/allocate space
   for (int i = 0; i < spxVector.dim(); ++i) {
      colVector.emplace_back(i,spxVector[i]);
   }
   return colVector;
}
double LPSolver::objective() {
   return m_soplex.objValueReal();
}
void LPSolver::clear() {
   m_soplex.clearLPReal();
}

} // namespace pcog
