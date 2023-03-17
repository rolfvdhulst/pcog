//
// Created by rolf on 18-11-22.
//

#include "pcog/LPSolver.hpp"
using namespace soplex;
namespace pcog {

LPSolverStatus convertStatus(SPxSolver::Status status){
   switch (status) {
   case SPxSolverBase<double>::Status::OPTIMAL:
      return LPSolverStatus::OPTIMAL;
   case SPxSolverBase<double>::Status::INFEASIBLE:
      return LPSolverStatus::INFEASIBLE;
//   case SPxSolverBase<double>::Status::UNBOUNDED:
//   case SPxSolverBase<double>::Status::INForUNBD:
//   case SPxSolverBase<double>::Status::OPTIMAL_UNSCALED_VIOLATIONS:
//   case SPxSolverBase<double>::Status::ERROR:
//   case SPxSolverBase<double>::Status::NO_RATIOTESTER:
//   case SPxSolverBase<double>::Status::NO_PRICER:
//   case SPxSolverBase<double>::Status::NO_SOLVER:
//   case SPxSolverBase<double>::Status::NOT_INIT:
//   case SPxSolverBase<double>::Status::ABORT_EXDECOMP:
//   case SPxSolverBase<double>::Status::ABORT_DECOMP:
//   case SPxSolverBase<double>::Status::ABORT_CYCLING:
//   case SPxSolverBase<double>::Status::ABORT_TIME:
//   case SPxSolverBase<double>::Status::ABORT_ITER:
//   case SPxSolverBase<double>::Status::ABORT_VALUE:
//   case SPxSolverBase<double>::Status::SINGULAR:
//   case SPxSolverBase<double>::Status::NO_PROBLEM:
//   case SPxSolverBase<double>::Status::REGULAR:
//   case SPxSolverBase<double>::Status::RUNNING:
//   case SPxSolverBase<double>::Status::UNKNOWN:
   default:
      return LPSolverStatus::ERROR; //TODO: create more cases
   }
}
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
   m_soplex.setIntParam(soplex::SoPlexBase<double>::VERBOSITY,0);
   auto status = m_soplex.optimize();

   return status == SPxSolver::Status::OPTIMAL;
}
RowVector LPSolver::getDualSolution() {
   assert(m_soplex.status() == SPxSolver::OPTIMAL);
   DVector spxVector(m_soplex.numRows()); //TODO: how do we know how much space to allocate?
   bool result = m_soplex.getDual(spxVector); //TODO: error handling
   assert(result);
   RowVector rowVector; //TODO reserve/allocate space
   for (int i = 0; i < spxVector.dim(); ++i) {
      rowVector.emplace_back(i,spxVector[i]);
   }
   return rowVector;
}
RowVector LPSolver::columnUpperBounds() {
   DVector spxVector(m_soplex.numCols());
   m_soplex.getUpperReal(spxVector); //TODO: error handling
   RowVector vector; //TODO reserve/allocate space
   for (int i = 0; i < spxVector.dim(); ++i) {
      vector.emplace_back(i,spxVector[i]);
   }
   return vector;
}
RowVector LPSolver::getPrimalSolution() {
   assert(m_soplex.status() == SPxSolver::OPTIMAL);
   DVector spxVector(m_soplex.numCols());
   m_soplex.getPrimal(spxVector);

   RowVector colVector; //TODO reserve/allocate space
   for (int i = 0; i < spxVector.dim(); ++i) {
      if(spxVector[i] != 0.0){ // 'soft' checking for 0.0 because we let users of this function decide what to do with very small numerical values
         colVector.emplace_back(i,spxVector[i]);
      }
   }
   return colVector;
}
double LPSolver::objective() {
   return m_soplex.objValueReal();
}
void LPSolver::clear() {
   m_soplex.clearLPReal();
}
LPSolverStatus LPSolver::status(){
return convertStatus(m_soplex.status());
}
std::size_t LPSolver::numRows() {
return m_soplex.numRows();
}
std::size_t LPSolver::numCols() {
   return m_soplex.numCols();
}
LPBasis LPSolver::getLPBasis() {
   LPBasis basis;
   basis.rowStatus.resize(m_soplex.numRows());
   basis.colStatus.resize(m_soplex.numCols());
   m_soplex.getBasis(basis.rowStatus.data(),basis.colStatus.data());

   return basis;
}
void LPSolver::setBasis(const LPBasis &basis) {

   m_soplex.setBasis(basis.rowStatus.data(),basis.colStatus.data());
}
void LPSolver::addColumns(const std::vector<ColVector> &t_columnElements,
                          std::vector<double> objective,
                          std::vector<double> t_lowerBound,
                          std::vector<double> t_upperBound) {
   if(t_columnElements.empty()){
      return;
   }
   LPColSetReal colSet;
   for(std::size_t i = 0; i < t_columnElements.size(); ++i){
      DSVector col(static_cast<int>(t_columnElements[i].size()));
      for(const auto& elem : t_columnElements[i]){
         col.add(elem.row,elem.value);
      }
      colSet.add(LPCol(objective[i],col,t_upperBound[i],t_lowerBound[i]));
   }
   m_soplex.addColsReal(colSet);
}

} // namespace pcog
