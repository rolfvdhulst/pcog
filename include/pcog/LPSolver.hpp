//
// Created by rolf on 18-11-22.
//

#ifndef PCOG_SRC_LPSOLVER_HPP
#define PCOG_SRC_LPSOLVER_HPP

#include "soplex.h"
#include "soplex/spxbasis.h"
#include <optional>

namespace pcog {

using RowIdx = int;
using ColIdx = int;

struct RowElem{
   RowElem() = default;
   RowElem(ColIdx t_col, double t_val) : column{t_col},value{t_val}{};
   ColIdx column;
   double value;
};
struct ColElem{
   ColElem() = default;
   ColElem(RowIdx t_row, double t_val) : row{t_row},value{t_val}{};
   RowIdx row;
   double value;
};
using RowVector = std::vector<RowElem>; //TODO: definitions are confusing when working with LP solutions
using ColVector = std::vector<ColElem>;

enum class ObjectiveSense { MINIMIZE, MAXIMIZE };

//TODO: how to store basis efficiently, and in particular how to change it when columns or rows are added in order to hot start?

enum class LPSolverStatus {INFEASIBLE, OPTIMAL, ERROR};
struct LPBasis{
   std::vector<soplex::SPxSolver::VarStatus> rowStatus;
   std::vector<soplex::SPxSolver::VarStatus> colStatus;
   friend class LPSolver;
};
/// \brief class which holds the LP solver
class LPSolver {
 public:
   /// Sets the LP solver to minimize/maximize.
   /// \param t_objective the desired objective sense which the LP solver is set
   /// to.
   void setObjectiveSense(ObjectiveSense t_objective);

   /// Adds a linear row with the form lhs <= a^t x <= rhs to the linear program.
   /// \param t_rowElements The a^t vector.
   /// \param t_lhs Left hand side of the equation. Use nullopt to indicate infinity
   /// \param t_rhs Right hand side of the equation. Use nullopt to indicate infinity
   void addRow(const RowVector& t_rowElements, std::optional<double> t_lhs, std::optional<double> t_rhs);

   /// Adds a column to the linear program with given objective and bounds.
   /// \param t_colElements The vector of column entries for this variable
   /// \param t_objective The objective coefficient for this variable
   /// \param t_lowerBound The lower bound on the variable. Use -soplex::infinity to indicate an unbounded variable
   /// \param t_upperBound The upper bound on the variable. Use soplex::infinity to indicate an unbounded variable
   void addColumn(const ColVector& t_colElements, double t_objective, double t_lowerBound, double t_upperBound);

   void addColumns(const std::vector<ColVector>& t_columnElements, std::vector<double> objective,
                   std::vector<double> t_lowerBound, std::vector<double> t_upperBound);
   /// Solves the linear program
   /// \return true if solved successfully, false if some error occurred.
   bool solve();

   //Following functions are unfortunately not const because of SoPlex
   /// Gets the dual solution vector of the solved LP
   /// \return
   RowVector getDualSolution();

   ///
   /// \return The primal solution of the solved LP
   RowVector getPrimalSolution();

   ///
   /// \return Vector containing upper bounds on the variables
   RowVector columnUpperBounds();

   /// Returns the current LP objective
   /// \return Returns the current LP objective
   double objective();

   void clear();
   //TODO: functions for setting/getting the basis

   LPSolverStatus status();
   std::size_t numRows();
   std::size_t numCols();

   /// \return returns the current LP basis
   LPBasis getLPBasis();

   void setBasis(const LPBasis& basis);
 private:
   soplex::SoPlex m_soplex;
};
} // namespace pcog

#endif // PCOG_SRC_LPSOLVER_HPP
