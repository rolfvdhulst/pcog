//
// Created by rolf on 18-11-22.
//

#include <gtest/gtest.h>
#include <pcog/LPSolver.hpp>

TEST(LPSolver, basicProblem) {
   //Solve a simple problem to test basic LPSolver functionalities
   pcog::LPSolver solver;
   solver.setObjectiveSense(pcog::ObjectiveSense::MINIMIZE);

   solver.addColumn({}, 3.0, 1.0, soplex::infinity);
   solver.addColumn({}, 2.0, 1.0, soplex::infinity);

   //   m_soplex.setIntParam(SoPlex::VERBOSITY,SoPlex::VERBOSITY_ERROR);

   solver.addRow({{0, 0.2}, {1, 1.0}}, 2.0, soplex::infinity);
   /* solve LP */
   bool good = solver.solve();
   EXPECT_TRUE(good);

   EXPECT_EQ(solver.objective(),6.6);
   pcog::RowVector dualSolution = solver.getDualSolution();
   pcog::ColVector primalSolution = solver.getPrimalSolution();
   ASSERT_EQ(dualSolution.size(),1);
   ASSERT_EQ(primalSolution.size(),2);
   EXPECT_DOUBLE_EQ(primalSolution[0].value,1.0);
   EXPECT_DOUBLE_EQ(primalSolution[1].value,1.8);
   EXPECT_DOUBLE_EQ(dualSolution[0].value,2.0);
}