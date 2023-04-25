//
// Created by rolf on 4-4-23.
//

#include "pcog/Statistics.hpp"


void Statistics::reset() {
   m_presolve_time =  std::chrono::duration<double>(0.0);
   m_branch_and_bound_time =  std::chrono::duration<double>(0.0);
   m_total_solve_time =  std::chrono::duration<double>(0.0);

}

