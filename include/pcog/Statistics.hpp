// Created by rolf on 4-4-23.
//

#ifndef PCOG_SRC_STATISTICS_HPP
#define PCOG_SRC_STATISTICS_HPP

#include <chrono>
#include <iostream>

struct Statistics {
   void reset();
   std::chrono::duration<double> m_presolve_time;
   std::chrono::duration<double> m_branch_and_bound_time;
   std::chrono::duration<double> m_total_solve_time;

};

#endif // PCOG_SRC_STATISTICS_HPP
