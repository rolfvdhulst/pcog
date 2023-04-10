// Created by rolf on 4-4-23.
//

#ifndef PCOG_SRC_STATISTICS_HPP
#define PCOG_SRC_STATISTICS_HPP

#include <chrono>
#include <iostream>

struct Statistics {
   [[nodiscard]] double timeSinceStart() const;
   void reset();
   void displayHeader(std::ostream& t_stream) const;
   void display(std::ostream& t_stream);
   double m_presolve_time;
   double m_branch_and_bound_time;
   double m_total_solve_time;

};

#endif // PCOG_SRC_STATISTICS_HPP
