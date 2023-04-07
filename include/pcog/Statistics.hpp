// Created by rolf on 4-4-23.
//

#ifndef PCOG_SRC_STATISTICS_HPP
#define PCOG_SRC_STATISTICS_HPP

#include <chrono>
struct Statistics {
   [[nodiscard]] double timeSinceStart() const;
   void reset();
   void displayHeader(std::ostream& t_stream) const;
   void display(std::ostream& t_stream);
   std::chrono::high_resolution_clock::time_point m_start_solve_time;
   double m_presolve_time;
   double m_branch_and_bound_time;
   double m_total_solve_time;
   std::size_t m_num_processed_nodes;
   std::size_t m_num_open_nodes;

   std::size_t m_lb;
   std::size_t m_ub;
   double m_fractional_lb;
   //printing options
   int m_printheader_counter;
};

#endif // PCOG_SRC_STATISTICS_HPP
