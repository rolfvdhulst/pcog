//
// Created by rolf on 4-4-23.
//

#include "pcog/Statistics.hpp"
#include <iostream>
#include <iomanip>


void Statistics::reset() {
   m_presolve_time = 0.0;
   m_branch_and_bound_time = 0.0;
   m_total_solve_time = 0.0;
   m_num_processed_nodes = 0;

   m_printheader_counter = 0;
}
void Statistics::displayHeader(std::ostream &t_stream) const {
   if(m_printheader_counter == 0){
      t_stream << "      System      |        Nodes        |      Bounds     |           Work              \n";
   }
   t_stream << "   time |  memory | explored |     open |  lower |  upper |   vars |  LP it | Pric it | \n";
}

void Statistics::display(std::ostream &t_stream) {
   if(m_printheader_counter % 20 == 0){
      displayHeader(t_stream);
   }
   // TODO: include fractional LB and # solutions, gap% and maximum tree depth?
   //<< std::setw(8) << std::setprecision(4) << timeSinceStart() <<"|"
   auto duration = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() -m_start_solve_time);
   auto mem = 0.1;
   //TODO: fix almost all fields being written to (OR pass in solver struct data somehow)

   t_stream << std::fixed << std::setw(7) << std::setprecision(1) << std::right << duration << " | "
            << std::fixed << std::setw(7) << std::setprecision(2) << std::right<< mem <<" | "
            << std::scientific << std::setw(8) << std::setprecision(2) << std::right << m_num_processed_nodes << " | "
            << std::scientific << std::setw(8) << std::setprecision(2) << std::right << m_num_open_nodes << " | "
            << std::scientific << std::setw(6) << std::setprecision(2) << std::right << m_lb << " | "
            << std::scientific << std::setw(6) << std::setprecision(2) << std::right << m_ub << " | "

            << std::endl;

   ++m_printheader_counter;

}
double Statistics::timeSinceStart() const {

   return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - m_start_solve_time).count();
}

