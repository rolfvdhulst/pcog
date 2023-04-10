//
// Created by rolf on 4-4-23.
//

#include "pcog/Statistics.hpp"
#include <iostream>
#include <iomanip>
#include <array>

#ifdef __linux__
#include <sys/sysinfo.h>
#include <unistd.h>
#endif

#ifdef __APPLE__
#include <mach/mach_init.h>
#include <mach/task.h>
#endif

#ifdef _WINDOWS
#include <windows.h>
#else
#include <sys/resource.h>
#endif

// Stackoverflow solution to get memory used by the process. Not perfect.
/// The amount of memory currently being used by this process, in bytes.
/// By default, returns the full virtual arena, but if resident=true,
/// it will report just the resident set in RAM (if supported on that OS).
size_t memory_used (bool resident=false)
{
#if defined(__linux__)
   // Ugh, getrusage doesn't work well on Linux.  Try grabbing info
   // directly from the /proc pseudo-filesystem.  Reading from
   // /proc/self/statm gives info on your own process, as one line of
   // numbers that are: virtual mem program size, resident set size,
   // shared pages, text/code, data/stack, library, dirty pages.  The
   // mem sizes should all be multiplied by the page size.
   size_t size = 0;
   FILE *file = fopen("/proc/self/statm", "r");
   if (file) {
      unsigned long vm = 0;
      fscanf (file, "%ul", &vm);  // Just need the first num: vm size
      fclose (file);
      size = (size_t)vm * getpagesize();
   }
   return size;

#elif defined(__APPLE__)
   // Inspired by:
   // http://miknight.blogspot.com/2005/11/resident-set-size-in-mac-os-x.html
   struct task_basic_info t_info;
   mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
   task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
   size_t size = (resident ? t_info.resident_size : t_info.virtual_size);
   return size;

#elif defined(_WINDOWS)
   // According to MSDN...
   PROCESS_MEMORY_COUNTERS counters;
   if (GetProcessMemoryInfo (GetCurrentProcess(), &counters, sizeof (counters)))
      return counters.PagefileUsage;
   else return 0;

#else
   // No idea what platform this is
   return 0;   // Punt
#endif
}

std::string bytesToString(std::size_t bytes){
   static constexpr std::array<const char*,5> extension = {"B","kB","MB","GB","TB"};
   std::size_t extension_ind = 0;
   double show = bytes;
   while(extension_ind < 5){
      if(show <  1000.0){
         break;
      }
      show /= 1000.0;
      ++extension_ind;
   }
   std::stringstream stream;
   stream << std::setprecision(3) << show;
   std::string str = stream.str() + extension[extension_ind];
   return str;
}

void Statistics::reset() {
   m_presolve_time = 0.0;
   m_branch_and_bound_time = 0.0;
   m_total_solve_time = 0.0;

}
void Statistics::displayHeader(std::ostream &t_stream) const {
   if(m_printheader_counter == 0){
      t_stream << "      System      |        Nodes        |      Bounds     |            Work              \n";
   }
   t_stream << "   time |  memory | explored |     open |  lower |  upper |  vars |   LP it | Pric it | \n";
}

void Statistics::display(std::ostream &t_stream) {
   if(m_printheader_counter % 20 == 0){
      displayHeader(t_stream);
   }
   // TODO: include fractional LB and # solutions, gap% and maximum tree depth?
   //<< std::setw(8) << std::setprecision(4) << timeSinceStart() <<"|"
   auto duration = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() -m_start_solve_time);
   std::string memory = bytesToString(memory_used());
   //TODO: fix almost all fields being written to (OR pass in solver struct data somehow)


   t_stream << std::fixed << std::setw(7) << std::setprecision(1) << std::right << duration << " | "
            << std::setw(7) << std::right<< memory <<" | "
            << std::scientific << std::setw(8) << std::setprecision(2) << std::right << m_num_processed_nodes << " | "
            << std::scientific << std::setw(8) << std::setprecision(2) << std::right << m_num_open_nodes << " | "
            << std::scientific << std::setw(6) << std::setprecision(2) << std::right << m_lb << " | "
            << std::scientific << std::setw(6) << std::setprecision(2) << std::right << m_ub << " | "
            << std::scientific << std::setw(5) << std::setprecision(2) << std::right << m_ub << " | "
            << std::scientific << std::setw(7) << std::setprecision(2) << std::right << m_ub << " | "
            << std::scientific << std::setw(7) << std::setprecision(2) << std::right << m_ub << " | "
            << std::endl;
   ++m_printheader_counter;

}
double Statistics::timeSinceStart() const {

   return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - m_start_solve_time).count();
}

