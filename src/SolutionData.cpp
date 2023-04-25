//
// Created by rolf on 7-4-23.
//

#include "pcog/SolutionData.hpp"

#include <iostream>
#include <iomanip>
#include <array>
#include <utility>

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


namespace pcog {
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
      fscanf (file, "%lu", &vm);  // Just need the first num: vm size
      fclose (file);
      size = static_cast<size_t>(vm) * getpagesize();
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
std::string itsToString(std::size_t iterations) {
   static constexpr std::array<const char *, 4> extension = {"", "k", "M", "B"};
   std::size_t extension_ind = 0;
   double show = iterations;
   while(extension_ind < 4 && show > 1'000.0){
      show /= 1000.0;
      ++extension_ind;
   }
   std::stringstream stream;
   stream << std::setprecision(3) << show;
   std::string str = stream.str() + extension[extension_ind];
   return str;
}

void SolutionData::addSolution(std::vector<std::size_t> t_stable_set_indices) {
#ifndef NDEBUG
   {
      // assert that coloring is indeed a valid coloring
      DenseSet coveredNodes(m_preprocessedGraph.numNodes());
      for (const auto &index : t_stable_set_indices) {
         assert(index < m_variables.size());
         coveredNodes.inplaceUnion(m_variables[index].set());
      }
      assert(coveredNodes.full());
   }
#endif

   std::size_t ub = t_stable_set_indices.size();
   if (ub < m_upperBound) {
      // prune redundant nodes from the tree
      m_tree.pruneUpperBound(ub);
      m_upperBound = ub;
      m_incumbent_index = m_colorings.size();
      display(std::cout);
   }
   m_colorings.emplace_back(std::move(t_stable_set_indices));
}
std::size_t SolutionData::upperBound() const { return m_upperBound; }
std::size_t SolutionData::lowerBound() const { return m_lowerBound; }
double SolutionData::fractionalLowerBound() const {
   return m_fractionalLowerBound;
}
bool SolutionData::isNewSet(const DenseSet &t_set) const {
   return std::all_of(m_variables.begin(), m_variables.end(),
                      [&](const StableSetVariable &variable) {
                         return variable.set() != t_set;
                      });
}

std::size_t SolutionData::addStableSet(pcog::DenseSet t_set) {
   assert(isNewSet(t_set));
   assert(m_preprocessedGraph.setIsStable(t_set));
   std::size_t index = m_variables.size();
   m_variables.emplace_back(std::move(t_set));
   return index;
}
std::size_t SolutionData::findOrAddStableSet(const DenseSet &t_set) {
   for (std::size_t i = 0; i < m_variables.size(); ++i) {
      if (m_variables[i].set() == t_set) {
         return i;
      }
   }
   return addStableSet(t_set);
}
SolutionData::SolutionData( Settings& t_settings)
    :m_incumbent_index{std::numeric_limits<std::size_t>::max()},
      m_upperBound{std::numeric_limits<std::size_t>::max()},
      m_lowerBound{0},
      m_fractionalLowerBound{0.0},
      m_printheader_counter{0},
      m_settings{t_settings},
      m_lpIterations{0},
      m_pricingIterations{0}
      {}

void SolutionData::reset(DenseGraph t_graph) {
   m_originalGraph = std::move(t_graph);
   m_preprocessedGraph.clear();
   m_preprocessedToOriginal.clear();
   m_variables.clear();
   m_colorings.clear();
   m_tree.clear();

   m_incumbent_index = std::numeric_limits<std::size_t>::max();
   m_upperBound = std::numeric_limits<std::size_t>::max();
   m_lowerBound = 0;
   m_fractionalLowerBound = 0.0 ;
   m_printheader_counter = 0;
   m_lpIterations = 0;
   m_pricingIterations = 0;
}
const DenseGraph &SolutionData::originalGraph() const {
   return m_originalGraph;
}
void SolutionData::displayHeader(std::ostream& t_stream) const{
   if(m_printheader_counter == 0){
      t_stream << "      System      |        Nodes        |          Bounds          |            Work              \n";
   }
   t_stream << "   time |  memory | explored |     open |   frac |  lower |  upper |  vars |   LP it | Pric it | \n";
}
void SolutionData::display(std::ostream& t_stream){
   if(m_printheader_counter % 20 == 0){
      displayHeader(t_stream);
   }
   // TODO: include fractional LB and # solutions, gap% and maximum tree depth?
   //<< std::setw(8) << std::setprecision(4) << timeSinceStart() <<"|"
   auto duration = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() -m_start_solve_time);
   std::string memory = bytesToString(memory_used());
   //TODO: fix almost all fields being written to (OR pass in solver struct data somehow)


   //TODO: ideally we would like to just use cout<< duration but it seems that this is not supported yet by all compilers
   t_stream << std::fixed << std::setw(7) << std::setprecision(1) << std::right << duration.count() << " | "
            << std::setw(7) << std::right<< memory <<" | "
            << std::scientific << std::setw(8) << std::setprecision(2) << std::right << numProcessedNodes() << " | "
            << std::scientific << std::setw(8) << std::setprecision(2) << std::right << numOpenNodes() << " | "
            << std::fixed << std::setw(6) << std::setprecision(2) << std::right << m_fractionalLowerBound << " | "
            << std::scientific << std::setw(6) << std::setprecision(2) << std::right << m_lowerBound << " | "
            << std::scientific << std::setw(6) << std::setprecision(2) << std::right << m_upperBound << " | "
            << std::scientific << std::setw(5) << std::setprecision(2) << std::right << m_variables.size() << " | "
            << std::scientific << std::setw(7) << std::setprecision(2) << std::right << itsToString(m_lpIterations) << " | "
            << std::scientific << std::setw(7) << std::setprecision(2) << std::right << itsToString(m_pricingIterations) << " | "
            << std::endl;
   ++m_printheader_counter;
}
double SolutionData::timeSinceStart() const {
   return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - m_start_solve_time).count();
}
bool SolutionData::checkTimelimitHit() const {
   double timeLimit = m_settings.timeLimit();
   if (timeLimit == NO_TIME_LIMIT) {
      return false;
   }
   return timeSinceStart() >= timeLimit;
}

void SolutionData::startSolveTime() {
   m_start_solve_time = std::chrono::high_resolution_clock::now();
}
void SolutionData::doPresolve() {
   // TODO: find an initial coloring using a simple greedy method.

   // This can help to judge if it's 'worth' to find a clique for the lower
   // bound; If the minimal degree of the graph is larger than the chromatic
   // number (which typically happens for dense graphs), then the low-degree
   // preprocessing rule will not find anything unless some other rule finds a
   // reduction which lowers a node degree of the graph below the chromatic
   // number
   auto result = preprocessOriginalGraph(m_originalGraph);
   m_preprocessedGraph = result.graph;
   m_preprocessedToOriginal = result.map;
   // Also add coloring variables as initial solution
}
bool SolutionData::checkNodeLimitHit() const {
   return m_tree.numProcessedNodes() >= m_settings.nodeLimit();
}
void SolutionData::initializeBBTree() {
   m_tree.createRootNode(m_preprocessedGraph.numNodes());
}
bool SolutionData::hasOpenNodes() const {
   return m_tree.hasOpenNodes();
}
std::size_t SolutionData::numOpenNodes() const {
   return m_tree.numOpenNodes();
}
std::size_t SolutionData::numProcessedNodes() const {
   return m_tree.numProcessedNodes();
}
const std::vector<StableSetVariable> &SolutionData::variables() const {
   return m_variables;
}
const Settings &SolutionData::settings() const { return m_settings; }
BBNode && SolutionData::popNextNode() {
   return m_tree.popNextNode();
}
void SolutionData::createChildren(const BBNode& t_node,
                                  ColorNodeWorker &t_nodeWorker) {
   m_tree.createChildren(t_node,t_nodeWorker);
}
const DenseGraph &SolutionData::preprocessedGraph() const {
   return m_preprocessedGraph;
}
const PreprocessedMap &SolutionData::preprocessingMap() const {
   return m_preprocessedToOriginal;
}
const NodeMap &SolutionData::preprocessedToOriginal() const {
   return m_preprocessedToOriginal.newToOldIDs;
}
const NodeMap &SolutionData::originalToPreprocessed() const {
   return m_preprocessedToOriginal.oldToNewIDs;
}
void SolutionData::addPricingIterations(std::size_t count) {
   m_pricingIterations += count;
}
void SolutionData::addLPIterations(std::size_t count) {
   m_lpIterations += count;
}
void SolutionData::updateLowerBound(std::size_t t_lb) {
   if(t_lb > m_lowerBound){
      m_lowerBound = t_lb;
      //TODO: are there any lb changes we need to trigger here?
      display(std::cout);
   }
}
void SolutionData::updateFractionalLowerBound(double t_fractional_lb) {
   if(t_fractional_lb > m_fractionalLowerBound){
      m_fractionalLowerBound = t_fractional_lb;
   }
}
void SolutionData::updateTreeBounds() {
   if(m_tree.hasOpenNodes()) {
      updateLowerBound(m_tree.lowerBound());
      updateFractionalLowerBound(m_tree.fractionalLowerBound());
   }else {
      updateLowerBound(m_upperBound);
      updateFractionalLowerBound(m_upperBound);
   }

}
}// namespace pcog