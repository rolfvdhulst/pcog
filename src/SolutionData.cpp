//
// Created by rolf on 7-4-23.
//

#include "pcog/SolutionData.hpp"
#include "pcog/GreedyColoring.hpp"
#include "pcog/StableSetMaximizer.hpp"
#include "pcog/TabuColoring.hpp"

#include <array>
#include <iomanip>
#include <iostream>
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


#include <thread>
#include "pcog/reduction/Reducer.hpp"

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

void SolutionData::addSolutionPresolve(std::vector<std::size_t> t_stable_set_indices) {
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
std::size_t SolutionData::lowerBoundUnscaled() {
   std::scoped_lock guard(m_lowerBound_mutex);
   return m_lowerBound;}
std::size_t SolutionData::upperBoundUnscaled() {
   std::scoped_lock guard(m_upperBound_mutex);
   return m_upperBound;
}

std::size_t SolutionData::upperBound() {
   return upperBoundUnscaled() + m_preprocessingStack.numFixedColors();
}
std::size_t SolutionData::lowerBound() {
   return lowerBoundUnscaled() + m_preprocessingStack.numFixedColors();
}
double SolutionData::fractionalLowerBound()  {
   std::scoped_lock guard(m_lowerBound_mutex);
   return m_fractionalLowerBound + static_cast<double>(m_preprocessingStack.numFixedColors());
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
   m_preprocessedToOriginalMap.clear();
   m_originalToPreprocessedMap.clear();
   m_preprocessingStack.clear();
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

//Estimates the difficulty of solving the pricing problem exactly by roughly
//estimating the b&b tree size of a simple heuristic (using quite a few assumptions)
double difficultyEstimation(const DenseGraph& graph){
   double invDensity = 1.0- graph.density();
   double N = graph.numNodes();
   double avDepth = - std::log(N) / std::log(invDensity);
   double sum = avDepth < 1.0 ? 1.0 : 0.0;
   for(std::size_t i = 0; i < (avDepth); ++i){
      sum += std::pow(invDensity, static_cast<double>(i * (i + 1)) /2.0);
   }
   sum *= N;

   //the function tends to slightly overestimate the difficulty of dense graphs; hence we 'correct' it
   sum = sum*(1+invDensity);
   return sum;
}

void SolutionData::displayHeader(std::ostream& t_stream) const{
   if(m_printheader_counter == 0){
      t_stream << "      System      |        Nodes        |          Bounds          |            Work              \n";
   }
   t_stream << "   time |  memory | explored |     open |   frac |  lower |  upper |  vars |   LP it | Pric it | \n";
}
void SolutionData::display(std::ostream& t_stream){
   std::scoped_lock guard(m_printing_mutex);
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
            << std::defaultfloat << std::setw(6) << std::setprecision(5) << std::right << fractionalLowerBound() << " | "
            << std::scientific << std::setw(6) << std::setprecision(2) << std::right << lowerBound() << " | ";
   if(upperBound() == std::numeric_limits<std::size_t>::max()) {
      t_stream << "    -  | ";
   }else {
      t_stream << std::scientific << std::setw(6) << std::setprecision(2) << std::right << upperBound() << " | ";
   }
   t_stream << std::scientific << std::setw(5) << std::setprecision(2) << std::right << m_variables.size() << " | "
            << std::scientific << std::setw(7) << std::setprecision(2) << std::right << itsToString(m_lpIterations) << " | "
            << std::scientific << std::setw(7) << std::setprecision(2) << std::right << itsToString(m_pricingIterations) << " | "
            << std::endl;
   ++m_printheader_counter;
}
std::chrono::duration<double> SolutionData::timeSinceStart() const {
   return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - m_start_solve_time);
}


void SolutionData::startSolveTime() {
   m_start_solve_time = std::chrono::high_resolution_clock::now();
}
void SolutionData::doPresolve() {
   // TODO: find an initial coloring using a simple greedy method.
   GreedyColoring greedy(m_originalGraph); //TODO: extract clique from saturation degree
   auto coloring = greedy.run_saturation_degree();
   std::cout<<"DSATUR found " << coloring.numColors()<<"-coloring\n";

   std::cout<<"Original graph has " <<m_originalGraph.numNodes()<<" nodes, density: "<< m_originalGraph.density()*100.0<<"%\n";
   auto result = reduceGraph(m_originalGraph);
   {
      std::size_t nSimplicial = 0;
      std::size_t nSimplicialNodes = 0;
      std::size_t nLowDegree = 0;
      std::size_t nDominated = 0;
      std::size_t nFoldDegreeTwo = 0;
      std::size_t nTwinDegreeThree = 0;
      std::size_t nTwinDegreeThreeNodes = 0;
      std::size_t nTwinDegreeThreeFold = 0;
      std::size_t nCrown = 0;
      std::size_t nCrownNodes = 0;
      std::size_t nCrownColors = 0;
      std::size_t nTwoFixing = 0;
      std::size_t nTwoFixingNodes = 0;

      //TODO: record this data during presolve itself rather than going through the stack here
      for(const auto& reduction : result.stack.reductionStack()) {
         if(std::holds_alternative<LowDegreeReduction>(reduction)) {
            nLowDegree += 1;
         }
         else if(std::holds_alternative<DominatedReduction>(reduction)) {
            nDominated += 1;
         }
         else if(std::holds_alternative<SimplicialReduction>(reduction)){
            nSimplicial += 1;
            auto reduc = std::get<SimplicialReduction>(reduction);
            nSimplicialNodes += reduc.set.size();
         }
         else if(std::holds_alternative<FoldDegreeTwoReduction>(reduction)){
            nFoldDegreeTwo += 1;
         }
         else if(std::holds_alternative<TwinDegreeThreeReduction>(reduction)){
            nTwinDegreeThree += 1;
            auto reduc = std::get<TwinDegreeThreeReduction>(reduction);
            nTwinDegreeThreeNodes += reduc.uSet.size();
            nTwinDegreeThreeNodes += reduc.vSet.size();
         }
         else if(std::holds_alternative<TwinDegreeThreeFoldReduction>(reduction)) {
            nTwinDegreeThreeFold += 1;
         }
         else if(std::holds_alternative<CrownReduction>(reduction)){
            nCrown += 1;
            const auto& crownReduction = std::get<CrownReduction>(reduction);
            nCrownColors += crownReduction.fixedSets.size();
            for(const auto& set : crownReduction.fixedSets) {
               nCrownNodes += set.size();
            }
         }else if(std::holds_alternative<TwoFixingReduction>(reduction)){
            nTwoFixing += 2;

            const auto& reduc = std::get<TwoFixingReduction>(reduction);
            nTwoFixingNodes += reduc.firstSet.size();
            nTwoFixingNodes += reduc.secondSet.size();
         }
      }
      std::cout << "Type          | Reductions | removed nodes | fixed colors\n";
      std::cout << "Low degree    | "<<std::setw(10) << nLowDegree <<" | "<< std::setw(13) << nLowDegree << " |            0\n";
      std::cout << "Dominated     | "<<std::setw(10) << nDominated <<" | "<< std::setw(13) << nDominated << " |            0\n";
      std::cout << "Simplicial    | "<<std::setw(10) << nSimplicial <<" | "<< std::setw(13) <<nSimplicialNodes<<" | "<< std::setw(12)<< nSimplicial<<"\n";
      std::cout << "2-folding     | "<<std::setw(10) << nFoldDegreeTwo <<" | "<< std::setw(13) <<2*nFoldDegreeTwo<<" | "<< std::setw(12)<< nFoldDegreeTwo<<"\n";
      std::cout << "Twin degree 3 | "<<std::setw(10) << nTwinDegreeThree <<" | "<< std::setw(13) << nTwinDegreeThreeNodes <<" | "<< std::setw(12)<< 2*nTwinDegreeThree<<"\n";
      std::cout << "Twin 3-folding| "<<std::setw(10) << nTwinDegreeThreeFold <<" | "<< std::setw(13) <<4*nTwinDegreeThreeFold << " |            0\n";
      std::cout << "Crown         | "<<std::setw(10) << nCrown <<" | "<< std::setw(13) <<nCrownNodes << " | "<< std::setw(12) << nCrownColors << "\n";
      std::cout << "2-fixing      | "<<std::setw(10) << nTwoFixing <<" | "<<std::setw(13) << nTwoFixingNodes <<" | "<<std::setw(12) << 2*nTwoFixing << "\n";

      using Reduction = std::variant<LowDegreeReduction,DominatedReduction,SimplicialReduction,
                                     FoldDegreeTwoReduction,TwinDegreeThreeReduction,TwinDegreeThreeFoldReduction,CrownReduction,
                                     TwoFixingReduction>;
   }
   std::cout<<"Presolved graph has "<<result.graph.graph.numNodes()<<" nodes, density: "<< result.graph.graph.density()*100.0<<"%\n";

   m_preprocessedGraph = result.graph.graph;
   m_preprocessingStack = result.stack;
   m_preprocessedToOriginalMap = result.graph.newToOld;
   m_originalToPreprocessedMap = result.graph.oldToNew;
   std::size_t lb = result.lowerBoundClique.size();
   assert(m_preprocessedGraph.setIsClique(result.lowerBoundClique));


   //Project coloring to preprocessed graph
   SetColoring projectedColoring;
   for(const auto& color : coloring.colors()) {
      DenseSet set(m_preprocessedGraph.numNodes());
      m_originalToPreprocessedMap.transform(color,set);
      if(!set.empty()) {
         projectedColoring.addColor(set);
      }
   }
   //TODO: re-run DSATUR, maybe?
   NodeColoring initialColoring(m_preprocessedGraph.numNodes(),projectedColoring);


   //Run tabu search to find a better initial solution.
   std::size_t numSearchColors = initialColoring.numColors()-1;
   NodeColoring searchColoring = initialColoring;
   NodeColoring bestColoring = initialColoring;
   TabuColoring tabuAlgorithm(m_preprocessedGraph);
   tabuAlgorithm.setMaxIterations(settings().numInitialTabuIterations());
   while( lb <= numSearchColors) {
      std::size_t removeColor = numSearchColors; //We remove the 'highest' color, but different strategies are possible
      searchColoring.setNumColors(numSearchColors);
      for(std::size_t i = 0; i < m_preprocessedGraph.numNodes(); ++i) {
         if(searchColoring[i] == removeColor) {
            searchColoring[i]--; // We set it to be the next color; again different strategies are possible.
         }
      }
      auto foundColoring = tabuAlgorithm.run(searchColoring);
      if(foundColoring.has_value()) {
         searchColoring = foundColoring.value();
         bestColoring = foundColoring.value();
         std::cout<<"Tabu coloring found "<<bestColoring.numColors()<<"-coloring after "<<tabuAlgorithm.numIterations()<<" iterations\n";
      }else {
         break;
      }
      --numSearchColors;
   }
   //Add best coloring to the LP
   StableSetMaximizer maximizer(42);
   SetColoring bestSetColoring(bestColoring);
   for(auto& color : bestSetColoring.colors()) {
      maximizer.maximizeRandomly(color,m_preprocessedGraph);
   }
   std::vector<std::size_t> indices;
   for(const auto& color : bestSetColoring.colors()) {
      std::size_t index = addStableSet(color);
      indices.push_back(index);
   }

   if(lb > 0) {
      auto f_lb =  static_cast<double>(lb) ; //TODO: rounding mode / safe rounding

      std::scoped_lock guard(m_lowerBound_mutex);
      assignLB(f_lb,lb);
   }
   addSolutionPresolve(indices);

}
bool SolutionData::checkNodeLimitHit() const {
   return m_tree.numProcessedNodes() >= m_settings.nodeLimit();
}
void SolutionData::initializeBBTree() {
   m_tree.createRootNode(m_preprocessedGraph.numNodes());
}

//Nodes count as explored when branch-and-price loop for that node was terminated
//Pruned nodes are nodes which are created but for which branch-and-price is never run
//Open nodes are all nodes which are not yet explored or pruned
std::size_t SolutionData::numOpenNodes() {
   std::scoped_lock guard(m_lowerBound_mutex);
   std::size_t numNodes = m_tree.numOpenNodes();
   for(const auto& lb : m_processing_node_lower_bounds) {
      if(lb.has_value()) {
         numNodes++;
      }
   }
   return numNodes;
}
std::size_t SolutionData::numProcessedNodes() {
   std::scoped_lock guard(m_lowerBound_mutex);
   std::size_t numNodes = m_tree.numProcessedNodes();
   for(const auto& lb : m_processing_node_lower_bounds) {
      if(lb.has_value()) {
         --numNodes;
      }
   }
   return numNodes;
}
const std::vector<StableSetVariable> &SolutionData::variables() const {
   return m_variables;
}
const Settings &SolutionData::settings() const { return m_settings; }
std::vector<node_id> SolutionData::createChildren(const BBNode& t_node,
                                  ColorNodeWorker &t_nodeWorker) {
   return m_tree.createChildren(t_node,t_nodeWorker);
}
const DenseGraph &SolutionData::preprocessedGraph() const {
   return m_preprocessedGraph;
}

const NodeMap &SolutionData::preprocessedToOriginal() const {
   return m_preprocessedToOriginalMap;
}
const NodeMap &SolutionData::originalToPreprocessed() const {
   return m_originalToPreprocessedMap;
}


SetColoring SolutionData::incumbentUnscaled() const {
   SetColoring coloring;
   for(const auto& index : m_colorings[m_incumbent_index]) {
      coloring.addColor(m_variables[index].set());
   }
   return coloring;
}
NodeColoring SolutionData::incumbent() const {
   SetColoring preprocessedSol = incumbentUnscaled();
   NodeColoring preprocessedCol(m_preprocessedGraph.numNodes(),preprocessedSol);

   NodeColoring originalColoring(m_originalToPreprocessedMap.size());
   for(std::size_t i = 0; i < preprocessedCol.numNodes(); ++i) {
      originalColoring[m_preprocessedToOriginalMap[i]] = preprocessedCol[i];
   }
   originalColoring.setNumColors(preprocessedCol.numColors());
   m_preprocessingStack.newToOldColoring(originalColoring);
#ifndef NDEBUG
   for(std::size_t i = 0; i < originalColoring.numNodes(); ++i) {
      assert(originalColoring[i] != INVALID_COLOR);
      for(const auto& neighbour : m_originalGraph.neighbourhood(i)) {
         assert(originalColoring[i] != originalColoring[neighbour]);
      }
   }
#endif
   return originalColoring;
}

void SolutionData::synchronizeLocalDataStatistics(LocalSolutionData &t_localData) {
   m_lpIterations += t_localData.m_lpIterations;
   t_localData.m_lpIterations = 0;

   m_pricingIterations += t_localData.m_pricingIterations;
   t_localData.m_pricingIterations = 0;
}

void SolutionData::writeLocalVariablesToGlobal(
    LocalSolutionData &t_localSolutionData) {
   //This function assumes that the user has locked the mutex; it is just there to avoid repetition
   assert(t_localSolutionData.m_lastGlobalAddedIndex >= t_localSolutionData.m_numStartGlobalVariables);
   assert(t_localSolutionData.m_lastGlobalAddedIndex - t_localSolutionData.m_numStartGlobalVariables ==
          t_localSolutionData.m_variable_mapping.size());

   if(t_localSolutionData.m_variables.size() == t_localSolutionData.m_numStartGlobalVariables){
      return;
   }
   assert(t_localSolutionData.m_variables.size() > t_localSolutionData.m_numStartGlobalVariables);
   t_localSolutionData.m_variable_mapping.resize(t_localSolutionData.m_variables.size() - t_localSolutionData.m_numStartGlobalVariables);
   for(std::size_t i = t_localSolutionData.m_lastGlobalAddedIndex; i < t_localSolutionData.m_variables.size(); ++i) {
      std::size_t index = findOrAddStableSet(t_localSolutionData.m_variables[i].set(),t_localSolutionData.m_numStartGlobalVariables);
      t_localSolutionData.m_variable_mapping[i- t_localSolutionData.m_numStartGlobalVariables] = index;
   }
   t_localSolutionData.m_lastGlobalAddedIndex = t_localSolutionData.m_variable_mapping.size() + t_localSolutionData.m_numStartGlobalVariables;
}
std::size_t SolutionData::findOrAddStableSet(const DenseSet &t_set,
                                             std::size_t checkNewFromIndex) {
#ifndef NDEBUG
   for(std::size_t i = 0; i < checkNewFromIndex; ++i) {
      assert(m_variables[i].set() != t_set);
   }
#endif
   for (std::size_t i = checkNewFromIndex; i < m_variables.size(); ++i) {
      if (m_variables[i].set() == t_set) {
         return i;
      }
   }
   return addStableSet(t_set);
}
void SolutionData::syncLocalVarsWithGlobal(
    LocalSolutionData &t_localSolutionData) {
   std::scoped_lock guard(m_variable_mutex);
   // First, add the local variables to the global pool.
   writeLocalVariablesToGlobal(t_localSolutionData);

// The variables up until t_localSolutionData.m_numStartGlobalVariables should
// be in the global data on the same indices already
#ifndef NDEBUG
   for (std::size_t i = 0; i < t_localSolutionData.m_numStartGlobalVariables;
        ++i) {
      assert(m_variables[i].set() == t_localSolutionData.m_variables[i].set());
   }
#endif
   // To prevent excessive copies, we copy only the variables which (may have)
   // changed
   t_localSolutionData.m_variables.resize(
       t_localSolutionData.m_numStartGlobalVariables);
   t_localSolutionData.m_variables.insert(
       t_localSolutionData.m_variables.end(),
       m_variables.begin() + t_localSolutionData.m_numStartGlobalVariables,
       m_variables.end());
// Check if the update was indeed correct
#ifndef NDEBUG
   assert(m_variables.size() == t_localSolutionData.m_variables.size());
   for (std::size_t i = t_localSolutionData.m_numStartGlobalVariables;
        i < m_variables.size(); ++i) {
      assert(m_variables[i].set() == t_localSolutionData.m_variables[i].set());
   }
#endif

   // Clear the residual data
   t_localSolutionData.m_numStartGlobalVariables = t_localSolutionData.m_variables.size();
   t_localSolutionData.m_lastGlobalAddedIndex =  t_localSolutionData.m_numStartGlobalVariables;
   t_localSolutionData.m_variable_mapping.clear();
}
void SolutionData::writeLocalVarsToGlobal(
    LocalSolutionData &t_localSolutionData, std::atomic_bool& stop) {
   // Synchronize variables
   {
      std::scoped_lock guard(m_variable_mutex);
      writeLocalVariablesToGlobal(t_localSolutionData);
   }
   // Transform local solutions to global ones and add them
   writeLocalSolutionsToGlobal(t_localSolutionData,stop);
}

void SolutionData::syncLocalLowerBound(LowerBoundInfo lbInfo, std::size_t worker_id, std::atomic_bool& stop) {
   bool improved = false;
   {
      std::scoped_lock guard(m_lowerBound_mutex);
      assert(m_processing_node_lower_bounds[worker_id].has_value());
      if(m_processing_node_lower_bounds[worker_id]->m_fractional_lb < lbInfo.m_fractional_lb) {
         m_processing_node_lower_bounds[worker_id]->m_fractional_lb = lbInfo.m_fractional_lb;
      }
      if(m_processing_node_lower_bounds[worker_id]->m_lb < lbInfo.m_lb) {
         m_processing_node_lower_bounds[worker_id]->m_lb = lbInfo.m_lb;
      }
      improved = recomputeLowerBound();
      std::scoped_lock ub_guard(m_upperBound_mutex);
      if(m_lowerBound == m_upperBound) {
         stop = true;
         for(auto& worker : m_workers) {
            worker.cancelNode(true);
         }
      }
   }
   if(improved) {
      display(std::cout);
   }
}
bool SolutionData::recomputeLowerBound() {
   // recompute the global lower bounds
   double m_tree_f_lb = m_tree.fractionalLowerBound();
   std::size_t m_tree_lb = m_tree.lowerBound();
   // also take into account the LB's of the nodes which are currently processing
   for (const auto &info : m_processing_node_lower_bounds) {
      if (info.has_value()) {
         if (info->m_fractional_lb < m_tree_f_lb) {
            m_tree_f_lb = info->m_fractional_lb;
         }
         if (info->m_lb < m_tree_lb) {
            m_tree_lb = info->m_lb;
         }
      }
   }
   return assignLB(m_tree_f_lb,m_tree_lb);
}
bool SolutionData::assignLB(double t_frac_lb, std::size_t t_lb) {
   bool improved = false;
   //B&B tree has no open nodes and no node is being processed
   if(t_lb == std::numeric_limits<std::size_t>::max()) {
      std::size_t upperBound = m_upperBound; //TODO: is this safe?
      improved = m_lowerBound < upperBound;
      m_fractionalLowerBound = upperBound ; //TODO: prevent overwriting fractional bound
      m_lowerBound = upperBound;
      return improved;
   }
   if (t_lb > m_lowerBound) {
      m_lowerBound = t_lb;
      improved = true;
   }
   m_fractionalLowerBound = std::max(m_fractionalLowerBound, t_frac_lb);
   return improved;
}
void SolutionData::runBranchAndBound() {
   initializeBBTree();
   std::size_t numHardWareThreads = std::min(settings().numThreads(),static_cast<std::size_t>(std::thread::hardware_concurrency()));
   m_workers.clear();
   for(std::size_t id = 0; id < numHardWareThreads; ++id) {
      m_workers.emplace_back(id);
   }

   m_processing_node_lower_bounds.resize(numHardWareThreads,std::nullopt);
   std::atomic_bool stop_jobs = false;

   std::vector<std::thread> threads;
   for(std::size_t id = 0; id < numHardWareThreads; ++id) {
      threads.emplace_back(&ColorNodeWorker::runLoop,&m_workers[id],std::ref(*this),std::ref(stop_jobs));
   }
   for(auto& thread : threads) {
      thread.join();
   }

}
std::optional<BBNode> SolutionData::popNextNode(ColorNodeWorker &t_nodeWorker) {
   std::scoped_lock guard(m_lowerBound_mutex);
   if(!m_tree.hasOpenNodes()) {
      return std::nullopt;
   }
   //As this worker is currently idling, its previous node must not have had any children (otherwise the B&B tree would have had nodes available for it)
   BBNode node = m_tree.popNextNode();
   assert( !m_processing_node_lower_bounds[t_nodeWorker.id()].has_value());
   m_processing_node_lower_bounds[t_nodeWorker.id()] = LowerBoundInfo(node.fractionalLowerBound(),node.lowerBound());
   t_nodeWorker.cancelNode(false);
   return node;


}

std::optional<BBNode>
SolutionData::branchAndPopNode(BBNode& t_node, ColorNodeWorker &t_nodeWorker) {
   std::vector<node_id> children;
   std::scoped_lock guard(m_lowerBound_mutex);
   auto& basis = t_node.mutableBasis();
   for(std::size_t& index : basis.basicCols){
      index = t_nodeWorker.mapToGlobalIndex(index);
   }
   if (t_node.status() == BBNodeStatus::BRANCHED) {
      children = createChildren(t_node,t_nodeWorker);
   }
   std::optional<BBNode> node = pickNextNode(t_node,t_nodeWorker,children);
   if(node.has_value()) {
      m_processing_node_lower_bounds[t_nodeWorker.id()] = LowerBoundInfo(node->fractionalLowerBound(),node->lowerBound());
   }else {
      m_processing_node_lower_bounds[t_nodeWorker.id()] = std::nullopt;
   }
   return node;
}
std::optional<BBNode> SolutionData::pickNextNode(BBNode &t_node,
                                    ColorNodeWorker &t_nodeWorker,
                                    const std::vector<node_id> &children) {
   if(!m_tree.hasOpenNodes()) {
      return std::nullopt;
   }
   switch (m_settings.nodeSelectionStrategy()) {
   default:
   case ::NodeSelectionStrategy::BOUND:{
      return m_tree.popNextNode();
   }
   case ::NodeSelectionStrategy::DFS_BOUND:{
      if(!children.empty() && m_lowerBound == t_node.lowerBound()) {
         std::size_t nodePos = t_nodeWorker.pickChildNode(m_settings.nodeChildSelectionStrategy(),children);
         return m_tree.popNodeWithID(nodePos);
      }
      return m_tree.popNextNode();
   }
   case ::NodeSelectionStrategy::DFS_RESTART:
   {
      if(!children.empty() && t_nodeWorker.successiveChildrenProcessed() < m_settings.dfsRestartFrequency()) {
         std::size_t nodePos = t_nodeWorker.pickChildNode(m_settings.nodeChildSelectionStrategy(),children);
         return m_tree.popNodeWithID(nodePos);
      }
      return m_tree.popNextNode();
   }
   }
}
bool SolutionData::doRecomputeLowerBound(std::atomic_bool& stop) {
   bool improved = false;
   {
      std::scoped_lock guard(m_lowerBound_mutex);
      improved = recomputeLowerBound();

      std::scoped_lock ub_guard(m_upperBound_mutex);
      if(improved && m_lowerBound == m_upperBound) {
         stopComputation(stop);
      }
   }
   if (improved) {
      display(std::cout);

   }
   return improved;
}
void SolutionData::pruneUpperBound(std::size_t upperBound, std::atomic_bool& stop) {
   std::scoped_lock guard(m_lowerBound_mutex);
   m_tree.pruneUpperBound(upperBound);
   for(std::size_t i = 0; i < m_processing_node_lower_bounds.size(); ++i) {
      if(m_processing_node_lower_bounds[i].has_value() &&
          m_processing_node_lower_bounds[i]->m_lb >= upperBound) {
         m_workers[i].cancelNode(true);
      }
   }
   if(m_lowerBound == upperBound) {
      stop = true;
   }

}
void SolutionData::writeLocalSolutionsToGlobal(
    LocalSolutionData &t_localSolutionData, std::atomic_bool &stop) {
   if (t_localSolutionData.m_solutions.empty()) {
      return;
   }
   for (auto &solution : t_localSolutionData.m_solutions) {
      for (auto &var_index : solution) {
         if (var_index >= t_localSolutionData.m_numStartGlobalVariables) {
            var_index = t_localSolutionData.m_variable_mapping[var_index - t_localSolutionData.m_numStartGlobalVariables];
         }
      }
   }
   // move solutions to global solution pool

   std::size_t minSize = t_localSolutionData.m_solutions[0].size();
   std::size_t minColIndex = m_colorings.size();
   bool incumbentChanged = false;
   {
      std::scoped_lock guard(m_upperBound_mutex);
      for (const auto &coloring : t_localSolutionData.m_solutions) {
         if (coloring.size() < minSize) {
            minSize = coloring.size();
            minColIndex = m_colorings.size();
         }
         m_colorings.push_back(coloring);
      }
      if (minSize < m_upperBound) {
         m_upperBound = minSize;
         m_incumbent_index = minColIndex;
         incumbentChanged = true;
         std::scoped_lock lb_guard(m_lowerBound_mutex);
         if(m_lowerBound == m_upperBound) {
            stopComputation(stop);
         }
      }
      t_localSolutionData.m_solutions.clear();
   }
   if (incumbentChanged) {
      pruneUpperBound(minSize,stop);
      display(std::cout);
   }

}
std::optional<std::chrono::duration<double>> SolutionData::timeLeft() const {
   if(m_settings.timeLimit() == NO_TIME_LIMIT) {
      return std::nullopt;
   }
   auto timeLimit = std::chrono::duration<double>(m_settings.timeLimit());
   return timeLimit - std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - m_start_solve_time);
}

bool SolutionData::checkTimeLimitHit(std::atomic_bool &stop) {
   auto durationLeft = timeLeft();
   bool timeLimitHit = durationLeft.has_value() && durationLeft.value().count() < TIMEOUT_SAFETY;
   if(timeLimitHit) {
      stopComputation(stop);
   }
   return timeLimitHit;
}
void SolutionData::stopComputation(std::atomic_bool &stop) {
   if(!stop) {
      stop = true;
      for (auto &worker : m_workers) {
         worker.cancelNode(true);
      }
   }
}
}// namespace pcog