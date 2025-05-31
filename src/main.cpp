//
// Created by rolf on 25-2-23.
//

#include "git_version.h"
#include <filesystem>
#include <pcog/ColorSolver.hpp>
#include <pcog/GraphIO.hpp>
#include <pcog/config.h>

using namespace pcog;

void printPcogVersion(){
   std::cout<<"PCOG version "<<PCOG_VERSION_MAJOR<<"."<<PCOG_VERSION_MINOR<<"."<<PCOG_VERSION_PATCH<<" [mode: "<<CMAKE_BUILD_TYPE <<"] [git hash: "<<kGitHash<<"]\n";
}

void printVersion(){
   printPcogVersion();
   std::cout<<"Compiled using: "<<__VERSION__<<"\n";
   std::cout<<"External libraries\n";
   std::cout<<"\t SoPlex version: "<<SOPLEX_VERSION_MAJOR<<"."<<SOPLEX_VERSION_MINOR<<"."<<SOPLEX_VERSION_PATCH<<"\n";
   std::cout<<"\t Boost version: "<<BOOST_VERSION / 100000 <<"." << BOOST_VERSION /100 %1000 <<"." << BOOST_VERSION % 100<<"\n";
}
void printUsage(){
   //TODO: make emphasis profiles
   Settings settings;
   std::cout<<"usage: pcog [inputfile] [option [value]]\n";
   std::cout<<"\t -v, --version\t\t\t\t: Prints version and build options, all other options are ignored\n";
   std::cout<<"\t -h, --help\t\t\t\t\t: Prints this help text, all other options are ignored\n";
   std::cout<<"\t -t, --threads\t\t\t\t: (>0, default: "<<settings.numThreads()<< ") The maximum number of threads used. The number of threads used is never greater than the reported number of CPU threads\n";
   std::cout<<"\t -l, --time-limit\t\t\t: (>0, default: none) Time limit in seconds\n";
   std::cout<<"\t -m, --node-limit\t\t\t: (>0, default: none) Stop after processing this number of nodes\n";
   std::cout<<"\t -s, --display-freq\t\t\t: (>=-1, default: "<< settings.nodeDisplayFrequency() << ") The frequency at which a display line is printed. Set to -1 to only print lower/upper bound improvements\n";

   std::cout<<"\t -b, --branch\t\t\t\t: ([0-13], default: "<<static_cast<int>(settings.branchingStrategy())<< ") The branching strategy used to select the branching node pair. See Settings.hpp for more info\n";
   std::cout<<"\t -n, --node-selection\t\t: ([0-2], default: "<<static_cast<int>(settings.nodeSelectionStrategy()) << ") What strategy to use to select the next open node. 0 = pick node with worst lower bound,"
                " 1 = Use DFS until integral lower bound is improved, then pick node with worst lower bound, 2 = similar to 1 but restart after a fixed number of nodes\n";
   std::cout<<"\t --backtrack-freq\t\t\t: ([>0], default: "<<settings.dfsRestartFrequency()<<") Backtracking frequency used when node-selection == 2 (DFS_RESTART)\n";
   std::cout<<"\t -c, --node-child-selection\t: ([0-2], default: "<<static_cast<int>(settings.nodeChildSelectionStrategy()) << ") How to choose child nodes in node selection. 0 = Choose SAME branch, 1 = choose DIFFER branch, 2 = choose randomly\n";
   std::cout<<"\t -p, --pricing \t: ([0-2], default: "<<static_cast<int>(settings.getPricingStrategy()) <<") What pricing algorithms to use. 0 = Only the exact combinatorial pricing algorithm. 1 = A few heuristics + Exact combinatorial algorithm, 2 = Many heuristics + Exact combinatorial algorithm, \n";

   std::cout<<"\t -r, --rounding-tolerance\t: ([0-0.5], default: "<<settings.roundingTolerance()<<") Internal tolerance used to round for heuristics. Does not affect the exact lower bound computation\n";
   std::cout<<"\t -d, --diving-freq\t\t\t: (>=-1, default: "<<settings.divingFrequency() <<") Perform diving at nodes with depth % frequency == 0. Set to -1 to disable, set to zero to only run at the root node\n";
   std::cout<<"\t -f, --diving-pricing-freq\t: (>=-1, default: "<<settings.divingPricingFrequency() <<") Perform pricing during diving at nodes with depth % frequency == 0. Set to -1 to disable, set to zero to only run at the root node\n";
   std::cout<<"\t -i, --init-tabu-iters\t\t: (>=0, default: "<<settings.numInitialTabuIterations() <<") Number of iterations to use during initial tabu search\n";

}
int main(int argc, char ** argv) {
   std::vector<std::string> arguments(argv, argv + argc);
   if (arguments.size() < 2) {
      printPcogVersion();
      printUsage();
      return EXIT_FAILURE;
   }
   if (arguments[1] == "-v" || arguments[1] == "--version") {
      printVersion();
      return EXIT_SUCCESS;
   } else if (arguments[1] == "-h" || arguments[1] == "--help") {
      printPcogVersion();
      printUsage();
      return EXIT_SUCCESS;
   }
   printPcogVersion();
   std::ifstream ifstream(arguments[1]);
   std::optional<DimacsFileHeader> header = readDimacsHeader(ifstream);
   if (!header.has_value()) {
      std::cout << "Could not open or read file!\n";
      return EXIT_FAILURE;
   }

   Settings settings;
   {
      std::size_t i = 2;
      try {
         for (; i < arguments.size(); ++i) {
            std::string &input = arguments[i];
            if (input == "-t" || input == "--threads") {
               ++i;
               std::size_t numThreads = std::stoull(arguments[i]);
               settings.setNumThreads(numThreads);
            } else if (input == "-b" || input == "--branch") {
               ++i;
               int branchInt = std::stoi(arguments[i]);
               if (branchInt < 0 || branchInt > 13)
                  throw std::invalid_argument(arguments[i]);
               settings.setBranchingStrategy(BranchingStrategy(branchInt));
            } else if (input == "-n" || input == "--node-selection") {
               ++i;
               int value = std::stoi(arguments[i]);
               if (value < 0 || value > 2)
                  throw std::invalid_argument(arguments[i]);
               settings.setNodeSelectionStrategy(NodeSelectionStrategy(value));
            }else if (input == "-c" || input == "--node-child-selection") {
               ++i;
               int value = std::stoi(arguments[i]);
               if (value < 0 || value > 2)
                  throw std::invalid_argument(arguments[i]);
               settings.setNodeChildSelectionStrategy(NodeChildSelectionStrategy(value));
            }else if (input == "-p" || input == "--pricing") {
               ++i;
               int value = std::stoi(arguments[i]);
               if (value < 0 || value > 2)
                  throw std::invalid_argument(arguments[i]);
               settings.setPricingAlgorithmStrategy(PricingAlgorithmStrategy(value));
            }else if (input == "-l" || input == "--time-limit") {
               ++i;
               double value = std::stod(arguments[i]);
               settings.setTimeLimit(value);
            }else if (input == "-m" || input == "--node-limit") {
               ++i;
               std::size_t value = std::stoull(arguments[i]);
               settings.setNodeLimit(value);
            }else if (input == "-s" || input == "--display-freq") {
               ++i;
               std::size_t value = std::stoull(arguments[i]);
               settings.setNodeDisplayFrequency(value);
            }else if (input == "-r" || input == "--rounding-tolerance") {
               ++i;
               double value = std::stod(arguments[i]);
               settings.setRoundingTolerance(value);
            }else if (input == "-d" || input == "--diving-freq") {
               ++i;
               int value = std::stoi(arguments[i]);
               if (value < -1)
                  throw std::invalid_argument(arguments[i]);
               settings.setDivingFrequency(value);
            }else if (input == "-f" || input == "--diving-pricing-freq") {
               ++i;
               int value = std::stoi(arguments[i]);
               if (value < -1)
                  throw std::invalid_argument(arguments[i]);
               settings.setDivingPricingFrequency(value);
            }else if(input == "-i" || input == "--init-tabu-iters") {
               ++i;
               std::size_t value = std::stoull(arguments[i]);
               settings.setNumInitialTabuIterations(value);
            }else if (input == "--backtrack-freq") {
               ++i;
               std::size_t value = std::stoull(arguments[i]);
               settings.setDfsRestartFrequency(value);
            }else {
               std::cout << "unrecognized command line option '" << input<< "`\n";
               printUsage();
               return EXIT_FAILURE;
            }
         }
      } catch (std::invalid_argument &arg) {
         std::cout << "Could not parse value: " << arguments[i] << "\n";
         printUsage();
         return EXIT_FAILURE;
      }
   }


   pcog::ColorSolver solver(settings);
   auto graph = header->ReadAsDense();
   if(!graph.has_value()){
      std::cout<<"Could not read graph from file!\n";
      return EXIT_FAILURE;
   }
   std::size_t loops = graph->numSelfLoops();
   if(loops != 0){
      std::cout<<"WARNING: Given graph contains "<<loops<<" loop edges. These are being ignored by the solver.\n";
      graph->removeSelfLoops();
   }
   assert(graph->numSelfLoops() == 0);
   std::string strippedName = std::filesystem::path(arguments[1]).stem();
   solver.setProblem(strippedName,graph.value());
   solver.solve();
   return EXIT_SUCCESS;
}