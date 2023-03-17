//
// Created by rolf on 25-2-23.
//

#include <pcog/GraphIO.hpp>
#include <pcog/ColorSolver.hpp>

using namespace pcog;
int main(int argc, char ** argv){
   std::vector<std::string> arguments(argv,argv+argc);
   std::ifstream ifstream(arguments[1]);
   std::optional<DimacsFileHeader> header = readDimacsHeader(ifstream);
   if(!header.has_value()){
      std::cout<<"Could not open or read file!\n";
      return EXIT_FAILURE;
   }
   pcog::ColorSolver solver;
   auto graph = header->ReadAsDense();
   if(!graph.has_value()){
      std::cout<<"Could not read graph from file!\n";
      return EXIT_FAILURE;
   }
   solver.setProblem("test",graph.value());
   solver.solve();
   solver.printStatistics(std::cout);
   return EXIT_SUCCESS;
}