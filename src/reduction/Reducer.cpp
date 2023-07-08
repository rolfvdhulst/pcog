//
// Created by rolf on 6-7-23.
//

#include "pcog/reduction/Reducer.hpp"

#include <iostream> //TODO: remove after testing
#include <chrono>
namespace pcog{

void reduceGraph(const DenseGraph& graph){
   auto start = std::chrono::high_resolution_clock::now();
   std::cout<<"Original graph has: "<<graph.numNodes()<<" nodes\n";
   DenseReductionGraph redGraph(graph);
   ReductionStack result;
   ReductionVertexQueue queue(graph.numNodes(),true);
   while(!queue.empty()){
      Node node = queue.pop();
      if(!redGraph.containsNode(node)){
         continue;
      }
      if(dominatedReduceNodeDense(node,redGraph,result,queue)){
         continue;
      }
      else if(simplicialReduceNode(node,redGraph,result,queue)){
         continue;
      }
   }

   std::cout<<"Preprocessed graph has: "<<redGraph.nodes().size()<<" nodes\n";
   auto end = std::chrono::high_resolution_clock::now();
   std::cout<<"Took: "<< std::chrono::duration_cast<std::chrono::milliseconds>(end-start)<<"\n";
   exit(1);
}

}