//
// Created by rolf on 6-7-23.
//

#include "pcog/reduction/Reducer.hpp"

#include <iostream> //TODO: remove after testing
#include <chrono>
namespace pcog{
DenseSet greedyClique(const DenseReductionGraph& t_graph){
   assert(!t_graph.nodes().empty());
   Node bestNode = t_graph.nodes().first();
   std::size_t bestDegree = 0;
   Node iterNode = bestNode;
   assert(iterNode != INVALID_NODE);
   while(iterNode != INVALID_NODE){
      std::size_t degree = t_graph.nodeDegree(iterNode);
      if(degree > bestDegree){
         bestNode = iterNode;
         bestDegree = degree;
      }
      iterNode = t_graph.nodes().find_next(iterNode);
   }

   DenseSet clique(t_graph.nodes().capacity());
   clique.add(bestNode);
   DenseSet candidates = t_graph.neighbourhood(bestNode);
   do{
      bestNode = INVALID_NODE;
      bestDegree = 0;
      for(const auto& candidate : candidates){
         std::size_t degree = t_graph.neighbourhood(candidate).intersection(candidates).size();
         if(degree > bestDegree){
            bestDegree = degree;
            bestNode = candidate;
         }
      }
      if(bestNode == INVALID_NODE){
         break;
      }
      assert(bestNode != INVALID_NODE);

      clique.add(bestNode);
      candidates.inplaceIntersection(t_graph.neighbourhood(bestNode));
   }while(!candidates.empty());

   return clique;
}

void reduceGraph(const DenseGraph& graph){
   auto start = std::chrono::high_resolution_clock::now();
   std::cout<<"Original graph has: "<<graph.numNodes()<<" nodes\n";
   DenseReductionGraph redGraph(graph);
   DenseSet clique = greedyClique(redGraph); //TODO: add more advanced way to find clique
   redGraph.setLowerBoundClique(clique,clique.size());
   ReductionStack result;
   ReductionVertexQueue queue(clique);
   while(!queue.empty()){
      Node node = queue.pop();
      if(!redGraph.containsNode(node)){
         continue;
      }
      if(simplicialReduceNode(node,redGraph,result,queue)){
         std::cout<<"Fixed set!\n";
         continue;
      }
      else if(lowDegreeReduceNode(node,redGraph,result,queue)){
         std::cout<<"Low degree!\n";
         continue;
      }
      else if(foldDegreeTwoReduceNode(node,redGraph,result,queue)){
         std::cout<<"Degree two reduced!\n";
         continue;
      }
      else if(twinDegreeThreeReduction(node,redGraph,result,queue)){
         std::cout<<"Twin degree 3 reduced!\n";
         continue;
      }
      else if(dominatedReduceNode(node,redGraph,result,queue)){
         std::cout<<"Dominated!\n";
         continue;
      }


   }

   std::cout<<"Preprocessed graph has: "<<redGraph.nodes().size()<<" nodes\n";
   auto end = std::chrono::high_resolution_clock::now();
   std::cout<<"Took: "<< std::chrono::duration_cast<std::chrono::milliseconds>(end-start)<<"\n";
}

}