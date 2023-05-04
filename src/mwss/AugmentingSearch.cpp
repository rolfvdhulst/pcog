//
// Created by rolf on 26-4-23.
//

#include "pcog/mwss/AugmentingSearch.hpp"

#include <stack>
#include <iostream>
#include <numeric>

namespace pcog{
AugmentingSearch::AugmentingSearch(const std::vector<SafeWeight> &greedyWeights, const DenseGraph &graph) :
                                                                                                                graph{graph}{
   std::size_t num_nodes = graph.numNodes();
   assert(greedyWeights.size() == num_nodes);
   table.resize(num_nodes);
   for (std::size_t i = 0; i < num_nodes; ++i) {
      table[i].status = GreedySolutionStatus::FREE;
      table[i].tightness = 0;
      table[i].greedy_weight = greedyWeights[i];
   }
   num_solution_nodes = 0;
   num_free_nodes = num_nodes;
   temporary_set = SparseSet(graph.numNodes());
}

void AugmentingSearch::addNodeToSolution(Node node) {
   assert(table[node].status == GreedySolutionStatus::FREE);
   table[node].status = GreedySolutionStatus::SOLUTION;
   num_free_nodes--;
   num_solution_nodes++;

   for (const auto &neighbour : graph.neighbourhood(node)) {
      assert(table[neighbour].status != GreedySolutionStatus::SOLUTION);
      if (table[neighbour].status == GreedySolutionStatus::FREE) {
         assert(table[neighbour].tightness == 0);
         table[neighbour].status = GreedySolutionStatus::NONFREE;
         num_free_nodes--;
      }
      table[neighbour].tightness++;
   }

}

void AugmentingSearch::removeNodeFromSolution(Node node) {
   assert(table[node].status == GreedySolutionStatus::SOLUTION);
   table[node].status = GreedySolutionStatus::FREE;
   num_free_nodes++;
   num_solution_nodes--;

   for (const auto &neighbour : graph.neighbourhood(node)) {
      assert(table[neighbour].status == GreedySolutionStatus::NONFREE);
      assert(table[neighbour].tightness > 0);
      table[neighbour].tightness--;
      if (table[neighbour].tightness == 0) {
         table[neighbour].status = GreedySolutionStatus::FREE;
         num_free_nodes++;
      }
   }

}

AugmentingSearch::Weight AugmentingSearch::getWeight(Node node) const {
   return table[node].greedy_weight;
}

std::size_t AugmentingSearch::getTightness(Node node) const {
   return table[node].tightness;
}

bool AugmentingSearch::nodeInSolution(Node node) const {
   assert( node < table.size());
   return table[node].status == GreedySolutionStatus::SOLUTION;
}

bool AugmentingSearch::solutionIsMaximal() const {
   return num_free_nodes == 0;
}

bool AugmentingSearch::nodeFree(Node node) const {
   return table[node].status == GreedySolutionStatus::FREE;
}

bool AugmentingSearch::searchTwoImprovement(Node node) {
   //try to see if a two improvement exists.
   //Try to find a stable set among the neighbours of this node with tightness 1, essentially
   temporary_set.clear();
   for(const auto& neighbour : graph.neighbourhood(node)){
      if(getTightness(neighbour) == 1){
         temporary_set.unsafe_add(neighbour);
      }
   }
   if(temporary_set.size() < 2){
      return false;
   }
   Weight current_weight = getWeight(node);
   Node best_v = INVALID_NODE;
   Node best_w = INVALID_NODE;

   auto begin = temporary_set.cbegin();
   auto end = temporary_set.cend();
   auto end_minone = end;
   end_minone--;

   for (auto v_it = begin; v_it != end_minone; v_it++) {
      Node v = *v_it;
      Weight v_weight = getWeight(v);
      for (auto w_it = v_it + 1; w_it != end; w_it++) {
         Node w = *w_it;
         Weight w_weight = getWeight(w);
         if (v_weight + w_weight > current_weight && !graph.isEdge(v, w)) {
            best_v = v;
            best_w = w;
            current_weight = v_weight + w_weight;
         }
      }
   }
   if (best_v != INVALID_NODE) {
      removeNodeFromSolution(node);
      addNodeToSolution(best_v);
      addNodeToSolution(best_w);
   }
   return best_v != INVALID_NODE;
}

bool AugmentingSearch::doTwoImprovements() {
   bool found_nodes = true;
   bool changed = false;
   while(found_nodes){
      found_nodes = false;
      for(Node node = 0; node < graph.numNodes(); ++node){
         if(nodeInSolution(node)){
            bool result = searchTwoImprovement(node);
            found_nodes |= result;
            changed |= result;
         }
      }
   }
   return changed;
}

void AugmentingSearch::setSolution(const DenseSet &set) {
   clearSolution();
   for(const auto& node : set){
      addNodeToSolution(node);
   }
}

void AugmentingSearch::clearSolution() {
   for(Node node = 0; node < graph.numNodes(); ++node){
      table[node].tightness = 0;
      table[node].status = GreedySolutionStatus::FREE;
   }
   num_free_nodes = graph.numNodes();
   num_solution_nodes = 0;

}

std::size_t AugmentingSearch::currentSolutionSize() const {
   return num_solution_nodes;
}

bool AugmentingSearch::searchImprovements() {
   bool did_improve = false;

   while(true){
      bool improved = false;
      improved |= doTwoImprovements();
      //std::cout<<" weight after 2-improvents: "<< totalWeight()<<std::endl;
      improved |= doTwoKImprovements();
      //std::cout<<" weight after k-improvents: "<< totalWeight()<<std::endl;

      if(!improved){
         break;
      }else{
         did_improve = true;
      }
   }
   return did_improve;
}

DenseSet AugmentingSearch::getDenseSolution() const {
   DenseSet set(table.size());
   for (Node node = 0; node < table.size(); ++node) {
      if(nodeInSolution(node)){
         set.add(node);
      }
   }
   return set;
}

void AugmentingSearch::updateWeights(const std::vector<SafeWeight> &greedyWeights) {
   for (std::size_t i = 0; i < table.size(); ++i) {
      table[i].greedy_weight = greedyWeights[i];
   }
}

bool AugmentingSearch::doTwoKImprovements() {
   bool changed = false;
   bool changed_this_it = true;
   while(changed_this_it){
      changed_this_it = false;
      for(Node node = 0; node < graph.numNodes(); ++node){
         if(getTightness(node) == 2 && getWeight(node) > 0){
            bool result = searchTwoKImprovement(node);
            changed_this_it |= result;
            changed |= result;
         }
      }
   }

   return changed;
}
std::size_t AugmentingSearch::doTwoKImprovementsWithGreedy(){
   std::size_t changes = 0;
   std::vector<Node> nodes(table.size());
   std::iota(nodes.begin(), nodes.end(),0);
   std::vector<Weight> weights(table.size());
   for(std::size_t i = 0; i < table.size(); ++i){
      if(getTightness(i) != 2){
         weights[i] = 0;
      }else{
         weights[i] = getWeight(i);

      }
   }

   std::sort(nodes.begin(), nodes.end(),[&](const Node& a, const Node& b){
      return weights[a] > weights[b];
   });
   for(std::size_t i = 0; i < num_solution_nodes; ++i){
      Node node = nodes[i];
      if (!(getTightness(node) == 2 && getWeight(node) > 0)) {
         break;
      }
      bool result = searchTwoKImprovement(node);
      if (result) {
         changes++;
      }
   }

   return changes;
}


bool AugmentingSearch::searchTwoKImprovement(Node node) {
   assert(getTightness(node)== 2);
   std::stack<Node> nodeStack;
   nodeStack.push(node);
   std::vector<Node> add_list;
   std::vector<Node> remove_list;

   std::vector<int> is_solution_restricted(table.size());
   for (std::size_t i = 0; i < table.size(); ++i) {
      is_solution_restricted[i] = 0;
   }
   long total_weight = 0; //need signed weight here
   add_list.push_back(node);
   total_weight+= getWeight(node);
   for(const Node& marked_node : graph.neighbourhood(node)){
      is_solution_restricted[marked_node]+=1;
   }

   while(!nodeStack.empty()){
      Node x = nodeStack.top();
      nodeStack.pop();
      if(nodeInSolution(x)){
         for(const Node& neighbour : graph.neighbourhood(x)){
            if(getTightness(neighbour) == 1 && getWeight(neighbour) > 0){
               if(!is_solution_restricted[neighbour]){
                  add_list.push_back(neighbour);
                  total_weight+= getWeight(neighbour);
                  for(const Node& mark_node: graph.neighbourhood(neighbour)){
                     is_solution_restricted[mark_node]+=1;
                  }
               }
            }
         }
      }else{
         for(const Node& neighbour : graph.neighbourhood(x)){
            if(nodeInSolution(neighbour)){
               remove_list.push_back(neighbour);
               nodeStack.push(neighbour);
               total_weight-= getWeight(neighbour);
            }
         }
      }
   }
   if(total_weight > 0){
      for(const Node& remove_node : remove_list){
         removeNodeFromSolution(remove_node);
      }
      for(const Node& add_node : add_list){
         addNodeToSolution(add_node);
      }
   }
   return total_weight > 0;
}

AugmentingSearch::Weight AugmentingSearch::totalWeight() const {
   Weight sum = 0;
   for(Node node = 0; node < table.size();++node ){
      if(nodeInSolution(node)){
         sum+=getWeight(node);
      }
   }
   return sum;
}

std::vector<Node> AugmentingSearch::decreasingWeightOrdering() const {
   std::vector<Node> nodes(table.size());
   std::iota(nodes.begin(),nodes.end(),0);
   std::sort(nodes.begin(),nodes.end(),[&](const Node& a, const Node& b){
      return getWeight(a) > getWeight(b);
   });
   return nodes;
}

std::size_t AugmentingSearch::greedyByWeight(long rotateOrderBy) {
   auto ordering = decreasingWeightOrdering();
   std::rotate(ordering.begin(), ordering.begin() + rotateOrderBy, ordering.end());
   std::size_t nodesAdded = 0;
   for (const auto &node: ordering) {
      if (nodeFree(node)) {
         addNodeToSolution(node);
         nodesAdded++;
      }
   }
   return nodesAdded;
}




std::size_t AugmentingSearch::greedyByDynamicSurplus(long rotateOrderBy) {

   //O(n^2), slow. Heap usage could lead to, but this is a pain to implement
   //as one needs to continuously update it when picking vertices

   return 0;
}

std::size_t AugmentingSearch::greedyByStaticSurplus(long rotateOrderBy) {
   std::vector<SafeWeight> neighbourhood_sum(table.size(),0);
   for(Node i = 0; i < table.size(); ++i){
      for(const auto& neighbour : graph.neighbourhood(i)){
         neighbourhood_sum[i] += table[neighbour].greedy_weight;
      }
   }
   std::vector<Node> nodes(table.size());
   std::iota(nodes.begin(),nodes.end(),0);
   std::sort(nodes.begin(),nodes.end(),[&](const Node& a, const Node& b){
      auto nodeWeightA = getWeight(a);
      auto nodeWeightB = getWeight(b);
      bool nodeAPositive = nodeWeightA >= neighbourhood_sum[a];
      bool nodeBPositive = nodeWeightB >= neighbourhood_sum[b];
      if(nodeAPositive){
         if(nodeBPositive){
            return (nodeWeightA-neighbourhood_sum[a]) > (nodeWeightB-neighbourhood_sum[b]);
         }
         return true; //Order is correcct
      }
      if(!nodeBPositive){
         return (neighbourhood_sum[a]-nodeWeightA) < (neighbourhood_sum[b]-nodeWeightB);
      }
      return false;
   });
   std::rotate(nodes.begin(), nodes.begin() + rotateOrderBy, nodes.end());
   std::size_t nodesAdded = 0;
   for (const auto &node: nodes) {
      if (nodeFree(node)) {
         addNodeToSolution(node);
         nodesAdded++;
      }
   }
   return nodesAdded;
}

}