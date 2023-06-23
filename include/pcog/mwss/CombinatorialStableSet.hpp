//
// Created by rolf on 27-2-23.
//

#ifndef PCOG_INCLUDE_PCOG_MWSS_COMBINATORIALSTABLESET_HPP
#define PCOG_INCLUDE_PCOG_MWSS_COMBINATORIALSTABLESET_HPP

// dimensions;
// weight function; unicost? linear (node weights)? or more complicated (cuts?)
// safe weights or not? -> function responsibility

// incremental computation of total set weight? For now total set weight,
// possibly support incremental in the future. Pruning strategies -> implicitly
// excluded nodes? Partition strategies

#include "WeightFunction.hpp"
#include "pcog/utilities/DenseGraph.hpp"
#include "pcog/utilities/DenseSet.hpp"
#include "pcog/utilities/SparseSet.hpp"
#include <chrono>
#include <utility>

#include <iostream>
namespace pcog {
template <WeightFunction F> class MaxWeightStableSetCombinatorial {
 public:
   explicit MaxWeightStableSetCombinatorial(F function,
                                            const DenseGraph &graph);
   using weight_type = typename F::weight_type;

   void run();
   typedef void SolutionCallback(const DenseSet &current_nodes,
                                 weight_type total_weight, void *user_data,
                                 bool first_solution, bool &stop_solving,
                                 bool &accepted_solution);

   const F &getWeightFunction() const;
   void updateWeightFunction(F function);
   void setLowerBound(weight_type lb);
   void setUpperBound(weight_type ub);
   void setInfiniteUpperBound();

   template<typename T>
   void setUserData(T *userData);
   void setCallback(SolutionCallback callback);
   void setNodeLimit(std::size_t num_nodes);
   void setTimeLimit(std::chrono::duration<double> seconds);
   void setTimeLimitCheckInterval(std::size_t interval_nodes);
   void setInterrupt(volatile bool * interrupt);

   [[nodiscard]] bool stoppedBeforeOptimality() const;
   [[nodiscard]] std::size_t numBranchAndBoundNodes() const;
   [[nodiscard]] std::chrono::high_resolution_clock::duration timeTaken() const;

 private:
   bool branch(std::size_t depth);
   struct CliqueCoveringResult {
      bool found_full_covering;
      weight_type weight;
   };
   CliqueCoveringResult
   clique_covering_heuristic_held(const DenseSet &free_nodes,
                                  SparseSet &branch_nodes, weight_type bound);
   bool executeCallback(bool improves, const DenseSet &set, weight_type weight);
   bool isPruned(const DenseSet &excluded_nodes, const DenseSet &free_nodes,
                 const DenseSet &currentNodes);
   void selectBranchingNodes(SparseSet &branch_nodes,
                             const DenseSet &excluded_nodes,
                             const DenseSet &free_nodes,
                             const DenseSet &currentNodes);
   struct BranchData {
      explicit BranchData(std::size_t capacity);
      DenseSet currentSet;
      DenseSet freeNodes;
      DenseSet excludedNodes;
      SparseSet branchNodes;
   };

   const DenseGraph &graph;

   F weightFunction;
   std::vector<weight_type> remaining_weights;

   std::vector<BranchData> memory_stack;

   weight_type global_lower_bound = 0;
   bool has_ub = false;
   weight_type global_upper_bound = 0;

   // user settings
   void *user_data = nullptr;
   SolutionCallback *callback = nullptr;
   volatile bool * interrupt = nullptr;

   std::chrono::duration<double> time_limit = std::chrono::duration<double>::max();
   std::size_t time_limit_check_interval = std::numeric_limits<std::size_t>::max();
   std::size_t node_limit = std::numeric_limits<std::size_t>::max();
   std::chrono::high_resolution_clock::time_point start_time;

   // statistics
   std::size_t num_branch_and_bound_nodes = 0;
   bool stopped_before_optimality = true;
   bool found_first_solution = false;
   std::chrono::high_resolution_clock::duration time_taken =
       std::chrono::high_resolution_clock::duration(0);
};
template <WeightFunction F>
void MaxWeightStableSetCombinatorial<F>::setInterrupt(
    volatile bool *t_interrupt) {
   interrupt = t_interrupt;
}

template <WeightFunction F>
MaxWeightStableSetCombinatorial<F>::BranchData::BranchData(std::size_t capacity)
    : currentSet(capacity), freeNodes(capacity), excludedNodes(capacity),
      branchNodes(capacity) {}

template <WeightFunction F>
void MaxWeightStableSetCombinatorial<F>::setLowerBound(weight_type lb) {
   global_lower_bound = lb;
}

template <WeightFunction F>
void MaxWeightStableSetCombinatorial<F>::setUpperBound(weight_type ub) {
   global_upper_bound = ub;
   has_ub = true;
}

template <WeightFunction F>
void MaxWeightStableSetCombinatorial<F>::setInfiniteUpperBound() {
   has_ub = false;
   global_upper_bound = std::numeric_limits<weight_type>::max();
}

template <WeightFunction F>
void MaxWeightStableSetCombinatorial<F>::setTimeLimit(
   std::chrono::duration<double> seconds) {
   time_limit = seconds;
}

template <WeightFunction F>
void MaxWeightStableSetCombinatorial<F>::setTimeLimitCheckInterval(
    std::size_t interval_nodes) {
   time_limit_check_interval = interval_nodes;
}

template <WeightFunction F>
bool MaxWeightStableSetCombinatorial<F>::stoppedBeforeOptimality() const {
   return stopped_before_optimality;
}

template <WeightFunction F>
std::size_t MaxWeightStableSetCombinatorial<F>::numBranchAndBoundNodes() const {
   return num_branch_and_bound_nodes;
}

template <WeightFunction F>
std::chrono::high_resolution_clock::duration
MaxWeightStableSetCombinatorial<F>::timeTaken() const {
   return time_taken;
}

template <WeightFunction F>
bool MaxWeightStableSetCombinatorial<F>::branch(std::size_t depth) {
   BranchData &current_data = memory_stack[depth];
   const DenseSet &current_set = current_data.currentSet;
   const DenseSet &free_nodes = current_data.freeNodes;
   const DenseSet &excluded_nodes = current_data.excludedNodes;
   SparseSet &branch_nodes = current_data.branchNodes;
   branch_nodes.clear();
   num_branch_and_bound_nodes++;

   weight_type current_weight = weightFunction.setCost(current_set);
   if (has_ub && current_weight > global_upper_bound) {
      return false; // we know through a global upper bound that the current
                    // branch cannot find any new stable set
   }

   bool improves = current_weight > global_lower_bound;
   if (free_nodes.empty()) {
      return executeCallback(improves, current_set, current_weight);
   }
   // check nodelimit and timelimit
   if (num_branch_and_bound_nodes >= node_limit) {
      stopped_before_optimality = true;
      return true;
   }
   if (num_branch_and_bound_nodes % time_limit_check_interval == 0) {
      auto end = std::chrono::high_resolution_clock::now();
      if ((interrupt != nullptr && *interrupt) || (end - start_time) > time_limit) {
         stopped_before_optimality = true;
         return true;
      }
   }

   weight_type bound = global_lower_bound - current_weight;
   if (current_weight >= global_lower_bound) { // preventing underflows
      bound = 0;
   }
   // partition the free nodes into two sets; one set is a set S such that W(S)
   // <= bound for every stable set in S The other nodes can then be used as
   // branching nodes. The first set is computed by using a clique-covering
   // argument, involving duality (there are more ways to do this) Clearly, if
   // we can conclude that for every stable set in F (the free nodes) W(F) <=
   // bound, we do not need to search this branch any further.

   CliqueCoveringResult result =
       clique_covering_heuristic_held(free_nodes, branch_nodes, bound);
   if (result.found_full_covering) {
      return false;
   }
   // check pruning rules
   if (isPruned(excluded_nodes, free_nodes, current_set)) {
      return false;
   }

   // select branching nodes based on results found in previous steps
   selectBranchingNodes(branch_nodes, excluded_nodes, free_nodes, current_set);

   auto max_weight_lambda = [&](const Node &a, const Node &b) -> bool {
      return remaining_weights[a] >
             remaining_weights[b]; // TODO: perhaps use degrees when tied
   };
   branch_nodes.sortBy(max_weight_lambda);

   // loop over branch nodes and recursively solve subproblems
   depth++;
   BranchData &next_data = memory_stack[depth];
   next_data.excludedNodes = excluded_nodes; // TODO: prevent reallocation
   for (const auto &branch_node : branch_nodes) {
      next_data.freeNodes = free_nodes.difference(
          graph.neighbourhood(branch_node)); // TODO: prevent reallocation
      next_data.freeNodes.inplaceDifference(
          next_data.excludedNodes); // implicitly also removes all previous
                                    // branch nodes
      // excluded nodes contain all previous branch nodes but the node chosen
      // now
      next_data.freeNodes.remove(branch_node);

      // remove neighbourhood(branch_node) and any nodes which were also
      // excluded update our current set; this could be made more efficient
      // probably
      next_data.currentSet = current_set; // TODO: prevent reallocation
      next_data.currentSet.add(branch_node);

      bool stop_solving = branch(depth);
      if (stop_solving) {
         return true;
      }
      next_data.excludedNodes.add(branch_node);
   };
   return false;
}

template <WeightFunction F>
bool MaxWeightStableSetCombinatorial<F>::isPruned(
    const DenseSet &excluded_nodes, const DenseSet &free_nodes,
    const DenseSet &currentNodes) {
   DenseSet temporary_nodes =
       free_nodes.setUnion(currentNodes); // TODO: prevent reallocations
   for (const auto &node : excluded_nodes) {
      const DenseSet &excluded_neighbourhood = graph.neighbourhood(node);
      auto intersection = temporary_nodes.intersection(
          excluded_neighbourhood); // TODO: prevent reallocations
      if (weightFunction.nodeCost(node) >=
          weightFunction.setCost(intersection)) {
         return true;
      }
   }
   return false;
}

template <WeightFunction F>
void MaxWeightStableSetCombinatorial<F>::selectBranchingNodes(
    SparseSet &branch_nodes, const DenseSet &excluded_nodes,
    const DenseSet &free_nodes, const DenseSet &currentNodes) {

   assert(!branch_nodes.empty());
   // Look for nodes which have a weight >= than their entire free
   // neighbourhood; these are always better to pick than their neighbourhood in
   // that case
   Node bestFreeNode = INVALID_NODE;
   weight_type bestWeight = 0;
   for (const auto &node : free_nodes) {
      const DenseSet &neighbourhood = graph.neighbourhood(node);
      auto temporarySet =
          neighbourhood.intersection(free_nodes); // TODO: prevent reallocations
      weight_type node_weight = weightFunction.nodeCost(node);
      weight_type set_weight = weightFunction.setCost(temporarySet);
      if (node_weight >= set_weight && node_weight - set_weight >= bestWeight) {
         bestFreeNode = node;
         bestWeight = node_weight - set_weight;
      }
   }
   if (bestFreeNode != INVALID_NODE) {
      branch_nodes.clear();
      branch_nodes.unsafe_add(bestFreeNode);
      return;
   }

   Node bestExcludedNode = INVALID_NODE;
   degree_type best_size =
       branch_nodes.size(); // TODO: by size, or by violation?
   for (const auto &node : excluded_nodes) {
      const DenseSet &neighbourhood = graph.neighbourhood(node);
      auto temporary1 =
          neighbourhood.intersection(free_nodes); // prevent reallocations TODO
      degree_type size = temporary1.size();
      if (size != 0 && size < best_size) {
         auto temporary2 = neighbourhood.intersection(
             excluded_nodes); // prevent reallocations TODO
         weight_type node_weight = weightFunction.nodeCost(node);
         weight_type set_weight = weightFunction.setCost(temporary2);
         if (node_weight >= set_weight) {
            best_size = size;
            bestExcludedNode = node;
         }
      }
   }
   if (bestExcludedNode != INVALID_NODE) {
      auto branchSetDense =
          graph.neighbourhood(bestExcludedNode)
              .intersection(free_nodes); // TODO: prevent reallcation
      branch_nodes.fromDenseSet(branchSetDense);
      // TODO sort according to some metric
   }
}

template <WeightFunction F>
void MaxWeightStableSetCombinatorial<F>::setNodeLimit(std::size_t num_nodes) {
   node_limit = num_nodes;
}

template <WeightFunction F>
bool MaxWeightStableSetCombinatorial<F>::executeCallback(
    bool weight_improves, const DenseSet &set, weight_type current_weight) {
   bool stop_solving = false;
   if (weight_improves) {
      bool accept_solution = false;
      callback(set, current_weight, user_data, !found_first_solution,
               stop_solving, accept_solution);
      if (accept_solution) {
         global_lower_bound = std::max(global_lower_bound, current_weight);
         // we only update the LB if the user accepts this solution; otherwise
         // we might run into issues later
         found_first_solution = true;
      }
   }
   if (num_branch_and_bound_nodes >= node_limit) {
      stop_solving = true; // override no stop signal in case user finds a
                           // solution at the exact stopping depth
   }
   if (stop_solving) {
      stopped_before_optimality = true;
   }
   return stop_solving;
}

template <WeightFunction F>
typename MaxWeightStableSetCombinatorial<F>::CliqueCoveringResult
MaxWeightStableSetCombinatorial<F>::clique_covering_heuristic_held(
    const DenseSet &free_nodes, SparseSet &branch_nodes, weight_type bound) {
   // TODO: figure out how to handle nodes weight zero in this covering

   weight_type total_weight = weightFunction.fixedDualCost();
   for (const Node &node : free_nodes) {
      remaining_weights[node] = weightFunction.dualNodeCost(node);
      // TODO: can be precomputed in an array and just copied maybe? Could safe
      // some time for complicated functions
   }

   if (total_weight > bound) {
      branch_nodes.fromDenseSet(free_nodes);
      // TODO sort branch_nodes by some metric
      return CliqueCoveringResult{.found_full_covering = false,
                                  .weight = total_weight};
   }
   std::size_t num_covered_nodes = 0;
   std::size_t num_free_nodes = free_nodes.size();
   DenseSet uncovered_nodes = free_nodes; // TODO: prevent reallocation
   while (num_covered_nodes != num_free_nodes) {
      // incrementally construct a clique covering
      Node best_node = INVALID_NODE;
      weight_type smallest_weight = std::numeric_limits<weight_type>::max();
      for (const Node &node : uncovered_nodes) {
         weight_type node_weight = remaining_weights[node];
         if (node_weight < smallest_weight) {
            smallest_weight = node_weight;
            best_node = node;
         }
      }
      if (best_node == INVALID_NODE) {
         assert(false); //TODO: sometimes hit, how?
         break; // all nodes should have been covered
      }
      assert(total_weight + smallest_weight >= total_weight); // overflow check

      if (total_weight + smallest_weight > bound) {
         branch_nodes.fromDenseSet(uncovered_nodes);
         // TODO sort branch_nodes by some metric?
         return CliqueCoveringResult{.found_full_covering = false,
                                     .weight = total_weight};
      }

      // greedily build clique around best_node
      // TODO: might be more efficient with a sparse set and a single sort
      // beforehand
      DenseSet clique_nodes(graph.numNodes()); // TODO: prevent reallocation
      clique_nodes.add(best_node);

      DenseSet searchable_nodes = uncovered_nodes.intersection(
          graph.neighbourhood(best_node)); // TODO: prevent reallocation

      while (true) {
         Node expand_clique_node = INVALID_NODE;
         weight_type smallest_cliquenode_weight =
             std::numeric_limits<weight_type>::max();
         for (const Node &node : searchable_nodes) {
            weight_type node_weight = remaining_weights[node];
            if (node_weight < smallest_cliquenode_weight) {
               smallest_cliquenode_weight = node_weight;
               expand_clique_node = node;
            }
         }
         if (expand_clique_node == INVALID_NODE) {
            break;
         }
         searchable_nodes.inplaceIntersection(
             graph.neighbourhood(expand_clique_node));
         clique_nodes.add(expand_clique_node);
      }

      for (const auto &node : clique_nodes) {
         weight_type &node_weight = remaining_weights[node];
         if (smallest_weight >= node_weight) {
            node_weight = 0;
            uncovered_nodes.remove(node);
            num_covered_nodes++;
         } else {
            node_weight -= smallest_weight;
         }
      }
      total_weight += smallest_weight;
   }
   return CliqueCoveringResult{.found_full_covering = true,
                               .weight = total_weight};
}

template <WeightFunction F>
MaxWeightStableSetCombinatorial<F>::MaxWeightStableSetCombinatorial(
    F function, const DenseGraph &graph)
    : graph{graph}, weightFunction{function} {
   memory_stack.reserve(graph.numNodes());
   for (std::size_t i = 0; i < graph.numNodes(); ++i) {
      memory_stack.emplace_back(BranchData(graph.numNodes()));
   }
   remaining_weights.resize(graph.numNodes());
   for (std::size_t i = 0; i < graph.numNodes(); ++i) {
      remaining_weights[i] = 0;
   }
}

template <WeightFunction F> void MaxWeightStableSetCombinatorial<F>::run() {
   start_time = std::chrono::high_resolution_clock::now();

   stopped_before_optimality = false;
   found_first_solution = false;
   num_branch_and_bound_nodes = 0;

   BranchData &first_call = memory_stack[0];
   first_call.freeNodes.setAll();
   // TODO: add option for zero-weight-node extension?
   // TODO: add option of initializing currentset with a certain set,
   // removing all neighbourhoods of that set and itself from the free set of
   // nodes

   for (Node i = 0; i < first_call.freeNodes.capacity(); ++i) {
      if (weightFunction.nodeCost(i) ==
          0) { // TODO: does not work for nonlinear; need to check if
               // setCost(node) == 0
         first_call.freeNodes.remove(i);
      }
   }
   first_call.excludedNodes.clear();
   first_call.currentSet.clear();
   first_call.branchNodes.clear();

   branch(0);
   auto end = std::chrono::high_resolution_clock::now();
   time_taken = (end - start_time);
}

template <WeightFunction F>
template<typename T>
void MaxWeightStableSetCombinatorial<F>::setUserData(T *userData) {
   user_data = reinterpret_cast<void *>(userData);
}

template <WeightFunction F>
void MaxWeightStableSetCombinatorial<F>::setCallback(
    MaxWeightStableSetCombinatorial<F>::SolutionCallback solutionCallback) {
   callback = solutionCallback;
}

template <WeightFunction F>
void MaxWeightStableSetCombinatorial<F>::updateWeightFunction(F function) {
   weightFunction = function;
}

template <WeightFunction F>
const F &MaxWeightStableSetCombinatorial<F>::getWeightFunction() const {
   return weightFunction;
}

struct UniformWeightFunction {
   using weight_type = std::size_t;
   static constexpr WeightFunctionType function_type =
       WeightFunctionType::UNIFORM_COST;
   static weight_type nodeCost(Node ) { return 1; };
   static weight_type setCost(const DenseSet &set) { return set.size(); };
   static weight_type dualNodeCost(Node ) { return 1; };
   static weight_type fixedDualCost() { return 0; };
};

struct PartialUniformWeightFunction {
   using weight_type = std::size_t;
   PartialUniformWeightFunction(DenseSet set) : mask{std::move(set)} {}
   static constexpr WeightFunctionType function_type =
       WeightFunctionType::LINEAR;
   weight_type nodeCost(Node node) const { return mask.contains(node); }
   weight_type setCost(const DenseSet &set) const {
      return set.intersection(mask).size();
   }
   weight_type dualNodeCost(Node node) const { return nodeCost(node); }
   static weight_type fixedDualCost() { return 0; }

 private:
   DenseSet mask;
};

class LinearWeightFunction {
 public:
   using weight_type = std::size_t;
   LinearWeightFunction(std::vector<std::size_t> nodes)
       : nodeWeights(std::move(nodes)){};
   static constexpr WeightFunctionType function_type =
       WeightFunctionType::LINEAR;
   weight_type nodeCost(Node node) const { return nodeWeights[node]; };
   weight_type setCost(const DenseSet &set) const {
      weight_type sum = 0;
      for (const auto &node : set) {
         sum += nodeCost(node);
      }
      return sum;
   };
   weight_type dualNodeCost(Node node) const { return nodeCost(node); };
   static weight_type fixedDualCost() { return 0; };

 private:
   std::vector<std::size_t> nodeWeights;
};
using MaxStableSetCombinatorial =
    MaxWeightStableSetCombinatorial<UniformWeightFunction>;
using MaxLinearWeightStableSetCombinatorial =
    MaxWeightStableSetCombinatorial<LinearWeightFunction>;
} // namespace pcog
#endif // PCOG_INCLUDE_PCOG_MWSS_COMBINATORIALSTABLESET_HPP
