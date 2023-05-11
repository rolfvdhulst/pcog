//
// Created by rolf on 6-5-23.
//

#include "pcog/TabuColoring.hpp"
namespace pcog {

std::optional<NodeColoring> TabuColoring::run(NodeColoring initialColoring) {

   NodeColoring coloring = initialColoring;
   std::size_t num_nodes = graph.numNodes();
   std::vector<std::vector<std::size_t>> tabu_list(
       num_nodes, std::vector<std::size_t>(coloring.numColors(), 0));
   std::vector<std::vector<std::size_t>> adjacency_list(
       num_nodes, std::vector<std::size_t>(coloring.numColors(), 0));
   long objective = 0;
   for (Node node = 0; node < num_nodes; ++node) {
      std::size_t color = coloring[node];
      const DenseSet &neighbourhood = graph.neighbourhood(node);
      for (const Node &other_node : neighbourhood) {
         std::size_t other_color = coloring[other_node];
         (adjacency_list[node][other_color])++;
         if (color == other_color) {
            objective++;
         }
      }
   }
   assert(objective % 2 == 0);
   objective /= 2;
   //assert(objective == numViolatedEdges(graph, coloring)); //TODO: commented out since it really slows down debug performance

   long bestObjective = objective;
   done_iterations = max_iterations;
   for (std::size_t iter = 0; iter < max_iterations; iter++) {

      Node bestNode = INVALID_NODE;
      std::size_t best_color = -1;
      std::size_t num_critical = 0;
      long min_value = graph.numNodes() * graph.numNodes();

      for (Node node = 0; node < graph.numNodes(); ++node) {
         bool aspiration = false;
         std::size_t node_color = coloring[node];
         assert(node_color < coloring.numColors());

         if (adjacency_list[node][node_color] > 0) {
            num_critical++;
            for (std::size_t color = 0; color < coloring.numColors(); color++) {
               if (color != node_color) {
                  /* change in # of violated edges */
                  long change = adjacency_list[node][color] -
                                adjacency_list[node][node_color];
                  if ((objective + change) == 0) {
                     bestNode = node;
                     best_color = color;
                     min_value = change;
                     aspiration = true;
                  }
                  if (tabu_list[node][color] < iter && change < min_value) {
                     bestNode = node;
                     best_color = color;
                     min_value = change;
                  }
               }
            }
         }
         if (aspiration) {
            break;
         }
      }
      // current tabu list is too restrictive/ e.g. no good choices exist; just skip to next iteration
      if (bestNode == INVALID_NODE) {
         continue;
      }
      assert(bestNode != INVALID_NODE);
      assert(best_color < coloring.numColors());
      assert(coloring[bestNode] != best_color);
      std::size_t oldColor = coloring[bestNode];
      coloring[bestNode] = best_color;
      objective += min_value;
      //assert(objective == numViolatedEdges(graph, coloring));//TODO: commented out since it really slows down debug performance
      if (objective < bestObjective) {
         bestObjective = objective;
      }
      assert(objective >= 0);
      if (objective == 0) {
         done_iterations = iter;
         break;
      }
      /*update tabu list*/
      std::size_t forbidden_iterations = tabu_base + num_critical * tabu_gamma;
      tabu_list[bestNode][oldColor] = iter + forbidden_iterations;

      for (const Node &node : graph.neighbourhood(bestNode)) {
         (adjacency_list[node][best_color])++;
         (adjacency_list[node][oldColor])--;
      }
   }

   if (objective == 0) {
      return coloring;
   } else {
      return std::nullopt;
   }
}

long TabuColoring::numViolatedEdges(const DenseGraph &graph,
                                  const NodeColoring &coloring) {
   long numViolated = 0;
   for (Node node = 0; node < graph.numNodes(); ++node) {
      for (const Node &other_node : graph.neighbourhood(node)) {
         if (other_node >= node) {
            break;
         }
         if (coloring[node] == coloring[other_node]) {
            numViolated++;
         }
      }
   }
   return numViolated;
}

TabuColoring::TabuColoring(const DenseGraph &graph) : graph{graph} {}

std::size_t TabuColoring::numIterations() const { return done_iterations; }
void TabuColoring::setMaxIterations(std::size_t iters) {
   max_iterations = iters;
}
void TabuColoring::setTabuBase(std::size_t numBase) {
   tabu_base = numBase;
}
void TabuColoring::setGamma(double gamma) {
   tabu_gamma = gamma;
}
}
