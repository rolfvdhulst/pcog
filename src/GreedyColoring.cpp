//
// Created by rolf on 6-5-23.
//


#include "pcog/GreedyColoring.hpp"
namespace pcog{
SetColoring GreedyColoring::run_degree_node() const {
   SetColoring coloring;
   std::vector<std::pair<std::size_t,std::size_t>> degrees(graph.numNodes());
   for (std::size_t i = 0; i < graph.numNodes() ; ++i) {
      degrees[i] = std::make_pair(i,graph.neighbourhood(i).size());
   }
   std::sort(degrees.begin(),degrees.end(),[](const std::pair<std::size_t,std::size_t>& a,
                                                const std::pair<std::size_t,std::size_t>& b){
      return a.second > b.second;
   });

   while(!degrees.empty()){
      std::size_t considered_node = degrees.back().first;
      const auto& neighbourhood = graph.neighbourhood(considered_node);
      bool none_found = true;
      for(auto& color : coloring.colors()){
         DenseSet intersection = neighbourhood.intersection(color);
         if(intersection.empty()){
            color.add(considered_node);
            none_found = false;
            break;
         }
      }
      if(none_found){
         DenseSet color(graph.numNodes());
         color.add(considered_node);
         coloring.addColor(color);
      }
      degrees.pop_back();
   }
#ifndef NDEBUG
   for(const auto& color : coloring.colors()){
      assert(graph.setIsStable(color));
   }
   DenseSet complete_covered(graph.numNodes());
   for(const auto& color : coloring.colors()){
      complete_covered.inplaceUnion(color);
   }
   assert(complete_covered.full());
#endif
   return coloring;
}
GreedyColoring::GreedyColoring(const DenseGraph &graph) : graph{graph}{

}

SetColoring GreedyColoring::run_saturation_degree() const {
   SetColoring coloring;
   std::vector<std::size_t> saturation(graph.numNodes(),0);
   std::vector<std::size_t> degrees(graph.numNodes());
   for (std::size_t i = 0; i < graph.numNodes() ; ++i) {
      degrees[i] = graph.neighbourhood(i).size();
   }
   auto max_degree = std::max_element(degrees.begin(),degrees.end());
   auto max_degree_node = std::distance(degrees.begin(),max_degree);
   DenseSet first_color(graph.numNodes());
   first_color.add(max_degree_node);
   coloring.addColor(first_color);
   DenseSet uncolored_nodes(graph.numNodes(),true);
   uncolored_nodes.remove(max_degree_node);

   std::vector<DenseSet> color_saturated_nodes;
   color_saturated_nodes.push_back(graph.neighbourhood(max_degree_node));
   std::size_t saturation_node = color_saturated_nodes.front().first();
   while(saturation_node != DenseSet::INVALID_ELEMENT){
      saturation[saturation_node]++;
      saturation_node = color_saturated_nodes.front().find_next(saturation_node);
   }

   std::size_t nodes = graph.numNodes() - 1;
   while(nodes != 0){
      //pick most saturated node (degree is tiebreaker)
      std::size_t uncolored_node = uncolored_nodes.first();
      std::size_t best_saturation = saturation[uncolored_node];
      std::size_t best_degree = degrees[uncolored_node];
      std::size_t best_node = uncolored_node;

      while(uncolored_node != DenseSet::INVALID_ELEMENT){
         if(saturation[uncolored_node] > best_saturation ||
             (saturation[uncolored_node] == best_saturation && degrees[uncolored_node] > best_degree)){
            best_node = uncolored_node;
            best_saturation = saturation[uncolored_node];
            best_degree = degrees[uncolored_node];
         }
         uncolored_node = uncolored_nodes.find_next(uncolored_node);
      }

      //color node using the first available color
      bool colored = false;
      auto& colors = coloring.colors();
      for(std::size_t i = 0; i < colors.size(); i++){
         if(!color_saturated_nodes[i].contains(best_node)){
            colors[i].add(best_node);
            colored = true;
            //update saturation of colored nodes
            DenseSet temporary = color_saturated_nodes[i];
            color_saturated_nodes[i].inplaceUnion(graph.neighbourhood(best_node));
            DenseSet iterate_set =  color_saturated_nodes[i];
            iterate_set.inplaceDifference(temporary);
            std::size_t node = iterate_set.first();
            while(node != DenseSet::INVALID_ELEMENT){
               saturation[node]++;
               node = iterate_set.find_next(node);
            }
            break;
         }
      }
      if(!colored){
         DenseSet newcolor(graph.numNodes());
         newcolor.add(best_node);
         coloring.addColor(newcolor);
         color_saturated_nodes.push_back(graph.neighbourhood(best_node));
         //update saturation of newly colored nodes
         std::size_t node = color_saturated_nodes.back().first();
         while(node != DenseSet::INVALID_ELEMENT){
            saturation[node]++;
            node = color_saturated_nodes.back().find_next(node);
         }
      }
      uncolored_nodes.remove(best_node);
      nodes--;
   }


#ifndef NDEBUG
   for(const auto& color : coloring.colors()){
      assert(graph.setIsStable(color));
   }
   DenseSet complete_covered(graph.numNodes());
   for(const auto& color : coloring.colors()){
      complete_covered.inplaceUnion(color);
   }
   assert(complete_covered.full());
#endif
   return coloring;
}

NodeColoring GreedyColoring::run_partial_sequential_coloring(const std::vector<Node> &ordering, std::size_t colors) const {
   NodeColoring coloring(graph.numNodes());
   DenseSet free_colors(colors);
   for(const auto& node : ordering){
      free_colors.setAll();
      for(const auto& neighbour : graph.neighbourhood(node)){
         std::size_t color = coloring[neighbour];
         if(color < colors){
            free_colors.remove(color);
         }
      }
      std::size_t first_free_color = free_colors.first();
      coloring[node] = first_free_color;
   }
   coloring.setNumColors(colors);
   return coloring;
}

NodeColoring
GreedyColoring::run_partial_saturation_coloring(const std::vector<Node> &ordering, std::size_t colors) const {

   std::vector<std::size_t> saturation(graph.numNodes(),0);
   std::vector<std::size_t> nodePositions(ordering.size(),0);
   for(std::size_t i = 0; i < ordering.size(); i++){
      nodePositions[ordering[i]] = i;
   }
   DenseSet uncolored_nodes(graph.numNodes(),true);

   NodeColoring coloring(graph.numNodes());
   coloring.setNumColors( colors) ;

   coloring[ordering.front()] = 0;
   std::size_t num_nodes = graph.numNodes()-1;
   DenseSet available_colors(colors);
   while ( num_nodes != 0){
      std::size_t uncolored_node = uncolored_nodes.first();
      std::size_t best_saturation = saturation[uncolored_node];
      std::size_t best_position = nodePositions[uncolored_node];
      std::size_t best_node = uncolored_node;

      while(uncolored_node != DenseSet::INVALID_ELEMENT){
         if(saturation[uncolored_node] > best_saturation ||
             (saturation[uncolored_node] == best_saturation && nodePositions[uncolored_node] > best_position)){
            best_node = uncolored_node;
            best_saturation = saturation[uncolored_node];
            best_position = nodePositions[uncolored_node];
         }
         uncolored_node = uncolored_nodes.find_next(uncolored_node);
      }
      if(best_saturation >= colors){
         coloring[best_node] = INVALID_NODE;
      }else{
         available_colors.setAll();
         for(const auto& node : graph.neighbourhood(best_node)){
            std::size_t color = coloring[node];
            if(color < colors){
               available_colors.remove(color);
            }
         }
         std::size_t best_color = available_colors.first();
         coloring[best_node] = best_color;
      }
      uncolored_nodes.remove(best_node);
      num_nodes--;

   }
   return coloring;
}
}