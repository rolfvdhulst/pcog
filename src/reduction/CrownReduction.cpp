//
// Created by rolf on 15-7-23.
//

#include "pcog/reduction/CrownReduction.hpp"
#include "pcog/utilities/DenseSet.hpp"
#include "pcog/reduction/ReductionStack.hpp"


namespace pcog {


bool coverDfs(std::size_t v,
         const std::vector<std::vector<std::size_t>>& adj,
         const std::vector<int>& x,
         const std::vector<std::size_t>& level,
         std::vector<int>& iter,
         std::vector<int>& in,
         std::vector<int>& out){
   while (iter[v] >= 0) {
      std::size_t u = adj[v][iter[v]--];
      if (x[u] >= 0) continue;
      int w = in[u];
      if (w < 0 || (level[v] < level[w] && iter[w] >= 0 && coverDfs(w,adj,x,level,iter,in,out))) {
         in[u] = static_cast<int>(v);
         out[v] = static_cast<int>(u);
         return true;
      }
   }
   return false;
}

std::vector<std::vector<std::size_t>> computeCrownCliques(std::vector<int>& x,
                        const std::vector<int>& out,
                        const std::vector<std::vector<std::size_t>>& adj,
                                                         const DenseSet& used){
   std::vector<std::vector<std::size_t>> crownCliques;
   std::size_t n = x.size();
   for (std::size_t v = 0; v <n; v++) {
      if (x[v] < 0 && used.contains(v) && !used.contains(n+ v)) {
         //set x
         x[v] = 0;
         for(std::size_t u : adj[v]){
            if(x[u] < 0){
               x[u] = 1;
            }
         }
         std::vector<std::size_t> clique;
         if (out[v] != -1){
            clique.push_back(v);
            clique.push_back(out[v]);
         }else {
            clique.push_back(v);
         }
         crownCliques.push_back(clique);
      }
   }

   return crownCliques;
}
std::vector<std::vector<std::size_t>> findCrownCliques(std::size_t n,
              const std::vector<std::vector<std::size_t>>& adj) {
   // solution: -1 : undetermined, 0: not in the cover, 1 in the cover
   std::vector<int> x(n, -1);
   std::vector<int> in(n, -1);
   std::vector<int> out(n, -1);
   std::vector<std::size_t> queue(2 * n, 0);
   std::vector<std::size_t> level(2 * n, 0);
   std::vector<int> iter(2 * n, 0);
   DenseSet used(2 * n, false);

   for (;;) {
      used.clear();
      std::size_t qs = 0;
      std::size_t qt = 0;
      for (std::size_t v = 0; v < n; v++) {
         if (x[v] < 0 && out[v] < 0) {
            level[v] = 0;
            used.add(v);
            queue[qt++] = v;
         }
      }
      bool ok = false;
      while (qs < qt) {
         std::size_t v = queue[qs++];
         iter[v] = static_cast<int>(adj[v].size()) - 1;
         for (std::size_t u : adj[v]){
            if(x[u] >= 0) continue;
            if(used.contains(n+u)) continue;
            used.add(n+u);
            int w = in[u];
            if (w < 0) {
               ok = true;
            } else {
               level[w] = level[v] + 1;
               used.add(w);
               queue[qt++] = w;
            }
         }
      }
      if (!ok) break;
      for (int v = static_cast<int>(n - 1); v >= 0; v--) {
         if (x[v] < 0 && out[v] < 0) {
            coverDfs(v,adj,x,level,iter,in,out);
         }
      }
   }

   return computeCrownCliques(x,out,adj,used);
}
bool findCrownReductions(DenseReductionGraph& graph,
                         ReductionStack& stack,
                         ReductionVertexQueue& queue
){
   //build the complement graph of the current graph
   std::vector<std::size_t> originalToNew(graph.nodes().capacity(),INVALID_NODE);
   std::vector<std::size_t> newToOriginal;
   for(Node node : graph.nodes()){
      newToOriginal.push_back(node);
   }
   for(std::size_t index = 0; index < newToOriginal.size(); ++index){
      originalToNew[newToOriginal[index]] = index;
   }
   std::size_t numNodes = newToOriginal.size();
   std::vector<std::vector<std::size_t>> adj(numNodes,std::vector<std::size_t>());
   std::size_t index = 0;
   for(Node node : graph.nodes()){
      for(const auto& neighbour : graph.complementNeighbourhood(node)){
         assert(originalToNew[neighbour] != INVALID_NODE);
         adj[index].push_back(originalToNew[neighbour]);
      }
      ++index;
   }

   std::vector<std::vector<std::size_t>> crownCliques = findCrownCliques(adj.size(),adj);
   if(crownCliques.empty()){
      return false;
   }

   CrownReduction reduction;
   for(const auto& clique : crownCliques){
      DenseSet cliqueSet(graph.nodes().capacity(),false);
      for(const auto& cliqueNode : clique){
         cliqueSet.add(newToOriginal[cliqueNode]);
      }
      reduction.fixedSets.push_back(cliqueSet);
   }
   //Assert that we have truly found a correct crown reduction:
#ifndef NDEBUG
   //asserts sets are indeed stable
   for(const auto& set : reduction.fixedSets){
      if(set.size() == 2){
         Node firstNode = set.first();
         Node secondNode = set.find_next(firstNode);
         assert(firstNode != INVALID_NODE && secondNode != INVALID_NODE);
         assert(!graph.neighbourhood(firstNode).contains(secondNode));
      }else{
         assert(set.size() == 1);
      }
   }
   //Find back head and crown
   DenseSet head(graph.nodes().capacity());
   DenseSet crown(graph.nodes().capacity());
   DenseSet crownNeighbourhood(graph.nodes().capacity());
   for(const auto& clique : crownCliques){
      if(clique.size() == 2) {
         Node crownNode = newToOriginal[clique[0]];
         Node headNode = newToOriginal[clique[1]];

         crown.add(crownNode);
         crownNeighbourhood.inplaceUnion(graph.complementNeighbourhood(crownNode));

         head.add(headNode);
      }else{
         Node node = newToOriginal[clique[0]];
         crown.add(node);
         crownNeighbourhood.inplaceUnion(graph.complementNeighbourhood(node));
      }
   }

   //check that the crown is a clique
   for(const auto& crownNode : crown){
      DenseSet set = crown;
      set.remove(crownNode);
      assert(set.isSubsetOf(graph.neighbourhood(crownNode)));
   }
   assert(crownNeighbourhood == head);
   assert(crown.size() >= head.size());
#endif


   stack.push(reduction);
   for(const auto& fixedSet : reduction.fixedSets){
      graph.removeStableSet(fixedSet);
   }
   for(Node node : graph.nodes()){
      queue.push(node,graph.lowerBoundNodes().contains(node));
   }
   return true;
}
void CrownReduction::transformStableSet(DenseSet &) const {
   //no-op, since the fixed sets are removed
}
void CrownReduction::newToOldColoring(NodeColoring &coloring) const {
   std::size_t count = coloring.numColors();
   for(const auto& set : fixedSets){
      for(Node node : set){
         assert(coloring[node] == INVALID_COLOR);
         coloring[node] = count;
      }
      ++count;
   }
   coloring.setNumColors(count);
}
}