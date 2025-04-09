//
// Created by rolf on 12-7-23.
//

#include "pcog/reduction/TwinDegreeThreeReduction.hpp"
#include "pcog/reduction/ReductionStack.hpp"

namespace pcog {
bool twinDegreeThreeReduction(Node node, DenseReductionGraph &graph,
                              ReductionStack &stack,
                              ReductionVertexQueue &queue) {
   if(graph.nodeDegree(node) != (graph.numNodes() - 4)){
      return false;
   }
   DenseSet nonNeighbours = graph.complementNeighbourhood(node);

   Node w = nonNeighbours.first();
   Node x = nonNeighbours.find_next(w);
   Node y = nonNeighbours.find_next(x);
   assert(w != INVALID_NODE && x != INVALID_NODE && y != INVALID_NODE);
   assert(nonNeighbours.find_next(y) == INVALID_NODE);

   DenseSet twinCandidates = graph.neighbourhood(node).difference(graph.neighbourhood(w));
   twinCandidates.inplaceDifference(graph.neighbourhood(x));
   twinCandidates.inplaceDifference(graph.neighbourhood(y));

   Node twin = INVALID_NODE;
   for(Node candidate : twinCandidates){
      if(graph.nodeDegree(candidate) == graph.numNodes()-4){
         twin = candidate;
         break;
      }
   }
   if(twin == INVALID_NODE){
      return false;
   }
   bool edgeWX = graph.neighbourhood(w).contains(x);
   bool edgeWY = graph.neighbourhood(w).contains(y);
   bool edgeXY = graph.neighbourhood(x).contains(y);
   if(!(edgeWX && edgeWY && edgeXY)){
      //There is at least one non-edge, we can remove node twin and wxy from the graph by fixing two stable sets
      DenseSet firstSet(graph.nodes().capacity());
      firstSet.add(node);
      firstSet.add(w);
      if(!edgeWX){
         firstSet.add(x);
      }
      if(!edgeWY){
         firstSet.add(y);
      }
      DenseSet secondSet(graph.nodes().capacity());
      secondSet.add(twin);
      if(!edgeWX){
         secondSet.add(y);
      }else if(!edgeWY){
         secondSet.add(x);
      }else{
         assert(!edgeXY);
         secondSet.add(x);
         secondSet.add(y);
      }

      assert(firstSet.setUnion(secondSet).size() >= 5);
      graph.removeStableSet(firstSet);
      graph.removeStableSet(secondSet);

      stack.push(TwinDegreeThreeReduction{.u= node,.v = twin,.uSet = firstSet,.vSet = secondSet});
      //TODO: what nodes to push to the queue? (perhaps intersection of neighbourhoods? or non-neighbourhoods)
      for(const auto& addNode : graph.nodes()){
         queue.push(addNode,graph.lowerBoundNodes().contains(addNode));
      }

      return true;
   }
   //Perform a folding reduction. Similar to degree 2 reduction, pick which nodes to remove to preserve clique lb;
   //If u and/or v are in clique, then no problems; folding does not destroy the clique. We arbitrarily pick from wxy
   //If one or two nodes of wxy are in the clique, then remove these clique nodes first (otherwise, clique lb might be invalidated)
   //If wxy are all in the clique, we can arbitrarily remove two of them without any problems.

   Node keepNode,removeNode1,removeNode2;
   if(!graph.lowerBoundNodes().contains(x)){
      //Fold neighbourhoods, keeping x
      keepNode = x;
      removeNode1 = w;
      removeNode2 = y;
   }else if(!graph.lowerBoundNodes().contains(y)) {
      //Fold neighbourhoods keeping y
      keepNode = y;
      removeNode1 = w;
      removeNode2 = x;
   }else{
      //keep w and fold neighbourhoods
      keepNode = w;
      removeNode1 = x;
      removeNode2 = y;
   }
   DenseSet keepRemoveNeighbours = graph.neighbourhood(keepNode);
   DenseSet node1NonNeighbours = graph.complementNeighbourhood(removeNode1);
   DenseSet node2NonNeighbours = graph.complementNeighbourhood(removeNode2);
   DenseSet possibleRemoveEdges = node1NonNeighbours.setUnion(node2NonNeighbours);
   keepRemoveNeighbours.inplaceIntersection(possibleRemoveEdges);

   DenseSet keepNonNeighbours = graph.complementNeighbourhood(keepNode);
   //TODO: remove as stable sets, or is this good enough?
   graph.removeNode(node);
   graph.removeNode(twin);
   graph.removeNode(removeNode1);
   graph.removeNode(removeNode2);
   graph.removeNodeEdges(keepNode,keepRemoveNeighbours);

   stack.push(TwinDegreeThreeFoldReduction{
       .u= node,.v = twin,.keep = keepNode,.remove1 = removeNode1, .remove2 = removeNode2,
       .keepNonNeighbours = keepNonNeighbours,.remove1NonNeighbours = node1NonNeighbours,
       .remove2NonNeighbours = node2NonNeighbours,
   });
   //TODO: what nodes to push to the queue? (perhaps intersection of neighbourhoods? or non-neighbourhoods)
   for(const auto& addNode : graph.nodes()){
      queue.push(addNode,graph.lowerBoundNodes().contains(addNode));
   }

   return true;
}
void TwinDegreeThreeReduction::transformStableSet(DenseSet &) const {
   //no-op
}
void TwinDegreeThreeReduction::newToOldColoring(NodeColoring &coloring) const {
   std::size_t count = coloring.numColors();
   for (Node fixed : uSet) {
      assert(coloring[fixed] == INVALID_COLOR);
      coloring[fixed] = count;
   }
   ++count;
   for (Node fixed : vSet) {
      assert(coloring[fixed] == INVALID_COLOR);
      coloring[fixed] = count;
   }
   ++count;
   coloring.setNumColors(count);
}
void TwinDegreeThreeFoldReduction::transformStableSet(DenseSet & /*set*/) const {
   //TODO: fix, not complete
}
void TwinDegreeThreeFoldReduction::newToOldColoring(NodeColoring &coloring) const {

   assert(coloring[u] == INVALID_COLOR);
   assert(coloring[v] == INVALID_COLOR);
   assert(coloring[remove1] == INVALID_COLOR);
   assert(coloring[remove2] == INVALID_COLOR);

   assert(coloring[keep] != INVALID_COLOR);

   Color keptNodeColor = coloring[keep];

   assert(keptNodeColor != INVALID_COLOR);

   bool inRemoveOneNonNeighbours = true;
//   bool inRemoveTwoNonNeighbours = true;
   bool inKeptNonNeighbours = true;
   for(std::size_t i = 0; i < coloring.numNodes(); ++i){
      if(i == keep) continue;
      if(keptNodeColor == coloring[i] ){
         if(!remove1NonNeighbours.contains(i)){
            inRemoveOneNonNeighbours = false;
         }
         if(!remove2NonNeighbours.contains(i)){
//            inRemoveTwoNonNeighbours = false;
         }
         if(!keepNonNeighbours.contains(i)){
            inKeptNonNeighbours = false;
         }
         //TODO: add these lines after extensive testing
         // they should provide a minor speedup but if there's a bug they also might cause us not to catch it as easily
//         if((inRemoveOneNonNeighbours ? 1 : 0) + (inRemoveTwoNonNeighbours ? 1 : 0) + (inKeptNonNeighbours ? 1 : 0) <= 1){
//            break;
//         }
      }
   }
   std::size_t count = coloring.numColors();

   if(inKeptNonNeighbours){
      coloring[u] = count;
      coloring[remove1] = count;
      coloring[v] = count+1;
      coloring[remove2] = count+1;
   }else if(inRemoveOneNonNeighbours){
      coloring[remove1] = coloring[keep];
      coloring[u] = count;
      coloring[keep] = count;
      coloring[v] = count+1;
      coloring[remove2] = count+1;
   }else{
      assert(inRemoveTwoNonNeighbours);
      coloring[remove2] = coloring[keep];
      coloring[u] = count;
      coloring[remove1] = count;
      coloring[v] = count+1;
      coloring[keep] = count+1;
   }

   coloring.setNumColors(count+2);
}
}