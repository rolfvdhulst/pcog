//
// Created by rolf on 26-11-22.
//

#include <gtest/gtest.h>
#include <pcog/DenseGraph.hpp>
using namespace pcog;
TEST(DenseGraph,constructor){
   DenseGraph defCons;
   EXPECT_TRUE(defCons.isConsistent());
   EXPECT_EQ(defCons.numNodes(),0);

   DenseGraph four(4);

   EXPECT_TRUE(four.isConsistent());
   EXPECT_EQ(four.numNodes(),4);
   EXPECT_EQ(four.numEdges(),0);

}
TEST(DenseGraph,modification){

   DenseGraph four(4);
   EXPECT_EQ(four.numEdges(),0);

   four.addEdge(1,2);
   EXPECT_EQ(four.numEdges(),1);
   EXPECT_TRUE(four.isEdge(1,2));
   EXPECT_TRUE(four.isEdge(2,1));
   four.addEdge(2,1);
   EXPECT_EQ(four.numEdges(),1);
   EXPECT_TRUE(four.isEdge(1,2));
   EXPECT_TRUE(four.isEdge(2,1));
   EXPECT_EQ(four.numSelfLoops(),0);

   four.addEdge(2,2);
   EXPECT_EQ(four.numSelfLoops(),1);
   EXPECT_TRUE(four.isEdge(2,2));
   four.removeEdge(2,2);
   EXPECT_FALSE(four.isEdge(2,2));
   EXPECT_EQ(four.numEdges(),1);
   EXPECT_EQ(four.numSelfLoops(),0);

   four.clearEdges();
   EXPECT_EQ(four.numEdges(),0);

   DenseSet fiveNeighbourhood(4);
   fiveNeighbourhood.add(0);
   fiveNeighbourhood.add(3);

   Node result = four.addNodeWithNeighbourhood(fiveNeighbourhood);
   EXPECT_EQ(result,4);
   EXPECT_TRUE(four.isEdge(0,4));
   EXPECT_FALSE(four.isEdge(1,4));
   EXPECT_FALSE(four.isEdge(2,4));
   EXPECT_TRUE(four.isEdge(0,4));
   EXPECT_TRUE(four.isConsistent());
}

TEST(DenseGraph, complement){
   DenseGraph graph(3);
   graph.addEdge(0,1);
   graph.addEdge(0,2);

   EXPECT_TRUE(graph.isEdge(0,1));
   EXPECT_TRUE(graph.isEdge(0,2));
   EXPECT_FALSE(graph.isEdge(1,2));
   graph.complement();
   EXPECT_FALSE(graph.isEdge(0,1));
   EXPECT_FALSE(graph.isEdge(0,2));
   EXPECT_TRUE(graph.isEdge(1,2));
}
TEST(DenseGraph,SetIsStable){
   DenseGraph graph(3);
   graph.addEdge(0,1);
   graph.addEdge(0,2);

   DenseSet nonMaximal(3);
   nonMaximal.add(1);

   DenseSet maximal(3);
   maximal.add(1);
   maximal.add(2);

   DenseSet nonStable(3);
   nonStable.add(0);
   nonStable.add(1);

   EXPECT_TRUE(graph.setIsStable(nonMaximal));
   EXPECT_TRUE(graph.setIsStable(maximal));
   EXPECT_FALSE(graph.setIsStable(nonStable));

   EXPECT_TRUE(graph.setIsStableMaximal(maximal));
   EXPECT_FALSE(graph.setIsStableMaximal(nonMaximal));
}

TEST(DenseGraph,InducedSubgraph){
   DenseGraph graph(4);
   graph.addEdge(0,1);
   graph.addEdge(0,3);
   graph.addEdge(1,2);
   graph.addEdge(2,3);

   DenseSet subset(4);
   subset.add(0);
   subset.add(1);
   subset.add(3);

   auto [inducedGraph,map] = graph.nodeInducedSubgraph(subset);

   EXPECT_EQ(map.size(),subset.size());
   EXPECT_EQ(map[0],0);
   EXPECT_EQ(map[1],1);
   EXPECT_EQ(map[2],3);
   EXPECT_TRUE(inducedGraph.isEdge(0,1));
   EXPECT_TRUE(inducedGraph.isEdge(0,2));
   EXPECT_FALSE(inducedGraph.isEdge(1,2));
}

TEST(DenseGraph,Dominate) {
   DenseGraph graph(4);
   graph.addEdge(0, 1);
   graph.addEdge(0, 3);
   graph.addEdge(2, 1);
   graph.addEdge(2, 3);

   EXPECT_TRUE(graph.nodeDominates(0, 2));
   EXPECT_TRUE(graph.nodeDominates(2, 0));
   EXPECT_FALSE(graph.nodeStrictlyDominates(0, 2));
   EXPECT_FALSE(graph.nodeStrictlyDominates(2, 0));
   EXPECT_EQ(graph.nodeDegree(2),2);

   graph.removeEdge(2, 3);

   EXPECT_TRUE(graph.nodeDominates(0, 2));
   EXPECT_FALSE(graph.nodeDominates(2, 0));
   EXPECT_TRUE(graph.nodeStrictlyDominates(0, 2));
   EXPECT_FALSE(graph.nodeStrictlyDominates(2, 0));
   EXPECT_EQ(graph.nodeDegree(2),1);

}
