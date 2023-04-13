//
// Created by rolf on 26-11-22.
//

#include <gtest/gtest.h>
#include <pcog/utilities/Coloring.hpp>
#include <pcog/utilities/DenseGraph.hpp>

using namespace pcog;
TEST(SetColoring,construction){
   SetColoring coloring;
   EXPECT_EQ(coloring.numColors(),0);
   EXPECT_TRUE(coloring.colors().empty());
}
TEST(SetColoring,addColor){
   SetColoring coloring;
   EXPECT_EQ(coloring.numColors(),0);
   EXPECT_TRUE(coloring.colors().empty());
   DenseGraph emptyGraph(3);
   DenseSet set(3);
   set.add(1);
   set.add(2);
   coloring.addColor(set);
   EXPECT_EQ(coloring.numColors(),1);
   EXPECT_EQ(coloring.colors().front(),set);
   EXPECT_FALSE(coloring.isValid(emptyGraph));
   DenseSet set2(3);
   set2.add(0);
   coloring.addColor(set2);
   EXPECT_EQ(coloring.numColors(),2);
   EXPECT_EQ(coloring.colors().back(),set2);
   EXPECT_TRUE(coloring.isValid(emptyGraph));
}
TEST(SetColoring, constructionNodeColoring){
   NodeColoring nodeCol(4);
   nodeCol[0] = Color(1);
   nodeCol[1] = Color(0);
   nodeCol[2] = Color(0);
   nodeCol[3] = Color(2);
   nodeCol.setNumColors(3);

   SetColoring setColoring(nodeCol);
   EXPECT_EQ(setColoring.numColors(),3);
   EXPECT_FALSE(setColoring.colors()[0].contains(0));
   EXPECT_TRUE(setColoring.colors()[0].contains(1));
   EXPECT_TRUE(setColoring.colors()[0].contains(2));
   EXPECT_FALSE(setColoring.colors()[0].contains(3));

   EXPECT_TRUE(setColoring.colors()[1].contains(0));
   EXPECT_FALSE(setColoring.colors()[1].contains(1));
   EXPECT_FALSE(setColoring.colors()[1].contains(2));
   EXPECT_FALSE(setColoring.colors()[1].contains(3));

   EXPECT_FALSE(setColoring.colors()[2].contains(0));
   EXPECT_FALSE(setColoring.colors()[2].contains(1));
   EXPECT_FALSE(setColoring.colors()[2].contains(2));
   EXPECT_TRUE(setColoring.colors()[2].contains(3));
}
TEST(NodeColoring, Construction){
   NodeColoring coloring(4);
   EXPECT_EQ(coloring.numUncolored(),4);
   EXPECT_EQ(coloring.numColors(),0);
}

TEST(NodeColoring,ConstructFromSetColoring){
   SetColoring coloring;
   DenseSet set(3);
   set.add(1);
   set.add(2);
   coloring.addColor(set);
   set.complement();
   coloring.addColor(set);

   NodeColoring nodeColoring(3,coloring);
   EXPECT_EQ(nodeColoring.numColors(),2);
   EXPECT_EQ(nodeColoring[0],1);
   EXPECT_EQ(nodeColoring[1],0);
   EXPECT_EQ(nodeColoring[2],0);
   EXPECT_EQ(nodeColoring.numUncolored(),0);
   EXPECT_TRUE(nodeColoring.hasNoInvalidNodes());
}
TEST(NodeColoring, setNumColors){
   NodeColoring coloring(4);
   EXPECT_EQ(coloring.numColors(),0);
   coloring.setNumColors(3);
   EXPECT_EQ(coloring.numColors(),3);

}
TEST(NodeColoring, NoInvalidNodes){
   NodeColoring coloring(2);
   coloring.setNumColors(2);
   EXPECT_FALSE(coloring.hasNoInvalidNodes());
   coloring[0] = 0;
   EXPECT_FALSE(coloring.hasNoInvalidNodes());
   coloring[1] = 1;
   EXPECT_TRUE(coloring.hasNoInvalidNodes());
   coloring[0] = INVALID_COLOR;
   EXPECT_FALSE(coloring.hasNoInvalidNodes());

}
TEST(NodeColoring,StableSetSizes){
   NodeColoring coloring(4);
   coloring.setNumColors(2);
   coloring[0] = 1;
   coloring[1] = 1;

   auto sizes =coloring.getStableSetSizes();
   EXPECT_EQ(sizes[0],0);
   EXPECT_EQ(sizes[1],2);
}
TEST(NodeColoring,NumUncolored){
   NodeColoring coloring(4);
   coloring.setNumColors(2);
   EXPECT_EQ(coloring.numUncolored(),4);
   coloring[0] = 1;
   EXPECT_EQ(coloring.numUncolored(),3);

}