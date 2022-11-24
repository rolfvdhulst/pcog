//
// Created by rolf on 24-11-22.
//

#include <pcog/NodeMap.hpp>
#include <pcog/DenseSet.hpp>
#include "gtest/gtest.h"


TEST(NodeMap, construction){
   NodeMap map;
   EXPECT_EQ(map.size(),0);
   NodeMap secondMap(4);
   EXPECT_EQ(secondMap.size(),4);
   EXPECT_EQ(secondMap[0],INVALID_NODE);
   EXPECT_EQ(secondMap[1],INVALID_NODE);
   EXPECT_EQ(secondMap[2],INVALID_NODE);
   EXPECT_EQ(secondMap[3],INVALID_NODE);
}
TEST(NodeMap,extend){
   NodeMap map;
   map.extend(1);
   map.extend(2);
   map.extend(0);
   EXPECT_EQ(map.size(),3);
   EXPECT_EQ(map[0],1);
   EXPECT_EQ(map[1],2);
   EXPECT_EQ(map[2],0);

   NodeMap inverse = NodeMap::inverse(map,3);
   EXPECT_EQ(inverse[1],0);
   EXPECT_EQ(inverse[2],1);
   EXPECT_EQ(inverse[0],2);

}

TEST(NodeMap,identity){
   NodeMap map;
   EXPECT_TRUE(map.isIdentity());
   map.extend(0);
   EXPECT_TRUE(map.isIdentity());
   map.extend(2);
   EXPECT_FALSE(map.isIdentity());
   map.setIdentity();
   EXPECT_TRUE(map.isIdentity());
   EXPECT_EQ(map[1],1);

}

TEST(NodeMap,transform){
   NodeMap map;
   map.extend(0);
   map.extend(2);
   map.extend(3);
   map.extend(1);

   DenseSet set(4);
   set.add(0);
   set.add(1);

   DenseSet output(4);

   map.transform(set,output);

   EXPECT_TRUE(output.contains(0));
   EXPECT_FALSE(output.contains(1));
   EXPECT_TRUE(output.contains(2));
   EXPECT_FALSE(output.contains(3));
}