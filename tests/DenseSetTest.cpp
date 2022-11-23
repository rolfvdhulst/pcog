//
// Created by rolf on 22-11-22.
//

#include <pcog/DenseSet.hpp>
#include "gtest/gtest.h"

TEST(DenseSet, construction1){
   DenseSet set;
   EXPECT_EQ(set.capacity(),0);
}

TEST(DenseSet,construction2){
   DenseSet otherSet(10);
   EXPECT_EQ(otherSet.capacity(),10);
   for (DenseSet::Element i = 0; i < 10; ++i) {
      EXPECT_FALSE(otherSet.contains(i));
   }
   DenseSet secondSet(0);
   EXPECT_EQ(secondSet.capacity(),0);
}

TEST(DenseSet,construction3){
   {
      DenseSet set(10,false);
      EXPECT_EQ(set.capacity(),10);
      for (DenseSet::Element i = 0; i < 10; ++i) {
         EXPECT_FALSE(set.contains(i));
      }
   }
   {
      DenseSet fullSet(10,true);
      EXPECT_EQ(fullSet.capacity(),10);
      for (DenseSet::Element i = 0; i < 10; ++i) {
         EXPECT_TRUE(fullSet.contains(i));
      }
   }

}
TEST(DenseSet,construction4){
   boost::dynamic_bitset bitset;
   bitset.resize(5,false);
   bitset.resize(10,true);
   DenseSet set(bitset);
   for (DenseSet::Element i = 0; i < 10; ++i) {
      if(i < 5){
         EXPECT_FALSE(set.contains(i));
      }else{
         EXPECT_TRUE(set.contains(i));
      }
   }
   //also tests the size functionality
   EXPECT_EQ(set.size(),5);
}

TEST(DenseSet, containmentFunctions){
   DenseSet emptyNoCapacitySet;
   DenseSet emptySet(10);
   DenseSet partialFilledSet(10);
   partialFilledSet.add(5);
   DenseSet fullSet(10,true);

   EXPECT_TRUE(emptyNoCapacitySet.empty());
   EXPECT_TRUE(emptySet.empty());
   EXPECT_FALSE(partialFilledSet.empty());
   EXPECT_FALSE(fullSet.empty());

   EXPECT_FALSE(emptyNoCapacitySet.any());
   EXPECT_FALSE(emptySet.any());
   EXPECT_TRUE(partialFilledSet.any());
   EXPECT_TRUE(fullSet.any());

   EXPECT_TRUE(emptyNoCapacitySet.full());
   EXPECT_FALSE(emptySet.full());
   EXPECT_FALSE(partialFilledSet.full());
   EXPECT_TRUE(fullSet.full());

   EXPECT_TRUE(partialFilledSet.contains(5));
   EXPECT_FALSE(partialFilledSet.contains(4));
}

TEST(DenseSet, setOperationsNonModifying) {
   DenseSet zeroOne(4);
   zeroOne.add(0);
   zeroOne.add(1);

   DenseSet oneTwo(4);
   oneTwo.add(1);
   oneTwo.add(2);

   DenseSet one = zeroOne.intersection(oneTwo);
   EXPECT_FALSE(one.contains(0));
   EXPECT_TRUE(one.contains(1));
   EXPECT_FALSE(one.contains(2));
   EXPECT_FALSE(one.contains(3));

   DenseSet zeroOneTwo = zeroOne.setUnion(oneTwo);
   EXPECT_TRUE(zeroOneTwo.contains(0));
   EXPECT_TRUE(zeroOneTwo.contains(1));
   EXPECT_TRUE(zeroOneTwo.contains(2));
   EXPECT_FALSE(zeroOneTwo.contains(3));

   DenseSet zeroTwo = zeroOne.symmetricDifference(oneTwo);
   EXPECT_TRUE(zeroTwo.contains(0));
   EXPECT_FALSE(zeroTwo.contains(1));
   EXPECT_TRUE(zeroTwo.contains(2));
   EXPECT_FALSE(zeroTwo.contains(3));

   DenseSet zeroOneComp = zeroOne;
   zeroOneComp.complement();
   EXPECT_FALSE(zeroOneComp.contains(0));
   EXPECT_FALSE(zeroOneComp.contains(1));
   EXPECT_TRUE(zeroOneComp.contains(2));
   EXPECT_TRUE(zeroOneComp.contains(3));

   DenseSet zero = zeroOneTwo.difference(oneTwo);

   EXPECT_TRUE(zero.contains(0));
   EXPECT_FALSE(zero.contains(1));
   EXPECT_FALSE(zero.contains(2));
   EXPECT_FALSE(zero.contains(3));
}

TEST(DenseSet, modifyingIntersection) {
   DenseSet zeroOne(4);
   zeroOne.add(0);
   zeroOne.add(1);

   DenseSet oneTwo(4);
   oneTwo.add(1);
   oneTwo.add(2);


   EXPECT_TRUE(zeroOne.contains(0));
   EXPECT_TRUE(zeroOne.contains(1));
   EXPECT_FALSE(zeroOne.contains(2));
   EXPECT_FALSE(zeroOne.contains(3));
   zeroOne.inplaceIntersection(oneTwo);
   EXPECT_FALSE(zeroOne.contains(0));
   EXPECT_TRUE(zeroOne.contains(1));
   EXPECT_FALSE(zeroOne.contains(2));
   EXPECT_FALSE(zeroOne.contains(3));

}

TEST(DenseSet, modifyingUnion) {
   DenseSet zeroOne(4);
   zeroOne.add(0);
   zeroOne.add(1);

   DenseSet oneTwo(4);
   oneTwo.add(1);
   oneTwo.add(2);


   EXPECT_TRUE(zeroOne.contains(0));
   EXPECT_TRUE(zeroOne.contains(1));
   EXPECT_FALSE(zeroOne.contains(2));
   EXPECT_FALSE(zeroOne.contains(3));
   zeroOne.inplaceUnion(oneTwo);
   EXPECT_TRUE(zeroOne.contains(0));
   EXPECT_TRUE(zeroOne.contains(1));
   EXPECT_TRUE(zeroOne.contains(2));
   EXPECT_FALSE(zeroOne.contains(3));

}
TEST(DenseSet, modifyingSymmetricDifference) {
   DenseSet zeroOne(4);
   zeroOne.add(0);
   zeroOne.add(1);

   DenseSet oneTwo(4);
   oneTwo.add(1);
   oneTwo.add(2);


   EXPECT_TRUE(zeroOne.contains(0));
   EXPECT_TRUE(zeroOne.contains(1));
   EXPECT_FALSE(zeroOne.contains(2));
   EXPECT_FALSE(zeroOne.contains(3));
   zeroOne.inplaceSymmetricDifference(oneTwo);
   EXPECT_TRUE(zeroOne.contains(0));
   EXPECT_FALSE(zeroOne.contains(1));
   EXPECT_TRUE(zeroOne.contains(2));
   EXPECT_FALSE(zeroOne.contains(3));

}
TEST(DenseSet, modifyingDifference) {
   DenseSet zeroOne(4);
   zeroOne.add(0);
   zeroOne.add(1);

   DenseSet oneTwo(4);
   oneTwo.add(1);
   oneTwo.add(2);


   EXPECT_TRUE(zeroOne.contains(0));
   EXPECT_TRUE(zeroOne.contains(1));
   EXPECT_FALSE(zeroOne.contains(2));
   EXPECT_FALSE(zeroOne.contains(3));
   zeroOne.inplaceDifference(oneTwo);
   EXPECT_TRUE(zeroOne.contains(0));
   EXPECT_FALSE(zeroOne.contains(1));
   EXPECT_FALSE(zeroOne.contains(2));
   EXPECT_FALSE(zeroOne.contains(3));

}

TEST(DenseSet, clear) {
   DenseSet zeroOne(4);
   zeroOne.add(0);
   zeroOne.add(1);

   zeroOne.clear();
   EXPECT_EQ(zeroOne.size(),0);
   EXPECT_EQ(zeroOne.capacity(),4);

   DenseSet empty;
   empty.clear();
   EXPECT_EQ(empty.size(),0);
   EXPECT_EQ(empty.capacity(),0);
}

TEST(DenseSet,Equality){
   DenseSet set1(4);
   set1.add(0);
   set1.add(1);

   DenseSet set2(4);
   set2.add(0);
   set2.add(1);
   set2.add(2);

   DenseSet set3(4);
   set3.add(0);
   set3.add(1);

   EXPECT_TRUE(set1 == set3);
   EXPECT_TRUE(set3 == set1);
   EXPECT_FALSE(set1 != set3);
   EXPECT_FALSE(set3 != set1);

   EXPECT_FALSE(set1 == set2);
   EXPECT_FALSE(set2 == set1);
   EXPECT_TRUE(set1 != set2);
   EXPECT_TRUE(set2 != set1);

   DenseSet empty1;
   DenseSet empty2;

   EXPECT_TRUE(empty1 == empty2);
}

TEST(DenseSet,modification){
   DenseSet set(3);
   EXPECT_FALSE(set.contains(0));
   EXPECT_FALSE(set.contains(1));
   EXPECT_FALSE(set.contains(2));
   set.add(1);
   EXPECT_FALSE(set.contains(0));
   EXPECT_TRUE(set.contains(1));
   EXPECT_FALSE(set.contains(2));
   set.flip(0);
   set.flip(1);
   EXPECT_TRUE(set.contains(0));
   EXPECT_FALSE(set.contains(1));
   EXPECT_FALSE(set.contains(2));

   set.setAll();
   EXPECT_TRUE(set.full());
   EXPECT_TRUE(set.contains(0));
   EXPECT_TRUE(set.contains(1));
   EXPECT_TRUE(set.contains(2));

   set.remove(1);
   EXPECT_TRUE(set.contains(0));
   EXPECT_FALSE(set.contains(1));
   EXPECT_TRUE(set.contains(2));

   DenseSet emptySet(3);

   emptySet.setRange(0,0);
   EXPECT_TRUE(emptySet.contains(0));
   EXPECT_FALSE(emptySet.contains(1));
   EXPECT_FALSE(emptySet.contains(2));
   emptySet.setRange(0,1);
   EXPECT_TRUE(emptySet.contains(0));
   EXPECT_TRUE(emptySet.contains(1));
   EXPECT_FALSE(emptySet.contains(2));
   emptySet.extend(false);
   EXPECT_EQ(emptySet.capacity(),4);
   EXPECT_FALSE(emptySet.contains(3));
   emptySet.setRange(1,3);
   EXPECT_TRUE(emptySet.full());

   emptySet.set(2,false);
   EXPECT_FALSE(emptySet.contains(2));
   emptySet.set(2,true);
   EXPECT_TRUE(emptySet.contains(2));

}

TEST(DenseSet,subset){
   DenseSet set1(3);
   set1.add(0);
   DenseSet set2(3);
   set2.add(0);
   set2.add(1);

   DenseSet set3(3);
   set3.add(0);

   EXPECT_TRUE(set1.isSubsetOf(set3));
   EXPECT_TRUE(set3.isSubsetOf(set1));
   EXPECT_TRUE(set1.isSubsetOf(set2));
   EXPECT_FALSE(set2.isSubsetOf(set1));

   EXPECT_TRUE(set1.hasAsSubset(set3));
   EXPECT_TRUE(set3.hasAsSubset(set1));
   EXPECT_FALSE(set1.hasAsSubset(set2));
   EXPECT_TRUE(set2.hasAsSubset(set1));

   EXPECT_FALSE(set1.isProperSubsetOf(set3));
   EXPECT_FALSE(set3.isProperSubsetOf(set1));
   EXPECT_TRUE(set1.isProperSubsetOf(set2));
   EXPECT_FALSE(set2.isProperSubsetOf(set1));

   EXPECT_FALSE(set1.hasAsProperSubset(set3));
   EXPECT_FALSE(set3.hasAsProperSubset(set1));
   EXPECT_FALSE(set1.hasAsProperSubset(set2));
   EXPECT_TRUE(set2.hasAsProperSubset(set1));

   DenseSet empty1;
   DenseSet empty2;
   //TODO: what should the behavior of the empty set be?
}

TEST(DenseSet,iteration){
   std::vector<DenseSet::Element> numbers = {2,3,5,7,8,9};

   DenseSet set(10);
   for(const auto& number : numbers){
      set.add(number);
   }

   EXPECT_EQ(set.first(),2);
   EXPECT_EQ(set.find_next(2),3);
   EXPECT_EQ(set.find_next(3),5);
   EXPECT_EQ(set.find_next(4),5);

   std::size_t index = 0;
   for(const auto& elem : set){
      EXPECT_EQ(elem,numbers[index]);
      ++index;
   }
   DenseSet empty;
   EXPECT_EQ(empty.begin(),empty.end());
}