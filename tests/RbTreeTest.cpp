//
// Created by rolf on 8-4-23.
//
#include <pcog/utilities/RbTree.hpp>
#include "gtest/gtest.h"
#include <numeric>
#include <random>

using namespace pcog;

struct Node{
   int64_t key;
   RbTreeLinks<int64_t> links;
};
class MyTree;
template<>
struct RbTreeTraits<MyTree>{
   using KeyType = int64_t;
   using LinkType = int64_t;
};

class MyTree : public RbTree<MyTree> {
 public:
   int64_t root = -1;
   std::vector<Node> nodes;
   MyTree() : RbTree<MyTree>(root) {};
   RbTreeLinks<int64_t>& getRbTreeLinks(int64_t node){
      return nodes[node].links;
   }
   [[nodiscard]] const RbTreeLinks<int64_t>& getRbTreeLinks(int64_t node) const{
      return nodes[node].links;
   }
   [[nodiscard]] int64_t getKey(int64_t node) const {return nodes[node].key;}

   bool test_insert(int64_t x) {
      std::pair<int64_t, bool> p = find(x);

      if (p.second)
         return false;

      int64_t newNode = nodes.size();
      nodes.emplace_back();
      nodes.back().key = x;
      insert(newNode, p.first);
      return true;
   }

   void test_erase(int64_t node) { erase(node); }

   bool contains(int64_t x) { return find(x).second; }

};
static void checkRbTree(MyTree& tree, int64_t* expectedKeys,
                        int64_t numExpectedKeys) {
   std::vector<int64_t> keys;
   keys.reserve(numExpectedKeys);
   if (tree.root != -1) EXPECT_TRUE(tree.nodes[tree.root].links.isBlack());

   int64_t x = tree.lower_bound();
   while (x != -1) {
      keys.push_back(tree.nodes[x].key);
      if (tree.nodes[x].links.isRed()) {
         int64_t lChild = tree.nodes[x].links.child[0];
         int64_t rChild = tree.nodes[x].links.child[1];
         if (lChild != -1) EXPECT_TRUE(tree.nodes[lChild].links.isBlack());
         if (rChild != -1) EXPECT_TRUE(tree.nodes[rChild].links.isBlack());
      }
      x = tree.successor(x);
      EXPECT_TRUE((int64_t)keys.size() <= numExpectedKeys);
   }

   EXPECT_TRUE((int64_t)keys.size() == numExpectedKeys);
   std::sort(expectedKeys, expectedKeys + numExpectedKeys);
   bool isOk = std::equal(keys.begin(), keys.end(), expectedKeys);
   EXPECT_TRUE(isOk);
}

TEST(RbTree, general) {
   std::vector<int64_t> keys;
   keys.resize(1000);
   std::iota(keys.begin(), keys.end(), 0);

   std::random_device device;
   std::shuffle(keys.begin(),keys.end(),device);
   MyTree rbTree;

   for (size_t i = 0; i < keys.size(); ++i) {
      int64_t x = keys[i];
      bool inserted = rbTree.test_insert(x);
      EXPECT_TRUE(inserted);
      checkRbTree(rbTree, keys.data(), i + 1);
   }

   // randomly delete half of the elements and check the tree after each deletion

   for (size_t i = keys.size() - 1; i > keys.size() / 2; --i) {
      std::uniform_int_distribution distribution(size_t(0),i+1);
      size_t k = distribution(device);
      std::swap(keys[k], keys[i]);
      int64_t x = keys[i];
      std::pair<int64_t, bool> node = rbTree.find(x);
      EXPECT_TRUE(node.second);
      rbTree.test_erase(node.first);
      checkRbTree(rbTree, keys.data(), i);
   }
}