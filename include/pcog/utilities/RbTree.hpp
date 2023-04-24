//
// Created by rolf on 8-4-23.
//

#ifndef PCOG_RBTREE_HPP
#define PCOG_RBTREE_HPP

#include <array>
#include <cstdint>
#include <cassert>

namespace pcog {

// For now we only allow T to be a signed integral type:
// In the future one could also allow T to be a pointer type, but this would
// complicate the code
template <typename T>
concept RbTreeIndex = std::is_signed<T>::value;

template <typename T> struct RbTreeTraits;

template <RbTreeIndex T> struct RbTreeLinks {

   std::array<T, 2> child = {noLink(),noLink()};

   static constexpr T noLink() { return -1; }
   using ParentStorageType = typename std::make_unsigned<T>::type;
   ParentStorageType parentAndColor = 0;

   static constexpr int colorBitPos() {
      return sizeof(ParentStorageType) * 8 - 1;
   }
   static constexpr ParentStorageType colorBitMask() {
      return ParentStorageType{1} << colorBitPos();
   }
   void makeRed() { parentAndColor |= colorBitMask(); }
   void makeBlack() { parentAndColor &= ~colorBitMask(); }
   ParentStorageType getColor() const {
      return parentAndColor >> colorBitPos();
   }
   void setColor(bool isRed){
      makeBlack();
      parentAndColor |= ParentStorageType(isRed) << colorBitPos();
   }
   [[nodiscard]] bool isBlack() const { return getColor() == 0; }
   [[nodiscard]] bool isRed() const { return getColor() == 1; }

   void setParent(T p) {
      assert(p >= -1);
      // keeping the color the same requires that we preserve the last bit
      parentAndColor = (parentAndColor & colorBitMask()) |
                       static_cast<ParentStorageType>(p + 1);
   }
   T getParent() const {
      return static_cast<T>(parentAndColor & ~colorBitMask()) - 1;
   }
};

template <typename T> class RbTree {
   ///All keys in left subtree are smaller, all keys in right subtree are larger
   enum Direction { LEFT = 0, RIGHT = 1 };
   using KeyType = typename RbTreeTraits<T>::KeyType;
   using LinkType = typename RbTreeTraits<T>::LinkType;
   using ColorType = typename RbTreeLinks<LinkType >::ParentStorageType;

   LinkType &rootNode;
   LinkType &lowerNode;

   static constexpr LinkType NO_LINK = RbTreeLinks<LinkType>::noLink();

   bool isRed(LinkType node) const {
      return node != NO_LINK &&
             static_cast<const T *>(this)->getRbTreeLinks(node).isRed();
   }

   bool isBlack(LinkType node) const {
      return node == NO_LINK ||
             static_cast<const T *>(this)->getRbTreeLinks(node).isBlack();
   }

   void makeRed(LinkType node) {
      static_cast<T *>(this)->getRbTreeLinks(node).makeRed();
   }

   void makeBlack(LinkType node) {
      static_cast<T *>(this)->getRbTreeLinks(node).makeBlack();
   }

   void setColor(LinkType node, bool isRed) {
      static_cast<T*>(this)->getRbTreeLinks(node).setColor(isRed);
   }

   ColorType getColor(LinkType node) const {
      return static_cast<const T*>(this)->getRbTreeLinks(node).getColor();
   }

   void setParent(LinkType node, LinkType parent) {
      static_cast<T*>(this)->getRbTreeLinks(node).setParent(parent);
   }

   LinkType getParent(LinkType node) const {
      return static_cast<const T*>(this)->getRbTreeLinks(node).getParent();
   }

   KeyType getKey(LinkType node) const {
      return static_cast<const T*>(this)->getKey(node);
   }

   static constexpr Direction opposite(Direction dir) { return Direction(1 - dir); }

   LinkType getChild(LinkType node, Direction dir) const {
      return static_cast<const T*>(this)->getRbTreeLinks(node).child[dir];
   }

   void setChild(LinkType node, Direction dir, LinkType child) {
      static_cast<T*>(this)->getRbTreeLinks(node).child[dir] = child;
   }

   ///Rotates the subtree around node x in the given direction
   void rotate(LinkType x, Direction dir) {
      LinkType y = getChild(x, opposite(dir));
      LinkType yDir = getChild(y, dir);
      setChild(x, opposite(dir), yDir);
      if (yDir != NO_LINK) setParent(yDir, x);

      LinkType pX = getParent(x);
      setParent(y, pX);

      if (pX == NO_LINK){
         rootNode = y;
      }else{
         setChild(pX, Direction((x != getChild(pX, dir)) ^ dir), y);
      }

      setChild(y, dir, x);
      setParent(x, y);
   }
   ///Fixes up the red-black property of the tree after a node has been inserted.
   void insertFixup(LinkType z) {
      LinkType pZ = getParent(z);
      while (isRed(pZ)) {
         LinkType zGrandParent = getParent(pZ);
         assert(zGrandParent != NO_LINK);

         Direction dir = Direction(getChild(zGrandParent, LEFT) == pZ);

         LinkType y = getChild(zGrandParent, dir);
         if (isRed(y)) {
            makeBlack(pZ);
            makeBlack(y);
            makeRed(zGrandParent);
            z = zGrandParent;
         } else {
            if (z == getChild(pZ, dir)) {
               z = pZ;
               rotate(z, opposite(dir));
               pZ = getParent(z);
               zGrandParent = getParent(pZ);
               assert(zGrandParent != NO_LINK);
            }

            makeBlack(pZ);
            makeRed(zGrandParent);
            rotate(zGrandParent, dir);
         }

         pZ = getParent(z);
      }

      makeBlack(rootNode);
   }

   void transplant(LinkType u, LinkType v, LinkType& nilParent) {
      LinkType p = getParent(u);

      if (p == NO_LINK) {
         rootNode = v;
      }else {
         setChild(p, Direction(u != getChild(p, LEFT)), v);
      }
      if (v == NO_LINK){
         nilParent = p;
      }else{
         setParent(v, p);
      }
   }

   void eraseFixup(LinkType x, const LinkType nilParent) {
      while (x != rootNode && isBlack(x)) {
         Direction dir;

         LinkType p = x == NO_LINK ? nilParent : getParent(x);
         dir = Direction(x == getChild(p, LEFT));
         LinkType w = getChild(p, dir);
         assert(w != NO_LINK);

         if (isRed(w)) {
            makeBlack(w);
            makeRed(p);
            rotate(p, opposite(dir));
            assert((x == NO_LINK && p == nilParent) ||
                   (x != NO_LINK && p == getParent(x)));
            w = getChild(p, dir);
            assert(w != NO_LINK);
         }

         if (isBlack(getChild(w, LEFT)) && isBlack(getChild(w, RIGHT))) {
            makeRed(w);
            x = p;
         } else {
            if (isBlack(getChild(w, dir))) {
               makeBlack(getChild(w, opposite(dir)));
               makeRed(w);
               rotate(w, dir);
               assert((x == NO_LINK && p == nilParent) ||
                      (x != NO_LINK && p == getParent(x)));
               w = getChild(p, dir);
            }
            setColor(w, getColor(p));
            makeBlack(p);
            makeBlack(getChild(w, dir));
            rotate(p, opposite(dir));
            x = rootNode;
         }
      }

      if (x != NO_LINK) makeBlack(x);
   }
 public:
   RbTree(LinkType &t_rootNode, LinkType& t_lowerNode) : rootNode{t_rootNode},
                                                         lowerNode{t_lowerNode}
                                                         {};

   ///Returns true if the RbTree is empty
   [[nodiscard]] bool empty() const { return rootNode == NO_LINK; }

   LinkType lower_bound(LinkType node) const {
      if(node == NO_LINK) return NO_LINK;
      while(true){
         LinkType child = getChild(node,LEFT);
         if(child == NO_LINK) return node;
         node = child;
      }
   }

   LinkType upper_bound(LinkType node) const{
      if(node == NO_LINK) return NO_LINK;
      while(true){
         LinkType child = getChild(node,RIGHT);
         if(child == NO_LINK) return node;
         node = child;
      }
   }

   LinkType lower_bound() const { return lowerNode; }

   LinkType upper_bound() const { return upper_bound(rootNode); }

   ///Gets the successor in the ordering
   LinkType successor(LinkType node) const {
      ///If the node has children, this is the smallest node which is greater (right)
      LinkType other_node = getChild(node, RIGHT);
      if (other_node != NO_LINK) return lower_bound(other_node);

      /// If not, we go back up the tree until our current node is a left child
      /// of some other node. This is the first node which is greater than the given node
      other_node = getParent(node);
      while (other_node != NO_LINK && node == getChild(other_node, RIGHT)) {
         node = other_node;
         other_node = getParent(node);
      }

      return other_node;
   }
   ///Gets predecessor in the ordering
   LinkType predecessor(LinkType node) const {
      ///If the node has children, this is the largest node which is smaller
      LinkType other_node = getChild(node, LEFT);
      if (other_node != NO_LINK) return upper_bound(other_node);

      /// If not, we go back up the tree until our current node is a right child
      /// of some other node. This is the first node which is greater than the given node
      other_node = getParent(node);
      while (other_node != NO_LINK && node == getChild(other_node, LEFT)) {
         node = other_node;
         other_node = getParent(node);
      }

      return other_node;
   }

   ///Inserts the node into the rb tree with the given parent
   void insert(LinkType node, LinkType parent) {
      // if we insert an element smaller than the parent, update the smallest element index
      if(lowerNode == parent &&
          (parent == RbTreeLinks<LinkType>::noLink() || getKey(node) < getKey(parent)) ){
         lowerNode = node;
      }
      setParent(node, parent);
      if (parent == NO_LINK){
         rootNode = node;
      }else{
         setChild(parent, Direction(getKey(parent) < getKey(node)), node);
      }

      setChild(node, LEFT, NO_LINK);
      setChild(node, RIGHT, NO_LINK);
      makeRed(node);
      insertFixup(node);
   }
   ///Inserts the new node into the rb-tree; first searches for an appropriate parent
   void insert(LinkType node) {
      LinkType y = NO_LINK;
      LinkType x = rootNode;
      while (x != NO_LINK) {
         y = x;
         x = getChild(y, Direction(getKey(x) < getKey(node)));
      }

      static_cast<T*>(this)->insert(node, y);
   }

   void erase(LinkType node) {
      if(node == lowerNode){
         lowerNode = successor(lowerNode);
      }
      LinkType nilParent = NO_LINK;
      LinkType y = node;
      bool yWasBlack = isBlack(y);
      LinkType x;

      if (getChild(node, LEFT) == NO_LINK) {
         x = getChild(node, RIGHT);
         transplant(node, x, nilParent);
      } else if (getChild(node, RIGHT) == NO_LINK) {
         x = getChild(node, LEFT);
         transplant(node, x, nilParent);
      } else {
         y = lower_bound(getChild(node, RIGHT));
         yWasBlack = isBlack(y);
         x = getChild(y, RIGHT);
         if (getParent(y) == node) {
            if (x == NO_LINK)
               nilParent = y;
            else
               setParent(x, y);
         } else {
            transplant(y, getChild(y, RIGHT), nilParent);
            LinkType zRight = getChild(node, RIGHT);
            setChild(y, RIGHT, zRight);
            setParent(zRight, y);
         }
         transplant(node, y, nilParent);
         LinkType zLeft = getChild(node, LEFT);
         setChild(y, LEFT, zLeft);
         setParent(zLeft, y);
         setColor(y, getColor(node));
      }

      if (yWasBlack) eraseFixup(x, nilParent);
   }
   ///Tries to find an element which has the given key in the given subtree
   ///Returns the closest?!? value to the key in the tree
   std::pair<LinkType,bool> find(const KeyType& key, LinkType treeRoot){
      LinkType y = NO_LINK;
      LinkType x = treeRoot;
      while(x != NO_LINK){
         int compare = 1 - (getKey(x) < key) + (key < getKey(x));
         switch (compare) {
         case 0:
            y = x;
            x = getChild(y, RIGHT);
            break;
         case 2:
            y = x;
            x = getChild(y, LEFT);
            break;
         default:
            return std::make_pair(x, true);
         }
      }
      return std::make_pair(y,false);
   }
   std::pair<LinkType,bool> find(const KeyType& key){
      return find(key,rootNode);
   }
};

} // namespace pcog
#endif // PCOG_RBTREE_HPP
