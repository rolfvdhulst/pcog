//
// Created by rolf on 8-4-23.
//

#ifndef PCOG_RBTREE_HPP
#define PCOG_RBTREE_HPP

#include <array>
#include <stdint.h>
namespace pcog {

// For now we only allow T to be a signed integral type:
// In the future one could also allow T to be a pointer type, but this would complicate the code
template <typename T>
concept RbTreeIndex = std::is_signed<T>::value;

template <typename T>
struct RbTreeTraits;

template <RbTreeIndex T> struct RbTreeLinks {

   std::array<T, 2> child;

   static constexpr T noLink() { return -1; }
   using ParentStorageType = std::make_unsigned<T>::type;
   ParentStorageType parentAndColor;

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
   bool isBlack() { return getColor() == 0; }
   bool isRed() { return getColor() == 1; }

   void setParent(T p) {
      assert(p >= 0);
      // keeping the color the same requires that we preserve the last bit
      parentAndColor = (parentAndColor & colorBitMask()) |
                       static_cast<ParentStorageType>(p + 1);
   }
   T getParent() const {
      return static_cast<T>(parentAndColor & ~colorBitMask()) - 1;
   }
};

template <typename T> class RbTree {
   enum Direction { LEFT = 0, RIGHT = 1 };
   using KeyType = typename RbTreeTraits<T>::KeyType;
   using LinkType = typename RbTreeTraits<T>::LinkType;

   LinkType& rootNode;

   static constexpr LinkType NO_LINK = RbTreeLinks<T>::noLink();

 public:


};

}
#endif // PCOG_RBTREE_HPP
