//
// Created by rolf on 27-2-23.
//

#ifndef PCOG_INCLUDE_PCOG_SPARSESET_HPP
#define PCOG_INCLUDE_PCOG_SPARSESET_HPP

#include <algorithm>
#include "Definitions.hpp"
namespace pcog {
class DenseSet;

class SparseSet {
 public:
   SparseSet();
   explicit SparseSet(std::size_t capacity);
   ~SparseSet();
   SparseSet(const SparseSet &other);
   SparseSet &operator=(const SparseSet &other);
   SparseSet(SparseSet &&other) noexcept;
   SparseSet &operator=(SparseSet &&other) noexcept;

   void clear();

   void fromDenseSet(const DenseSet &set);
   void unsafe_add(Node node);
   [[nodiscard]] bool empty() const;

   std::size_t size();

   Node *begin();
   Node *end();
   [[nodiscard]] const Node *begin() const;
   [[nodiscard]] const Node *end() const;
   [[nodiscard]] const Node *cbegin() const;
   [[nodiscard]] const Node *cend() const;

   template <typename Func> void sortBy(Func f) {
      std::sort(begin(), end(), f);
   }

 private:
   Node *storage;
   std::size_t _size;
   std::size_t _capacity;
};
}
#endif // PCOG_INCLUDE_PCOG_SPARSESET_HPP
