//
// Created by rolf on 27-2-23.
//

#include <cassert>
#include "pcog/utilities/SparseSet.hpp"
#include "pcog/utilities/DenseSet.hpp"
namespace pcog {
SparseSet::SparseSet() : storage{nullptr}, _size{0}, _capacity{0} {}

SparseSet::SparseSet(std::size_t capacity)
    : storage{new Node[capacity]}, _size{0}, _capacity{capacity} {}

SparseSet::~SparseSet() { delete[] storage; }

SparseSet::SparseSet(const SparseSet &other)
    : storage{new Node[other._capacity]}, _size{other._size},
      _capacity{other._capacity} {
   for (std::size_t i = 0; i < other._size; ++i) {
      storage[i] = other.storage[i];
   }
}

SparseSet &SparseSet::operator=(const SparseSet &other) {
   if (this == &other) {
      return *this;
   }
   if (_capacity != other._capacity) {
      delete[] storage;
      storage = new Node[other._capacity];
   }
   for (std::size_t i = 0; i < other._size; ++i) {
      storage[i] = other.storage[i];
   }
   _size = other._size;
   _capacity = other._capacity;
   return *this;
}

SparseSet::SparseSet(SparseSet &&other) noexcept
    : storage{other.storage}, _size{other._size}, _capacity{other._capacity} {
   other.storage = nullptr;
   other._size = 0;
   other._capacity = 0;
}

SparseSet &SparseSet::operator=(SparseSet &&other) noexcept {
   if (this == &other) {
      return *this;
   }
   delete[] storage;
   storage = other.storage;
   _size = other._size;
   _capacity = other._capacity;

   other.storage = nullptr;
   other._size = 0;
   other._capacity = 0;
   return *this;
}

void SparseSet::clear() { _size = 0; }

Node *SparseSet::begin() { return &storage[0]; }

Node *SparseSet::end() { return &storage[_size]; }

const Node *SparseSet::begin() const { return &storage[0]; }

const Node *SparseSet::end() const { return &storage[_size]; }

const Node *SparseSet::cbegin() const { return &storage[0]; }

const Node *SparseSet::cend() const { return &storage[_size]; }

std::size_t SparseSet::size() { return _size; }

bool SparseSet::empty() const { return _size == 0; }

void SparseSet::unsafe_add(Node node) {
   storage[_size] = node;
   _size++;
   assert(_size <= _capacity);
}

void SparseSet::fromDenseSet(const DenseSet &set) {
   assert(set.capacity() == _capacity); // probably a bug if not
   std::size_t index = 0;
   for (const auto &node : set) {
      storage[index] = node;
      index++;
   }
   _size = index;
   assert(_size <= _capacity);
}
}