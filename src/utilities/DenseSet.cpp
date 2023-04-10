//
// Created by rolf on 22-11-22.
//

#include "pcog/utilities/DenseSet.hpp"
#include <iostream>
namespace pcog {
degree_type DenseSet::size() const noexcept { return m_bitset.count(); }

degree_type DenseSet::capacity() const noexcept { return m_bitset.size(); }

DenseSet::DenseSet(degree_type t_capacity) : m_bitset(t_capacity) {}

bool DenseSet::empty() const noexcept { return !m_bitset.any(); }

bool DenseSet::full() const { return m_bitset.all(); }

bool DenseSet::any() const { return m_bitset.any(); }

bool DenseSet::contains(Element t_node) const {
   assert(t_node < m_bitset.size());
   return m_bitset[t_node];
}

DenseSet &DenseSet::inplaceIntersection(const DenseSet &t_other) {
   m_bitset &= t_other.m_bitset;
   return *this;
}

DenseSet &DenseSet::inplaceUnion(const DenseSet &t_other) {
   m_bitset |= t_other.m_bitset;
   return *this;
}

DenseSet &DenseSet::inplaceSymmetricDifference(const DenseSet &t_other) {
   m_bitset ^= t_other.m_bitset;
   return *this;
}

DenseSet &DenseSet::inplaceDifference(const DenseSet &t_other) {
   m_bitset -= t_other.m_bitset;
   return *this;
}

DenseSet &DenseSet::complement() {
   m_bitset.flip();
   return *this;
}

DenseSet DenseSet::intersection(const DenseSet &t_other) const {
   return DenseSet(m_bitset & t_other.m_bitset);
}

DenseSet DenseSet::setUnion(const DenseSet &t_other) const {
   return DenseSet(m_bitset | t_other.m_bitset);
}

DenseSet DenseSet::symmetricDifference(const DenseSet &t_other) const {
   return DenseSet(m_bitset ^ t_other.m_bitset);
}

DenseSet DenseSet::difference(const DenseSet &t_other) const {
   return DenseSet(m_bitset - t_other.m_bitset);
}

DenseSet::DenseSet(boost::dynamic_bitset<> t_other)
    : m_bitset(std::move(t_other)) {}

DenseSet::Element DenseSet::first() const { return m_bitset.find_first(); }

DenseSet::Element DenseSet::find_next(Element node) const {
   return m_bitset.find_next(node);
}

DenseSet &DenseSet::add(Element node) {
   m_bitset.set(node, true);
   return *this;
}

DenseSet &DenseSet::flip(Element node) {
   m_bitset.flip(node);
   return *this;
}

DenseSet &DenseSet::remove(Element node) {
   m_bitset.set(node, false);
   return *this;
}

void DenseSet::clear() { m_bitset.reset(); }

DenseSet::ConstIterator DenseSet::begin() const {
   return {*this, m_bitset.find_first()};
}
DenseSet::ConstIterator DenseSet::end() const { return {*this, INVALID_NODE}; }

bool DenseSet::operator==(const DenseSet &t_other) const {
   return m_bitset == t_other.m_bitset;
}

bool DenseSet::operator!=(const DenseSet &t_other) const {
   return m_bitset != t_other.m_bitset;
}

bool DenseSet::hasAsSubset(const DenseSet &t_other) const {
   return t_other.m_bitset.is_subset_of(m_bitset);
}

bool DenseSet::isSubsetOf(const DenseSet &t_other) const {
   return m_bitset.is_subset_of(t_other.m_bitset);
}

bool DenseSet::hasAsProperSubset(const DenseSet &t_other) const {
   return t_other.m_bitset.is_proper_subset_of(m_bitset);
}

bool DenseSet::isProperSubsetOf(const DenseSet &t_other) const {
   return m_bitset.is_proper_subset_of(t_other.m_bitset);
}

DenseSet::DenseSet(degree_type capacity, bool setFull) : m_bitset(capacity) {
   m_bitset.set(0, capacity, setFull);
}

DenseSet &DenseSet::set(Element node, bool value) {
   m_bitset.set(node, value);
   return *this;
}

DenseSet &DenseSet::setRange(Element first, Element last) {
   assert(last >= first);
   m_bitset.set(first, last - first + 1, true);
   return *this;
}

DenseSet &DenseSet::setAll() { return setRange(0, capacity() - 1); }

void DenseSet::print() const {
   for (const auto &elem : *this) {
      std::cout << elem << ",";
   }
   std::cout << "\n";
}

DenseSet &DenseSet::extend(bool value) {
   m_bitset.push_back(value);
   return *this;
}

DenseSet::ConstIterator::ConstIterator(const DenseSet &t_set,
                                       Element t_position)
    : m_node(t_position), m_denseSet(t_set) {}

DenseSet::Element DenseSet::ConstIterator::operator*() const { return m_node; }

bool DenseSet::ConstIterator::operator==(
    const DenseSet::ConstIterator &t_other) const {
   assert(this->m_denseSet ==
          t_other.m_denseSet); // we should be iterating over the same set;
                               // other things do not make sense
   return m_node == t_other.m_node;
}

bool DenseSet::ConstIterator::operator!=(
    const DenseSet::ConstIterator &t_other) const {
   return !(*this == t_other);
}
} // namespace pcog