//
// Created by rolf on 22-11-22.
//

#ifndef PCOG_SRC_DENSESET_HPP
#define PCOG_SRC_DENSESET_HPP

#include "Definitions.hpp"
#include <boost/dynamic_bitset.hpp>

class DenseSet {
 private:
   using BlockType = std::size_t;
   using Element = std::size_t;
   boost::dynamic_bitset<BlockType> m_bitset;

 public:
   constexpr static auto INVALID_ELEMENT = boost::dynamic_bitset<>::npos;
   static_assert(INVALID_ELEMENT == INVALID_NODE);

   ///Default constructed sets initially are empty and have zero capacity.
   DenseSet() = default; // Rule of zero
   /// Constructs an empty denseSet with the given capacity
   /// \param t_capacity The capacity to set this set to.
   explicit DenseSet(degree_type t_capacity);

   /// Constructs a dense set from a dynamic_bitset, the underlying type
   explicit DenseSet(boost::dynamic_bitset<> t_other);

   /// Constructs a dense set with the given capacity and sets all values if t_setFull is true.
   /// \param t_capacity
   /// \param t_setFull
   explicit DenseSet(degree_type t_capacity, bool t_setFull);

   /// Returns the number of elements in the current dense set.
   /// Note this is not equal to the buffer size, which can be found by querying
   /// the capacity. \return the number of elements currently in the set
   [[nodiscard]] degree_type size() const noexcept;

   /// The maximum number of elements which fit in this dense set
   /// \return The maximum number of elements which fit in this dense set
   [[nodiscard]] degree_type capacity() const noexcept;

   /// Checks if the set contains any element.
   /// \return Returns true if this set is non-empty, false otherwise
   [[nodiscard]] bool any() const;

   ///
   /// \return True if the set is empty, false other wise
   [[nodiscard]] bool empty() const noexcept;

   /// Returns true if all elements are set, e.g. all elements until capacity()
   /// are in the set. \return True if all elements are set, false otherwise
   [[nodiscard]] bool full() const;
   ///
   /// \return Returns true if the element is contained in the set, false
   /// otherwise.
   [[nodiscard]] bool contains(Element node) const;

   /// Flips all the bits in this set, setting all contained elements to false
   /// and all elements which are not contained to true. \return a reference to
   /// this
   DenseSet &complement();

   /// Intersects this set with an other set and stores the result of the
   /// intersection in this. \param other The other set to intersect this set
   /// with \return A reference to this
   DenseSet &inplaceIntersection(const DenseSet &other);
   /// Takes the union of this set and another set and stores the result in this
   /// set. \param other The other set to take the union with \return A
   /// reference to this
   DenseSet &inplaceUnion(const DenseSet &other);

   /// Takes the symmetric difference of this set and another set and stores the
   /// result in this set. \param other The other set to take the symmetric
   /// difference with \return A reference to this
   DenseSet &inplaceSymmetricDifference(const DenseSet &other);

   /// Computes *this - other and stores the result in this
   /// \param other The other set to be substracted from this
   /// \return A reference to this
   DenseSet &inplaceDifference(const DenseSet &other);


   /// Computes the intersection of this with other without modifying either.
   /// \param other the other set
   /// \return A DenseSet which is the intersection of this and other
   [[nodiscard]] DenseSet intersection(const DenseSet &other) const;
   /// Computes the union of this with other without modifying either.
   /// \param other the other set
   /// \return A DenseSet which is the union of this and other
   [[nodiscard]] DenseSet setUnion(const DenseSet &other) const;
   /// Computes the symmetric difference of this with other without modifying either.
   /// \param other the other set
   /// \return A DenseSet which is the symmetric difference of this and other
   [[nodiscard]] DenseSet symmetricDifference(const DenseSet &other) const;
   /// Computes this - other without modifying either set.
   /// \param other the other set
   /// \return A new DenseSet set such that set = this - other
   [[nodiscard]] DenseSet difference(const DenseSet &other) const;

   /// Clears all the elements from this set.
   void clear();

   ///Prints the contents of this set, useful for debugging
   void print() const;

   /// Checks if two sets are equal. Note that this requires both sets to have the same capacity to work correctly.
   /// \return true if the sets are equal, false otherwise
   bool operator==(const DenseSet &other) const;
   /// Checks if two sets are not equal. Note that this requires both sets to have the same capacity to work correctly.
   /// \return true if the sets are not equal, false otherwise
   bool operator!=(const DenseSet &other) const;
   /// Checks if the given set is a subset of this.
   /// \param other the other set.
   /// \return true if other is a subset of this, false otherwise
   [[nodiscard]] bool hasAsSubset(const DenseSet &other) const;
   /// Checks if this is a subset of other
   /// \param other the other set.
   /// \return true if this is a subset of other, false otherwise
   [[nodiscard]] bool isSubsetOf(const DenseSet &other) const;
   /// Checks if the given set is a proper subset of this.
   /// \param other the other set.
   /// \return true if other is a proper subset of this, false otherwise
   [[nodiscard]] bool hasAsProperSubset(const DenseSet &other) const;
   /// Checks if this is a proper subset of other
   /// \param other the other set.
   /// \return true if this is a proper subset of other, false otherwise
   [[nodiscard]] bool isProperSubsetOf(const DenseSet &other) const;

   //Modifying functions

   /// Increases the capacity of this set by 1 and sets the newly added bit to true or false
   /// \param value true if the new bit is set.
   /// \return a reference to this
   DenseSet &extend(bool value);

   /// Sets the given element to true or false
   /// \param node the given element
   /// \param value the value the element is set to
   /// \return A reference to this
   DenseSet &set(Element node, bool value);

   /// Adds the element to this. Does nothing if the element is already in the set.
   /// \param node The element to be added
   /// \return A reference to this
   DenseSet &add(Element node);
   /// Flips the given element, removing it if it is in the set and adding it if it is not.
   /// \param node The element to be flipped
   /// \return A reference to this
   DenseSet &flip(Element node);
   /// Removes the given node from the set
   /// \param node The given node
   /// \return A reference to this
   DenseSet &remove(Element node);
   /// Adds all elements in the range
   /// \param first The first element of the range to be added
   /// \param last The last element of the range to be added
   /// \return A reference to this
   DenseSet &setRange(Element first, Element last);

   /// Sets all bits in the set to true
   /// \return A reference to this
   DenseSet &setAll();

   ///
   /// \return The lowest element which is set in this dense set. If none is
   /// set, INVALID_ELEMENT is returned.
   [[nodiscard]] Element first() const;

   ///
   /// \return The lowest element which is set after node in this dense
   /// set. If there is no element after node, INVALID_ELEMENT is
   /// returned
   [[nodiscard]] Element find_next(Element node) const;

   //TODO: document iterator type?
   class ConstIterator {
    public:
      ConstIterator(const DenseSet &t_set, Element t_position);
      [[nodiscard]] Element operator*() const;
      bool operator==(const ConstIterator &t_other) const;
      bool operator!=(const ConstIterator &t_other) const;
      ConstIterator &operator++() {
         m_node = m_denseSet.find_next(m_node);
         return *this;
      }

    private:
      Node m_node;
      const DenseSet &m_denseSet;
   };

   /// Returns an iterator to iterate over the set.
   /// \return The start iterator for the set.
   [[nodiscard]] ConstIterator begin() const;

   /// Returns an iterator to iterate over the set.
   /// \return The end iterator for the set.
   [[nodiscard]] ConstIterator end() const;
};

#endif // PCOG_SRC_DENSESET_HPP
