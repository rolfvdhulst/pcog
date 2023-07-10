//
// Created by rolf on 5-7-23.
//

#ifndef PCOG_SRC_REDUCTION_REDUCTIONQUEUE_HPP
#define PCOG_SRC_REDUCTION_REDUCTIONQUEUE_HPP
#include "pcog/utilities/DenseSet.hpp"
namespace pcog {
class ReductionVertexQueue {
 public:
   ReductionVertexQueue(const DenseSet& lowPriorityNodes);
   [[nodiscard]] std::size_t size() const;
   [[nodiscard]] bool empty() const;
   void push(Node node, bool isLowPriority);
   Node pop();
 private:
   DenseSet inQueue;
   std::vector<Node> queueArray;
   std::size_t head;
   std::size_t tail;
   std::size_t queueSize;

   std::vector<Node> lowPriorityQueueArray;
   std::size_t lowPriorityHead;
   std::size_t lowPriorityTail;
   std::size_t lowPriorityQueueSize;

};
}
#endif // PCOG_SRC_REDUCTION_REDUCTIONQUEUE_HPP
