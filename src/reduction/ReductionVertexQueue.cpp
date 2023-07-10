//
// Created by rolf on 5-7-23.
//

#include "pcog/reduction/ReductionVertexQueue.hpp"

namespace pcog {

std::size_t ReductionVertexQueue::size() const {
   return queueSize + lowPriorityQueueSize;
}
bool ReductionVertexQueue::empty() const {
   return queueSize == 0 && lowPriorityQueueSize == 0;
}

void ReductionVertexQueue::push(Node node, bool isLowPriority) {
   if (!inQueue.contains(node)) {
      inQueue.add(node);
      if (isLowPriority) {
         assert(lowPriorityQueueSize < lowPriorityQueueArray.size());
         lowPriorityQueueArray[lowPriorityTail] = node;

         ++lowPriorityQueueSize;
         ++lowPriorityTail;
         if (lowPriorityTail == lowPriorityQueueArray.size()) {
            lowPriorityTail = 0;
         }
      } else {
         assert(queueSize < queueArray.size());
         queueArray[tail] = node;

         ++queueSize;
         ++tail;
         if (tail == queueArray.size()) {
            tail = 0;
         }
      }
   }
}
Node ReductionVertexQueue::pop() {
   assert(!empty());
   // First try high priority nodes
   if (queueSize != 0) {
      Node node = queueArray[head];
      ++head;
      if (head == queueArray.size()) {
         head = 0;
      }
      --queueSize;

      inQueue.remove(node);
      return node;
   }
   Node node = lowPriorityQueueArray[lowPriorityHead];
   ++lowPriorityHead;
   if (lowPriorityHead == lowPriorityQueueArray.size()) {
      lowPriorityHead = 0;
   }
   --lowPriorityQueueSize;

   inQueue.remove(node);
   return node;
}
ReductionVertexQueue::ReductionVertexQueue(const DenseSet &lowPriorityNodes)
    : inQueue(lowPriorityNodes.capacity(),true),
      queueArray(lowPriorityNodes.capacity(), 0), head{0}, tail{0},
      queueSize{0}, lowPriorityQueueArray(lowPriorityNodes.capacity(), 0),
      lowPriorityHead{0}, lowPriorityTail{0}, lowPriorityQueueSize{0} {

   for(std::size_t i = 0; i < lowPriorityNodes.capacity(); ++i){
      if(lowPriorityNodes.contains(i)){
         lowPriorityQueueArray[lowPriorityTail] = i;
         ++lowPriorityTail;
         ++lowPriorityQueueSize;
      }else{
         queueArray[tail] = i;
         ++tail;
         ++queueSize;
      }
   }

}
} // namespace pcog