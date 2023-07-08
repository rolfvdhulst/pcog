//
// Created by rolf on 5-7-23.
//

#include "pcog/reduction/ReductionVertexQueue.hpp"

namespace pcog {

std::size_t ReductionVertexQueue::size() const{
   return queueSize;
}
bool ReductionVertexQueue::empty() const{
   return queueSize == 0;
}

void ReductionVertexQueue::push(Node node){
   if(!inQueue.contains(node)){
      inQueue.add(node);
      assert(queueSize < queueArray.size());
      queueArray[tail] = node;

      ++queueSize;
      ++tail;
      if(tail == queueArray.size()){
         tail = 0;
      }

   }
}
Node ReductionVertexQueue::pop(){
   assert(!empty());
   assert(inQueue.contains(queueArray[head]));
   Node node = queueArray[head];
   ++head;
   if(head == queueArray.size()){
      head = 0;
   }
   --queueSize;

   inQueue.remove(node);
   return node;
}
ReductionVertexQueue::ReductionVertexQueue(std::size_t numNodes, bool allInitiallyPresent) :inQueue(numNodes,allInitiallyPresent), queueArray(numNodes,0){

   //items in queue are [head,tail)
   //If head == queueArray.size() head is reset to zero
   if(allInitiallyPresent){
      for (std::size_t i = 0; i < numNodes; ++i) {
         queueArray[i] = i;
      }
      head = 0;
      tail = 0;
      queueSize = numNodes;
   }else{
      head = 0;
      tail = 0;
      queueSize = 0;
   }

}
}