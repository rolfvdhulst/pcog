//
// Created by rolf on 5-7-23.
//

#ifndef PCOG_SRC_REDUCTION_REDUCTIONQUEUE_HPP
#define PCOG_SRC_REDUCTION_REDUCTIONQUEUE_HPP
#include "pcog/utilities/DenseSet.hpp"
namespace pcog {
class ReductionVertexQueue {
 public:
   ReductionVertexQueue(std::size_t numNodes, bool allInitiallyPresent);
   [[nodiscard]] std::size_t size() const;
   [[nodiscard]] bool empty() const;
   void push(Node node);
   Node pop();
 private:
   DenseSet inQueue;
   std::vector<Node> queueArray;
   std::size_t head;
   std::size_t tail;
   std::size_t queueSize;

};
}
#endif // PCOG_SRC_REDUCTION_REDUCTIONQUEUE_HPP
