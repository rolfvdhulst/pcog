//
// Created by rolf on 5-7-23.
//

#include "pcog/reduction/ReductionStack.hpp"
namespace pcog {
void ReductionStack::push(const Reduction& reduction) {
   reductions.push_back(reduction);
}
}