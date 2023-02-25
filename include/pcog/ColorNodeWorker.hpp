//
// Created by rolf on 25-2-23.
//

#ifndef PCOG_SRC_COLORNODEWORKER_HPP
#define PCOG_SRC_COLORNODEWORKER_HPP

#include "DenseGraph.hpp"
#include "BBNode.hpp"
#include "LPSolver.hpp"

namespace pcog {
/// This class is responsible for solving the branch-and-bound nodes e.g.
/// performing the pricing loop. It contains the information which is local to
/// each branch-and-bound node, such as the LP and graph representations.
class ColorNodeWorker {
 public:
   ColorNodeWorker() : m_focusNode{INVALID_BB_NODE}{};

   void processNode(BBNode& node);
 private:
   LPSolver m_lpSolver;
   node_id m_focusNode;
   DenseGraph m_focusGraph;
   DenseGraph m_completeFocusGraph;
};
} // namespace pcog

#endif // PCOG_SRC_COLORNODEWORKER_HPP
