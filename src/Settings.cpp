//
// Created by rolf on 16-3-23.
//

#include "pcog/Settings.hpp"

Settings::Settings()
    : m_node_limit{std::numeric_limits<std::size_t>::max()},
      m_time_limit{NO_TIME_LIMIT}, m_absgap_limit{0}, m_relgap_limit{0.0},
      m_branchingStrategy{BranchingStrategy::TRIANGLES_ADDED},
      m_branchCandidateSelectionStrategy{
          CandidateSelectionStrategy::VIOLATED_IN_BOTH},
      m_nodeSelectionStrategy{NodeSelectionStrategy::DFS_BOUND},
      m_nodeChildSelectionStrategy{NodeChildSelectionStrategy::PREFER_DIFFER},
      m_dfsRestartFrequency{10}, m_rounding_tolerance{1e-8},
      m_diving_frequency{0}, m_diving_pricing_frequency{-1},
      m_nodeDisplayFrequency{25}, m_numMaxThreads{16}
      {};