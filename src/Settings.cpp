//
// Created by rolf on 16-3-23.
//

#include "pcog/Settings.hpp"

Settings::Settings()
    : m_node_limit{std::numeric_limits<std::size_t>::max()},
      m_time_limit{NO_TIME_LIMIT}, m_absgap_limit{0}, m_relgap_limit{0.0},
      m_branchingStrategy{BranchingStrategy::INTERSECTION_UNION_SIZE},
      m_branchCandidateSelectionStrategy{
          CandidateSelectionStrategy::VIOLATED_IN_BOTH},
      m_nodeSelectionStrategy{NodeSelectionStrategy::DFS_BOUND},
      m_nodeChildSelectionStrategy{NodeChildSelectionStrategy::PREFER_DIFFER},
      m_dfsRestartFrequency{10}, m_rounding_tolerance{1e-8},
      m_diving_frequency{5}, m_diving_pricing_frequency{0},
      m_nodeDisplayFrequency{50}, m_numMaxThreads{32}
      {};