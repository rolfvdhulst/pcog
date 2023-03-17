//
// Created by rolf on 16-3-23.
//

#ifndef PCOG_INCLUDE_PCOG_SETTINGS_HPP
#define PCOG_INCLUDE_PCOG_SETTINGS_HPP

#include <limits>
static constexpr double NO_TIME_LIMIT = 1e20;
class Settings {
 public:
   Settings() : m_node_limit{std::numeric_limits<std::size_t>::max()},
                m_time_limit{1e20},
                m_absgap_limit{0},
                m_relgap_limit{0.0}{};

   void setNodeLimit(std::size_t t_node_limit) { m_node_limit = t_node_limit;}
   void setTimeLimit(double t_time_limit){m_time_limit = t_time_limit;}
   void setAbsGapLimit(std::size_t t_gap_limit) {m_absgap_limit = t_gap_limit;}
   void setRelGapLimit(double t_gap_limit ){m_relgap_limit = t_gap_limit;}
   [[nodiscard]] std::size_t nodeLimit() const {return m_node_limit;}
   [[nodiscard]] double timeLimit() const {return m_time_limit;}
   [[nodiscard]] std::size_t absGapLimit() const { return m_absgap_limit;}
   [[nodiscard]] double relGapLimit() const {return m_relgap_limit;}
 private:
   std::size_t m_node_limit; //maximal number of nodes to process
   double m_time_limit; // maximal time in seconds to run
   std::size_t m_absgap_limit; //solving stops as soon as upperBound-lowerBound <= m_absgap_limit is proven
   double m_relgap_limit; //solving stops as soon as (upperBound-lowerBound)/lowerBound <= m_relgap_limit
};

#endif // PCOG_INCLUDE_PCOG_SETTINGS_HPP
