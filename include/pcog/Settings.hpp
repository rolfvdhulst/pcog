//
// Created by rolf on 16-3-23.
//

#ifndef PCOG_INCLUDE_PCOG_SETTINGS_HPP
#define PCOG_INCLUDE_PCOG_SETTINGS_HPP

#include <limits>
static constexpr double NO_TIME_LIMIT = 1e20;

enum class BranchingStrategy : int {
   INTERSECTION_SIZE = 0,
   UNION_SIZE = 1,
   INTERSECTION_UNION_SIZE = 2,
   SYMMETRIC_DIFFERENCE_SIZE = 3,
   TRIANGLES_ADDED = 4,
   TRIANGLES_ADDED_SCALED = 5,
   HELDS_RULE = 6,
   RANDOMLY = 7,
   DUAL_MAX = 8,
   DUAL_MIN = 9,
   FRACTIONAL = 10,
   REMOVAL_SIZE = 11,
   MIN_REMOVAL_SIZE = 12,
};

enum class CandidateSelectionStrategy {
   VIOLATED_IN_BOTH = 0,
   VIOLATED_IN_SAME = 1,
   VIOLATED_IN_DIFFER = 2,
   VIOLATED_IN_ONE = 3,
   FIRST = 4
};

class Settings {
 public:
   Settings() : m_node_limit{std::numeric_limits<std::size_t>::max()},
                m_time_limit{NO_TIME_LIMIT},
                m_absgap_limit{0},
                m_relgap_limit{0.0},
                m_branchingStrategy{BranchingStrategy::INTERSECTION_UNION_SIZE},
                m_branchCandidateSelectionStrategy{CandidateSelectionStrategy::VIOLATED_IN_BOTH},
                m_rounding_tolerance{1e-8},
                m_diving_frequency{-1},
                m_diving_pricing_frequency{-1}
                {};

   ///Settings e.g. methods that will typically not affect performance (except for determining the starting/ending point)
   void setNodeLimit(std::size_t t_node_limit) { m_node_limit = t_node_limit; }
   void setTimeLimit(double t_time_limit) { m_time_limit = t_time_limit; }
   void setAbsGapLimit(std::size_t t_gap_limit) {
      m_absgap_limit = t_gap_limit;
   }
   void setRelGapLimit(double t_gap_limit) { m_relgap_limit = t_gap_limit; }
   [[nodiscard]] std::size_t nodeLimit() const { return m_node_limit; }
   [[nodiscard]] double timeLimit() const { return m_time_limit; }
   [[nodiscard]] std::size_t absGapLimit() const { return m_absgap_limit; }
   [[nodiscard]] double relGapLimit() const { return m_relgap_limit; }

   void setBranchingStrategy(BranchingStrategy t_strategy) {
      m_branchingStrategy = t_strategy;
   }
   [[nodiscard]] BranchingStrategy branchingStrategy() const {return m_branchingStrategy;}

   void setCandidateSelectionStrategy(CandidateSelectionStrategy t_strategy) {
      m_branchCandidateSelectionStrategy = t_strategy;
   }
   [[nodiscard]] CandidateSelectionStrategy candidateSelectionStrategy() const {
      return m_branchCandidateSelectionStrategy;
   }

   void setRoundingTolerance(double t_tolerance){
      assert(t_tolerance >= 0.0 && t_tolerance < 0.5);
      m_rounding_tolerance = t_tolerance;
   }
   [[nodiscard]] double roundingTolerance() const{
      return m_rounding_tolerance;
   }
   [[nodiscard]] bool isFeasOne(double t_value) const{
      return fabs(t_value-1.0) <= m_rounding_tolerance;
   }
   [[nodiscard]] bool isFeasZero(double t_value) const{
      return fabs(t_value) <= m_rounding_tolerance;
   }

   void setDivingFrequency(int t_frequency) { m_diving_frequency = t_frequency; }
   [[nodiscard]] int divingFrequency() const { return m_diving_frequency; }
   void setDivingPricingFrequency(int t_frequency) {
      m_diving_pricing_frequency = t_frequency;
   }
   [[nodiscard]] int divingPricingFrequency() const { return m_diving_pricing_frequency; }
 private:
   std::size_t m_node_limit; /// maximal number of nodes to process
   double m_time_limit; /// maximal time in seconds to run
   std::size_t m_absgap_limit; /// solving stops as soon as upperBound-lowerBound <= absolute_gap_limit is proven
   double m_relgap_limit; /// solving stops as soon as (upperBound-lowerBound)/lowerBound <= relative_gap_limit is proven

   BranchingStrategy m_branchingStrategy;
   CandidateSelectionStrategy m_branchCandidateSelectionStrategy;

   double m_rounding_tolerance; /// The rounding tolerance is used in various places to decide if a double value is integral or not
   /// Changing this value does not affect correctness of our algorithms since we use methods which protect against numerical errors

   /// Frequencies: -1 : never, 0 = only at root node, otherwise only at depths divisible by frequency
   int m_diving_frequency; /// Frequency for executing the diving heuristic
   int m_diving_pricing_frequency; /// Frequency for finding columns
};

#endif // PCOG_INCLUDE_PCOG_SETTINGS_HPP
