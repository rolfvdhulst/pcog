//
// Created by rolf on 16-3-23.
//

#ifndef PCOG_INCLUDE_PCOG_SETTINGS_HPP
#define PCOG_INCLUDE_PCOG_SETTINGS_HPP

#include <limits>
#include <cassert>
#include <cmath>

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
   SMALL_DIFFERENCE = 13
};

enum class CandidateSelectionStrategy {
   VIOLATED_IN_BOTH = 0,
   VIOLATED_IN_SAME = 1,
   VIOLATED_IN_DIFFER = 2,
   VIOLATED_IN_ONE = 3,
   FIRST = 4
};

enum class NodeSelectionStrategy {
   BOUND,     /// Choose the b&b node with smallest lower bound.
   DFS_BOUND, /// Pick child nodes of the current node until a child node has a
              /// lower bound greater than the current global lower bound.
              /// In this case, it backtracks to the node with smallest lower bound.
              /// This method is useful for hard instances to quickly improve the lower bound
   DFS_RESTART /// Pick a child node of the current node if possible.
               /// Restart every x nodes where x is equal to restartDFSFrequency,
               /// choosing the node with smallest bound.
};

///When picking a child node
enum class NodeChildSelectionStrategy{
   PREFER_SAME,
   PREFER_DIFFER,
   RANDOMLY
};
class Settings {
 public:
   Settings();

   void setNodeDisplayFrequency(std::size_t frequency) {m_nodeDisplayFrequency = frequency;}
   [[nodiscard]] std::size_t nodeDisplayFrequency() const {return m_nodeDisplayFrequency;}

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

   void setNodeSelectionStrategy(NodeSelectionStrategy t_strategy){
      m_nodeSelectionStrategy = t_strategy;
   }
   [[nodiscard]] NodeSelectionStrategy nodeSelectionStrategy() const {return m_nodeSelectionStrategy;}

   void setDfsRestartFrequency(std::size_t t_frequency){
      m_dfsRestartFrequency = t_frequency;
   }
   [[nodiscard]] std::size_t dfsRestartFrequency() const {return m_dfsRestartFrequency;}

   void setNodeChildSelectionStrategy(NodeChildSelectionStrategy t_strategy){
      m_nodeChildSelectionStrategy = t_strategy;
   }
   [[nodiscard]] NodeChildSelectionStrategy nodeChildSelectionStrategy() const {return m_nodeChildSelectionStrategy;}
 private:
   std::size_t m_node_limit; /// maximal number of nodes to process
   double m_time_limit; /// maximal time in seconds to run
   std::size_t m_absgap_limit; /// solving stops as soon as upperBound-lowerBound <= absolute_gap_limit is proven
   double m_relgap_limit; /// solving stops as soon as (upperBound-lowerBound)/lowerBound <= relative_gap_limit is proven

   BranchingStrategy m_branchingStrategy;
   CandidateSelectionStrategy m_branchCandidateSelectionStrategy;


   NodeSelectionStrategy m_nodeSelectionStrategy;
   NodeChildSelectionStrategy m_nodeChildSelectionStrategy;
   std::size_t m_dfsRestartFrequency;

   double m_rounding_tolerance; /// The rounding tolerance is used in various places to decide if a double value is integral or not
   /// Changing this value does not affect correctness of our algorithms since we use methods which protect against numerical errors

   /// Frequencies: -1 : never, 0 = only at root node, otherwise only at depths divisible by frequency
   int m_diving_frequency; /// Frequency for executing the diving heuristic
   int m_diving_pricing_frequency; /// Frequency for finding columns

   std::size_t m_nodeDisplayFrequency;
};

#endif // PCOG_INCLUDE_PCOG_SETTINGS_HPP
