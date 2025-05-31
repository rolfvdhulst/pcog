//
// Created by rolf on 16-3-23.
//

#ifndef PCOG_INCLUDE_PCOG_SETTINGS_HPP
#define PCOG_INCLUDE_PCOG_SETTINGS_HPP

#include <limits>
#include <cassert>
#include <cmath>

static constexpr double NO_TIME_LIMIT = 1e20;
static constexpr std::size_t NO_NODE_LIMIT = std::numeric_limits<std::size_t>::max();
static constexpr double TIMEOUT_SAFETY = 1e-3;
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

enum class CandidateSelectionStrategy : int {
   VIOLATED_IN_BOTH = 0,
   VIOLATED_IN_SAME = 1,
   VIOLATED_IN_DIFFER = 2,
   VIOLATED_IN_ONE = 3,
   FIRST = 4
};

enum class NodeSelectionStrategy : int{
   BOUND = 0,     /// Choose the b&b node with smallest lower bound.
   DFS_BOUND = 1, /// Pick child nodes of the current node until a child node has a
              /// lower bound greater than the current global lower bound.
              /// In this case, it backtracks to the node with smallest lower bound.
              /// This method is useful for hard instances to quickly improve the lower bound
   DFS_RESTART = 2 /// Pick a child node of the current node if possible.
               /// Restart every x nodes where x is equal to restartDFSFrequency,
               /// choosing the node with smallest bound.
};

///When picking a child node, decide whether to pick the left/right node for the current relaxation.
///Typically, one wants to pick the weakest node, which is often DIFFER, as this leads to less time wasted reloading the LP.
enum class NodeChildSelectionStrategy : int{
   PREFER_SAME = 0,
   PREFER_DIFFER = 1,
   RANDOMLY = 2
};

enum class PricingAlgorithmStrategy : int{
   EXACT = 0, ///Only use the exact combinatorial algorithm. This is only recommended for graphs where the pricing problem is very easy.
   HEURISTIC_LOW = 1, ///Call the greedy search heuristic.
   HEURISTIC_HIGH = 2, //Call the extensive greedy search heuristic with many different orderings. Recommended for graphs where the pricing problem becomes very hard.
};

///Control whether the dual weights for the pricing problem are altered.
///This may make the pricing problem easier and can avoid the ´tailing off' effect and reduce the number of pricing iterations.
///The downside is that not completly solving the pricing loop typically leads to worse branching decisions
enum class DualWeightReductionStrategy : int {
   NONE = 0,          ///Do not alter the dual weights for the pricing problem
   NEIGHBOURHOOD = 1, ///Alter the dual weights for the pricing problem by reducing the weights around the largest-value node.
   UNIFORM = 2        ///Alter the dual weights for the pricing problem
};

class Settings {
 public:
   Settings();

   void setNodeDisplayFrequency(std::size_t frequency) {m_nodeDisplayFrequency = frequency;}
   [[nodiscard]] std::size_t nodeDisplayFrequency() const {return m_nodeDisplayFrequency;}

   ///Settings e.g. methods that will typically not affect performance (except for determining the starting/ending point)
   void setNodeLimit(std::size_t t_node_limit) { m_node_limit = t_node_limit; }
   void setTimeLimit(double t_time_limit) { m_time_limit = t_time_limit; }
   [[nodiscard]] std::size_t nodeLimit() const { return m_node_limit; }
   [[nodiscard]] double timeLimit() const { return m_time_limit; }

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

   void setNumThreads(std::size_t t_numThreads){
      m_numMaxThreads = t_numThreads;
   }
   [[nodiscard]] std::size_t numThreads() const { return m_numMaxThreads;}

   void setNumInitialTabuIterations(std::size_t numIters) { m_numInitialTabuIterations = numIters;}
   [[nodiscard]] std::size_t numInitialTabuIterations() const {return m_numInitialTabuIterations;}

   [[nodiscard]] PricingAlgorithmStrategy getPricingStrategy() const { return m_pricingAlgorithmStrategy;}
   void setPricingAlgorithmStrategy(PricingAlgorithmStrategy t_strategy) { m_pricingAlgorithmStrategy = t_strategy;}

   void setDualWeightReductionStrategy(DualWeightReductionStrategy t_strategy) { m_dualWeightReductionStrategy = t_strategy; }
   [[nodiscard]] DualWeightReductionStrategy dualWeightReductionStrategy() const {return m_dualWeightReductionStrategy;}

   void setSolveNodesCompletely(bool t_solveCompletely) { m_solveNodesCompletely = t_solveCompletely; }
   [[nodiscard]] bool solveNodesCompletely() const { return m_solveNodesCompletely;}
 private:
   std::size_t m_node_limit; /// maximal number of nodes to process
   double m_time_limit; /// maximal time in seconds to run

   BranchingStrategy m_branchingStrategy;
   CandidateSelectionStrategy m_branchCandidateSelectionStrategy;
   PricingAlgorithmStrategy m_pricingAlgorithmStrategy;
   DualWeightReductionStrategy m_dualWeightReductionStrategy;
   bool m_solveNodesCompletely;

   NodeSelectionStrategy m_nodeSelectionStrategy;
   NodeChildSelectionStrategy m_nodeChildSelectionStrategy;
   std::size_t m_dfsRestartFrequency;

   double m_rounding_tolerance; /// The rounding tolerance is used in various places to decide if a double value is integral or not
   /// Changing this value does not affect correctness of our algorithms since we use methods which protect against numerical errors

   /// Frequencies: -1 : never, 0 = only at root node, otherwise only at depths divisible by frequency
   int m_diving_frequency; /// Frequency for executing the diving heuristic
   int m_diving_pricing_frequency; /// Frequency for finding columns

   std::size_t m_nodeDisplayFrequency;

   std::size_t m_numMaxThreads; /// Maximum # of threads to run the program on.
                                /// We use less if the system that we are running has less threads

   std::size_t m_numInitialTabuIterations;
};

#endif // PCOG_INCLUDE_PCOG_SETTINGS_HPP
