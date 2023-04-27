#pragma once

#include "common.hpp"
#include "formula.hpp"
#include "graph.hpp"

namespace dpve{
class JoinNode { // abstract
public:
  static Int nodeCount;
  static Int terminalCount;
  static Set<Int> nonterminalIndices;

  static Int backupNodeCount;
  static Int backupTerminalCount;
  static Set<Int> backupNonterminalIndices;

  // static Cnf cnf; // this field must be set exactly once before any JoinNode object is constructed

  Int nodeIndex = MIN_INT; // 0-indexing (equal to clauseIndex for JoinTerminal)
  vector<JoinNode*> children; // empty for JoinTerminal
  Set<Int> projectionVars; // empty for JoinTerminal
  Set<Int> preProjectionVars; // set by constructor

  void* dd = 0; // for sampling. void* type so as to avoid circular dependency with class Dd definition in dmc.hh

  static void resetStaticFields(); // backs up and re-initializes static fields
  static void restoreStaticFields(); // from backup

  virtual Int getWidth(const Assignment& assignment = Assignment()) const = 0; // of subtree

  virtual void updateVarSizes(
    Map<Int, size_t>& varSizes, // var x |-> size of biggest node containing x
    const Cnf& cnf
  ) const = 0;

  Set<Int> getPostProjectionVars() const;
  Int chooseClusterIndex(
    Int clusterIndex, // of this node
    const vector<Set<Int>>& projectableVarSets, // Z_1, ..., Z_m
    string clusteringHeuristic
  ); // target = |projectableVarSets| if projectableVars \cap postProjectionVars = \emptyset else clusterIndex < target < |projectableVarSets|
  Int getNodeRank(
    const vector<Int>& restrictedVarOrder,
    string clusteringHeuristic
  ); // rank = |restrictedVarOrder| if restrictedVarOrder \cap postProjectionVars = \emptyset else 0 \le rank < |restrictedVarOrder|
  bool isTerminal() const;
};

class JoinTerminal : public JoinNode {
public:
  Int getWidth(const Assignment& assignment = Assignment()) const override;

  void updateVarSizes(Map<Int, size_t>& varSizes, const Cnf& cnf) const override;

  JoinTerminal(const Cnf& cnf);
};

class JoinNonterminal : public JoinNode {
public:
  Int verboseSolving = 0;
  void printNode(const string& startWord) const; // 1-indexing
  void printSubtree(const string& startWord = "") const; // post-order traversal

  Int getWidth(const Assignment& assignment = Assignment()) const override;

  void updateVarSizes(Map<Int, size_t>& varSizes, const Cnf& cnf) const override;
  vector<Int> getBiggestNodeVarOrder(const Cnf& cnf) const;
  vector<Int> getHighestNodeVarOrder() const;
  vector<Int> getLexPVarOrder(const Cnf& cnf) const;
  vector<Int> getLexPVarRanking(Graph fullPrimalGraph,Set<Int>& processedVars, Map<Int, Int> tiebreaker) const;
  vector<Int> getVarOrder(Int varOrderHeuristic, const Cnf& cnf) const;

  vector<Assignment> getOuterAssignments(Int varOrderHeuristic, Int sliceVarCount, const Cnf& cnf) const;

  JoinNonterminal(
    const vector<JoinNode*>& children,
    const Set<Int>& projectionVars = Set<Int>(),
    Int requestedNodeIndex = MIN_INT
  );
};

class JoinTree { // for JoinTreeProcessor
public:
  Int declaredVarCount = MIN_INT;
  Int declaredClauseCount = MIN_INT;
  Int declaredNodeCount = MIN_INT;

  Map<Int, JoinTerminal*> joinTerminals; // 0-indexing
  Map<Int, JoinNonterminal*> joinNonterminals; // 0-indexing

  Int width = MIN_INT; // width of latest join tree
  Float plannerDuration = 0; // cumulative time for all join trees, in seconds

  JoinNode* getJoinNode(Int nodeIndex) const; // 0-indexing
  JoinNonterminal* getJoinRoot() const;
  void printTree() const;

  JoinTree(Int declaredVarCount, Int declaredClauseCount, Int declaredNodeCount);
};

class JoinTreeProcessor {
public:
  static Int plannerPid;
  static JoinTree* joinTree;
  static JoinTree* backupJoinTree;
  static TimePoint toolStartPoint;
  static Int verboseJoinTree;
  
  const Cnf& cnf;

  Int lineIndex = 0;
  Int problemLineIndex = MIN_INT;
  Int joinTreeEndLineIndex = MIN_INT;

  static void killPlanner(); // sends SIGKILL

  /* timer: */
  static void handleSigAlrm(int signal); // kills planner after receiving SIGALRM
  static bool hasDisarmedTimer();
  static void setTimer(Float seconds); // arms or disarms timer
  static void armTimer(Float seconds); // schedules SIGALRM
  static void disarmTimer(); // in case stdin ends before timer expires

  const JoinNonterminal* getJoinTreeRoot() const;

  void processCommentLine(const vector<string>& words);
  void processProblemLine(const vector<string>& words);
  void processNonterminalLine(const vector<string>& words);

  void finishReadingJoinTree();
  void readInputStream();

  JoinTreeProcessor(Float plannerWaitDuration, const Cnf& cnf);
};
} //end namespace dpve
