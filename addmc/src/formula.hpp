#pragma once

#include "types.hpp"
#include "graph.hpp"

#include <vector>
#include <string>

namespace dpve{

class Clause : public Set<Int> {
public:
  bool xorFlag;

  Clause(bool xorFlag);

  void insertLiteral(Int literal);

  void printClause() const;
  Set<Int> getClauseVars() const;
};

class Cnf {
public:
  Int declaredVarCount = 0;
  Set<Int> outerVars;
  Map<Int, Number> literalWeights; // for outer and inner vars
  vector<Clause> clauses;
  Int xorClauseCount = 0;

  const Int verboseCnf;
  const Int randomSeed;
  const bool weightedCounting;
  const bool projectedCounting;

  Set<Int> apparentVars; // as opposed to hidden vars that are declared but appear in no clause
  Map<Int, Set<Int>> varToClauses; // apparent var |-> clause indices

  Set<Int> getInnerVars() const;
  Map<Int, Number> getUnprunableWeights() const;

  static void printLiteralWeight(Int literal, const Number& weight);
  void printLiteralWeights() const;
  void printClauses() const;

  void addClause(const Clause& clause);
  void setApparentVars();
  Graph getPrimalGraph() const;
  vector<Int> getRandomVarOrder() const;
  vector<Int> getDeclarationVarOrder() const;
  vector<Int> getMostClausesVarOrder() const;
  vector<Int> getMinFillVarOrder() const;
  vector<Int> getMcsVarOrder() const;
  vector<Int> getLexPVarOrder() const;
  vector<Int> getLexMVarOrder() const;
  vector<Int> getCnfVarOrder(Int cnfVarOrderHeuristic) const;
  vector<Int> getColAMDVarOrder() const;

  bool isMc21ShowLine(const vector<string> &words) const; // c p show <vars> [0]
  bool isMc21WeightLine(const vector<string> &words) const; // c p weight <literal> <weight> [0]

  void completeLiteralWeights();

  void printStats() const;

  void readCnfFile(const string& filePath);

  Cnf(const Int verboseCnf, const Int randomSeed, const bool weightedCounting, const bool projectedCounting); // empty conjunction
  Cnf();
  private:
  // Cnf();
};

/* classes for join trees =================================================== */

class Assignment : public Map<Int, bool> { // partial var assignment
public:
  Assignment();
  Assignment(Int var, bool val);
  Assignment(const string& bitString);

  bool getValue(Int var) const; // returns `true` if `var` is unassigned
  void printAssignment() const;
  string getLongFormat(Int declaredVarCount) const;
  string getShortFormat(Int declaredVarCount) const;
  static vector<Assignment> getExtendedAssignments(const vector<Assignment>& assignments, Int var);
};
} //end namespace dpve