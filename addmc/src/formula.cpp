#include "formula.hpp"
#include "common.hpp"
#include "../libraries/colamd/colamd.h"
#include "io.hpp"
#include "util.hpp"
#include <fstream>
#include <iomanip>
#include <random>
/* class Clause ============================================================= */

// using dpve::util::MyError;
// using dpve::Clause;
// using dpve::Assignment;
// using dpve::Cnf;

using dpve::Number;
using std::cout;
using std::min;
using std::max;
using std::next;
using std::right;
using std::setw;
using std::to_string;

extern ostream& dpve::operator<<(ostream& stream, const Number& n);

namespace dpve{
using util::MyError;
using util::EmptyClauseException;
using io::printRow;

Clause::Clause(bool xorFlag) {
  this->xorFlag = xorFlag;
}

void Clause::insertLiteral(Int literal) {
  if (xorFlag && contains(literal)) {
    erase(literal);
  }
  else {
    insert(literal);
  }
}

void Clause::printClause() const {
  cout << (xorFlag ? " x" : "  ");
  for (auto it = begin(); it != end(); it++) {
    cout << " " << right << setw(5) << *it;
  }
  cout << "\n";
}

Set<Int> Clause::getClauseVars() const {
  Set<Int> vars;
  for (Int literal : *this) {
    vars.insert(abs(literal));
  }
  return vars;
}

/* class Cnf ================================================================ */

Set<Int> Cnf::getInnerVars() const {
  Set<Int> innerVars;
  for (Int var = 1; var <= declaredVarCount; var++) {
    if (!outerVars.contains(var)) {
      innerVars.insert(var);
    }
  }
  return innerVars;
}

Map<Int, Number> Cnf::getUnprunableWeights() const {
  Map<Int, Number> unprunableWeights;
  for (const auto& [literal, weight] : literalWeights) {
    if (weight > Number("1")) {
      unprunableWeights[literal] = weight;
    }
  }
  return unprunableWeights;
}

void Cnf::printLiteralWeight(Int literal, const Number& weight) {
  cout << "c  weight " << right << setw(5) << literal << ": " << weight << "\n";
}

void Cnf::printLiteralWeights() const {
  cout << "c literal weights:\n";
  for (Int var = 1; var <= declaredVarCount; var++) {
    printLiteralWeight(var, literalWeights.at(var));
    printLiteralWeight(-var, literalWeights.at(-var));
  }
}

void Cnf::printClauses() const {
  cout << "c CNF formula:\n";
  for (Int i = 0; i < clauses.size(); i++) {
    cout << "c  clause " << right << setw(5) << i + 1 << ":";
    clauses.at(i).printClause();
  }
}

void Cnf::addClause(const Clause& clause) {
  Int clauseIndex = clauses.size();
  clauses.push_back(clause);
  for (Int literal : clause) {
    Int var = abs(literal);
    auto it = varToClauses.find(var);
    if (it != varToClauses.end()) {
      it->second.insert(clauseIndex);
    }
    else {
      varToClauses[var] = {clauseIndex};
    }
  }
}

void Cnf::setApparentVars() {
  for (const pair<Int, Set<Int>>& kv : varToClauses) {
    apparentVars.insert(kv.first);
  }
}

Graph Cnf::getPrimalGraph() const {
  Graph graph(apparentVars);
  for (const Clause& clause : clauses) {
    for (auto literal1 = clause.begin(); literal1 != clause.end(); literal1++) {
      for (auto literal2 = next(literal1); literal2 != clause.end(); literal2++) {
        Int var1 = abs(*literal1);
        Int var2 = abs(*literal2);
        graph.addEdge(var1, var2);
      }
    }
  }
  return graph;
}

vector<Int> Cnf::getRandomVarOrder() const {
  vector<Int> varOrder(apparentVars.begin(), apparentVars.end());
  std::mt19937 generator;
  generator.seed(randomSeed);
  shuffle(varOrder.begin(), varOrder.end(), generator);
  return varOrder;
}

vector<Int> Cnf::getDeclarationVarOrder() const {
  vector<Int> varOrder;
  for (Int var = 1; var <= declaredVarCount; var++) {
    if (apparentVars.contains(var)) {
      varOrder.push_back(var);
    }
  }
  return varOrder;
}

vector<Int> Cnf::getMostClausesVarOrder() const {
  multimap<Int, Int, greater<Int>> m; // clause count |-> var
  for (const auto& [var, clauseIndices] : varToClauses) {
    m.insert({clauseIndices.size(), var});
  }

  vector<Int> varOrder;
  for (const auto& [clauseCount, var] : m) {
    varOrder.push_back(var);
  }
  return varOrder;
}

vector<Int> Cnf::getMinFillVarOrder() const {
  vector<Int> varOrder;

  Graph graph = getPrimalGraph();
  while (!graph.vertices.empty()) {
    Int vertex = graph.getMinFillVertex();
    graph.fillInEdges(vertex);
    graph.removeVertex(vertex);
    varOrder.push_back(vertex);
  }

  return varOrder;
}

vector<Int> Cnf::getMcsVarOrder() const {
  Graph graph = getPrimalGraph();

  auto startVertex = graph.vertices.begin();
  if (startVertex == graph.vertices.end()) {
    return vector<Int>();
  }

  Map<Int, Int> rankedNeighborCounts; // unranked vertex |-> number of ranked neighbors
  for (auto it = next(startVertex); it != graph.vertices.end(); it++) {
    rankedNeighborCounts[*it] = 0;
  }

  Int bestVertex = *startVertex;
  Int bestRankedNeighborCount = MIN_INT;

  vector<Int> varOrder;
  do {
    varOrder.push_back(bestVertex);

    rankedNeighborCounts.erase(bestVertex);

    for (Int n : graph.adjacencyMap.at(bestVertex)) {
      auto entry = rankedNeighborCounts.find(n);
      if (entry != rankedNeighborCounts.end()) {
        entry->second++;
      }
    }

    bestRankedNeighborCount = MIN_INT;
    for (pair<Int, Int> entry : rankedNeighborCounts) {
      if (entry.second > bestRankedNeighborCount) {
        bestRankedNeighborCount = entry.second;
        bestVertex = entry.first;
      }
    }
  }
  while (bestRankedNeighborCount != MIN_INT);

  return varOrder;
}

vector<Int> Cnf::getLexPVarOrder() const {
  Map<Int, Label> unnumberedVertices;
  for (Int vertex : apparentVars) {
    unnumberedVertices[vertex] = Label();
  }
  vector<Int> numberedVertices; // whose alpha numbers are decreasing
  Graph graph = getPrimalGraph();
  for (Int number = apparentVars.size(); number > 0; number--) {
    Int vertex = max_element(unnumberedVertices.begin(), unnumberedVertices.end(), Label::hasSmallerLabel)->first; // ignores label
    numberedVertices.push_back(vertex);
    unnumberedVertices.erase(vertex);
    for (Int neighbor : graph.adjacencyMap.at(vertex)) {
      auto unnumberedNeighborIt = unnumberedVertices.find(neighbor);
      if (unnumberedNeighborIt != unnumberedVertices.end()) {
        Int unnumberedNeighbor = unnumberedNeighborIt->first;
        unnumberedVertices.at(unnumberedNeighbor).addNumber(number);
      }
    }
  }
  return numberedVertices;
}

vector<Int> Cnf::getLexMVarOrder() const {
  Map<Int, Label> unnumberedVertices;
  for (Int vertex : apparentVars) {
    unnumberedVertices[vertex] = Label();
  }
  vector<Int> numberedVertices; // whose alpha numbers are decreasing
  Graph graph = getPrimalGraph();
  for (Int i = apparentVars.size(); i > 0; i--) {
    Int v = max_element(unnumberedVertices.begin(), unnumberedVertices.end(), Label::hasSmallerLabel)->first; // ignores label
    numberedVertices.push_back(v);
    unnumberedVertices.erase(v);

    /* updates numberedVertices: */
    Graph subgraph = getPrimalGraph(); // will only contain v, w, and unnumbered vertices whose labels are less than w's label
    for (pair<const Int, Label>& wAndLabel : unnumberedVertices) {
      Int w = wAndLabel.first;
      Label& wLabel = wAndLabel.second;

      /* removes numbered vertices except v: */
      for (Int numberedVertex : numberedVertices) {
        if (numberedVertex != v) {
          subgraph.removeVertex(numberedVertex);
        }
      }

      /* removes each non-w unnumbered vertex whose label is at least w's labels */
      for (const pair<Int, Label>& kv : unnumberedVertices) {
        Int unnumberedVertex = kv.first;
        const Label& label = kv.second;
        if (unnumberedVertex != w && label >= wLabel) {
          subgraph.removeVertex(unnumberedVertex);
        }
      }

      if (subgraph.hasPath(v, w)) {
        wLabel.addNumber(i);
      }
    }
  }
  return numberedVertices;
}

vector<Int> Cnf::getColAMDVarOrder() const {
  /* create matrix A, p of entries of constraints 
      
  */
  // int A[10], p[10];
  vector<uint64_t> p, A;
  uint64_t nNZ = 0;
  uint64_t nRows = clauses.size();
  uint64_t nCols = apparentVars.size();
  p.push_back(0);
  vector<Int> colToVar;
  for (const auto& [var, clauseIndices] : varToClauses) {
    uint32_t nCls = clauseIndices.size();
    p.push_back(p.back()+nCls);
    for (auto ind: clauseIndices){
      A.push_back(ind);
    }
    nNZ += nCls;
    colToVar.push_back(var);
  }
  size_t ALEN = colamd_recommended(nRows,nCols,nNZ) ;
  int* Aa = new int[ALEN];
  int* pp = new int[p.size()];
  std::copy(A.begin(),A.end(),Aa);
  std::copy(p.begin(),p.end(),pp);
  int stats [COLAMD_STATS] ;	/* for colamd and symamd output statistics */
  int ok;
  ok = colamd (nRows, nCols, ALEN, Aa, pp, (double *) NULL, stats) ;
  // colamd_report (stats) ;
  if (!ok){
	  cout<<"c colamd Error! Exiting..\n";
	  exit (1) ;
  }
  // pp[0] = j means that column j of original matrix A should be in 0th position after permutation
  // colToVar[j] = k means that jth column of original matrix A belonged to the apparentVar k
  // so colToVar[pp[0]] will give the apparentVar to place in 0th positon in varOrder
  vector<Int> varOrder;
  for(uint64_t i=0; i<apparentVars.size();i++) {
    varOrder.push_back(colToVar[pp[i]]);
  }
  delete Aa;
  delete pp;
  return(varOrder);
}

vector<Int> Cnf::getCnfVarOrder(Int cnfVarOrderHeuristic) const {
  vector<Int> varOrder;
  switch (abs(cnfVarOrderHeuristic)) {
    case RANDOM_HEURISTIC:
      varOrder = getRandomVarOrder();
      break;
    case DECLARATION_HEURISTIC:
      varOrder = getDeclarationVarOrder();
      break;
    case MOST_CLAUSES_HEURISTIC:
      varOrder = getMostClausesVarOrder();
      break;
    case MIN_FILL_HEURISTIC:
      varOrder = getMinFillVarOrder();
      break;
    case MCS_HEURISTIC:
      varOrder = getMcsVarOrder();
      break;
    case LEX_P_HEURISTIC:
      varOrder = getLexPVarOrder();
      break;
    case COLAMD_HEURISTIC:
      varOrder = getColAMDVarOrder();
      break;
    default:
      assert(abs(cnfVarOrderHeuristic) == LEX_M_HEURISTIC);
      varOrder = getLexMVarOrder();
  }

  if (cnfVarOrderHeuristic < 0) {
    reverse(varOrder.begin(), varOrder.end());
  }

  return varOrder;
}

bool Cnf::isMc21ShowLine(const vector<string> &words) const {
  return words.size() >= 4 && words.front() == "c" && words.at(1) == "p" && words.at(2) == "show";
}

bool Cnf::isMc21WeightLine(const vector<string> &words) const {
  bool b = words.size() >= 3 && words.front() == "c" && words.at(1) == "p" && words.at(2) == "weight";
  switch (words.size()) {
    case 5:
      return b;
    case 6:
      return b && words.back() == "0";
    default:
      return false;
  }
}

void Cnf::completeLiteralWeights() {
  if (weightedCounting) {
    for (Int var = 1; var <= declaredVarCount; var++) {
      if (!literalWeights.contains(var) && !literalWeights.contains(-var)) {
        literalWeights[var] = Number("1");
        literalWeights[-var] = Number("1");
      }
      else if (!literalWeights.contains(var)) {
        assert(literalWeights.at(-var) < Number("1"));
        literalWeights[var] = Number("1") - literalWeights.at(-var);
      }
      else if (!literalWeights.contains(-var)) {
        assert(literalWeights.at(var) < Number("1"));
        literalWeights[-var] = Number("1") - literalWeights.at(var);
      }
    }
  }
  else {
    for (Int var = 1; var <= declaredVarCount; var++) {
      literalWeights[var] = Number("1");
      literalWeights[-var] = Number("1");
    }
  }
}

void Cnf::printStats() const {
  Float clauseSizeSum = 0;
  Int clauseSizeMax = MIN_INT;
  Int clauseSizeMin = MAX_INT;

  for (const Clause& clause : clauses) {
    Int clauseSize = clause.size();
    clauseSizeSum += clauseSize;
    clauseSizeMax = max(clauseSizeMax, clauseSize);
    clauseSizeMin = min(clauseSizeMin, clauseSize);
  }

  dpve::io::printRow("clauseSizeMean", clauseSizeSum / clauses.size());
  dpve::io::printRow("clauseSizeMax", clauseSizeMax);
  dpve::io::printRow("clauseSizeMin", clauseSizeMin);
}

void Cnf::readCnfFile(const string& filePath) {
  cout << "c processing CNF formula...\n";

  std::ifstream inputFileStream(filePath);
  if (!inputFileStream.is_open()) {
    throw MyError("unable to open file '", filePath, "'");
  }

  Int declaredClauseCount = MIN_INT;

  Int lineIndex = 0;
  Int problemLineIndex = MIN_INT;

  string line;
  while (getline(inputFileStream, line)) {
    lineIndex++;
    std::istringstream inStringStream(line);

    if (verboseCnf >= 3) {
      io::printInputLine(line, lineIndex);
    }

    vector<string> words = util::splitInputLine(line);
    if (words.empty()) {
      continue;
    }
    string& frontWord = words.front();
    if (Set<string>{"s", "INDETERMINATE"}.contains(frontWord)) { // preprocessor pmc
      throw MyError("unexpected output from preprocessor pmc | line ", lineIndex, ": ", line);
    }
    else if (frontWord == "p") { // problem line
      if (problemLineIndex != MIN_INT) {
        throw MyError("multiple problem lines: ", problemLineIndex, " and ", lineIndex);
      }

      problemLineIndex = lineIndex;

      if (words.size() != 4) {
        throw MyError("problem line ", lineIndex, " has ", words.size(), " words (should be 4)");
      }

      declaredVarCount = stoll(words.at(2));
      declaredClauseCount = stoll(words.at(3));
    }
    else if (frontWord == "c") { // possibly show or weight line
      if (projectedCounting && isMc21ShowLine(words)) {
        if (problemLineIndex == MIN_INT) {
          throw MyError("no problem line before outer vars | line ", lineIndex, ": ", line);
        }

        for (Int i = 3; i < words.size(); i++) {
          Int num = stoll(words.at(i));
          if (num == 0) {
            if (i != words.size() - 1) {
              throw MyError("outer vars terminated prematurely by '0' | line ", lineIndex);
            }
          }
          else if (num < 0 || num > declaredVarCount) {
            throw MyError("var '", num, "' inconsistent with declared var count '", declaredVarCount, "' | line ", lineIndex);
          }
          else {
            outerVars.insert(num);
          }
        }
      }
      else if (weightedCounting && isMc21WeightLine(words)) {
        if (problemLineIndex == MIN_INT) {
          throw MyError("no problem line before literal weight | line ", lineIndex, ": ", line);
        }

        Int literal = stoll(words.at(3));
        assert(literal != 0);

        if (abs(literal) > declaredVarCount) {
          throw MyError("literal '", literal, "' inconsistent with declared var count '", declaredVarCount, "' | line ", lineIndex);
        }

        Number weight(words.at(4));
        if (weight <= Number()) {
          throw MyError("weight must be positive | line ", lineIndex);
        }
        literalWeights[literal] = weight;
      }
    }
    else if (!frontWord.starts_with("c")) { // clause line
      if (problemLineIndex == MIN_INT) {
        throw MyError("no problem line before clause | line ", lineIndex);
      }

      bool xorFlag = false;
      if (frontWord.starts_with("x")) {
        xorFlag = true;
        xorClauseCount++;

        if (frontWord == "x") {
          words.erase(words.begin());
        }
        else {
          frontWord.erase(frontWord.begin());
        }
      }
      Clause clause(xorFlag);

      for (Int i = 0; i < words.size(); i++) {
        Int num = stoll(words.at(i));

        if (abs(num) > declaredVarCount) {
          throw MyError("literal '", num, "' inconsistent with declared var count '", declaredVarCount, "' | line ", lineIndex);
        }

        if (num == 0) {
          if (i != words.size() - 1) {
            throw MyError("clause terminated prematurely by '0' | line ", lineIndex);
          }

          if (clause.empty()) {
            throw EmptyClauseException(lineIndex, line);
          }

          addClause(clause);
        }
        else { // literal
          if (i == words.size() - 1) {
            throw MyError("missing end-of-clause indicator '0' | line ", lineIndex);
          }
          clause.insertLiteral(num);
        }
      }
    }
  }

  if (problemLineIndex == MIN_INT) {
    throw MyError("no problem line before CNF file ends on line ", lineIndex);
  }

  setApparentVars();

  if (!projectedCounting) {
    for (Int var = 1; var <= declaredVarCount; var++) {
      outerVars.insert(var);
    }
  }

  completeLiteralWeights();

  if (verboseCnf >= 1) {
    printRow("declaredVarCount", declaredVarCount);
    printRow("apparentVarCount", apparentVars.size());

    printRow("declaredClauseCount", declaredClauseCount);
    printRow("apparentClauseCount", clauses.size());
    printRow("xorClauseCount", xorClauseCount);

    printStats();

    if (verboseCnf >= 2) {
      if (projectedCounting) {
        cout << "c outer vars: { ";
        for (Int var : util::getSortedNums(outerVars)) {
          cout << var << " ";
        }
        cout << "}\n";
      }

      if (weightedCounting) {
        printLiteralWeights();
      }

      printClauses();
    }
  }

  cout << "\n";
}

Cnf::Cnf(const Int verboseCnf, const Int randomSeed, const bool weightedCounting, const bool projectedCounting):
verboseCnf(verboseCnf), randomSeed(randomSeed), weightedCounting(weightedCounting), projectedCounting(projectedCounting)
{}

/* classes for join trees =================================================== */

/* class Assignment ========================================================= */

Assignment::Assignment() {}

Assignment::Assignment(Int var, bool val) {
  insert({var, val});
}

Assignment::Assignment(const string& bitString) {
  for (Int i = 0; i < bitString.size(); i++) {
    char bit = bitString.at(i);
    assert(bit == '0' or bit == '1');
    insert({i + 1, bit == '1'});
  }
}

bool Assignment::getValue(Int var) const {
  auto it = find(var);
  return (it != end()) ? it->second : true;
}

void Assignment::printAssignment() const {
  for (auto it = begin(); it != end(); it++) {
    cout << right << setw(5) << (it->second ? it->first : -it->first);
    if (next(it) != end()) {
      cout << " ";
    }
  }
}

vector<Assignment> Assignment::getExtendedAssignments(const vector<Assignment>& assignments, Int var) {
  vector<Assignment> extendedAssignments;
  if (assignments.empty()) {
    extendedAssignments.push_back(Assignment(var, false));
    extendedAssignments.push_back(Assignment(var, true));
  }
  else {
    for (Assignment assignment : assignments) {
      assignment[var] = false;
      extendedAssignments.push_back(assignment);
      assignment[var] = true;
      extendedAssignments.push_back(assignment);
    }
  }
  return extendedAssignments;
}

string Assignment::getShortFormat(Int declaredVarCount) const{
  string s;
  for (Int cnfVar = 1; cnfVar <= declaredVarCount; cnfVar++) {
    s += to_string(getValue(cnfVar));
  }
  return s;
}

string Assignment::getLongFormat(Int declaredVarCount) const{
  string s;
  for (Int cnfVar = 1; cnfVar <= declaredVarCount; cnfVar++) {
    s += (getValue(cnfVar) ? " " : " -") + to_string(cnfVar);
  }
  return s;
}
} //end namespace dpve