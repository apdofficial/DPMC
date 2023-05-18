#include "jointrees.hpp"
#include "common.hpp"
#include "formula.hpp"
#include "io.hpp"
#include "util.hpp"
#include <iomanip>
#include <queue>
#include <signal.h>
#include <sys/time.h>

/* class JoinTree =========================================================== */

const string JOIN_TREE_WORD = "jt";
const string ELIM_VARS_WORD = "e";

using dpve::io::WARNING;
using dpve::io::DASH_LINE;
using dpve::TimePoint;
using dpve::util::MyError;
using dpve::JoinNode;
using dpve::JoinTree;
using dpve::JoinNonterminal;
using dpve::JoinTreeProcessor;
using dpve::JoinTerminal;
using dpve::Assignment;
using dpve::Label;
using dpve::Cnf;
using dpve::Int;
using dpve::Set;
using std::cout;
using std::multimap;

using std::max;
using std::min;
using std::setw;
using std::to_string;

JoinNode* JoinTree::getJoinNode(Int nodeIndex) const {
  if (joinTerminals.contains(nodeIndex)) {
    return joinTerminals.at(nodeIndex);
  }
  return joinNonterminals.at(nodeIndex);
}

JoinNonterminal* JoinTree::getJoinRoot() const {
  return joinNonterminals.at(declaredNodeCount - 1);
}

void JoinTree::printTree() const {
  cout << "c p " << JOIN_TREE_WORD << " " << declaredVarCount << " " << declaredClauseCount << " " << declaredNodeCount << "\n";
  getJoinRoot()->printSubtree("c ");
}

JoinTree::JoinTree(Int declaredVarCount, Int declaredClauseCount, Int declaredNodeCount) {
  this->declaredVarCount = declaredVarCount;
  this->declaredClauseCount = declaredClauseCount;
  this->declaredNodeCount = declaredNodeCount;
}

/* class JoinTreeProcessor ================================================== */

Int JoinTreeProcessor::plannerPid = MIN_INT;
JoinTree* JoinTreeProcessor::joinTree = nullptr;
JoinTree* JoinTreeProcessor::backupJoinTree = nullptr;
TimePoint JoinTreeProcessor::toolStartPoint;
Int JoinTreeProcessor::verboseJoinTree;

void JoinTreeProcessor::killPlanner() {
  if (plannerPid == MIN_INT) {
    cout << WARNING << "found no pid for planner process\n";
  }
  else if (kill(plannerPid, SIGKILL) == 0) {
    cout << "c killed planner process with pid " << plannerPid << "\n";
  }
  else {
    // cout << WARNING << "failed to kill planner process with pid " << plannerPid << "\n";
  }
}

void JoinTreeProcessor::handleSigAlrm(int signal) {
  assert(signal == SIGALRM);
  cout << "c received SIGALRM after " << dpve::util::getDuration(toolStartPoint) << "s\n";

  if (joinTree == nullptr && backupJoinTree == nullptr) {
    cout << "c found no join tree yet; will wait for first join tree then kill planner\n";
  }
  else {
    cout << "c found join tree; killing planner\n";
    killPlanner();
  }
}

bool JoinTreeProcessor::hasDisarmedTimer() {
  struct itimerval curr_value;
  getitimer(ITIMER_REAL, &curr_value);

  return curr_value.it_value.tv_sec == 0 && curr_value.it_value.tv_usec == 0 && curr_value.it_interval.tv_sec == 0 && curr_value.it_interval.tv_usec == 0;
}

void JoinTreeProcessor::setTimer(Float seconds) {
  assert(seconds >= 0);

  Int secs = seconds;
  Int usecs = (seconds - secs) * 1000000;

  struct itimerval new_value;
  new_value.it_value.tv_sec = secs;
  new_value.it_value.tv_usec = usecs;
  new_value.it_interval.tv_sec = 0;
  new_value.it_interval.tv_usec = 0;

  setitimer(ITIMER_REAL, &new_value, nullptr);
}

void JoinTreeProcessor::armTimer(Float seconds) {
  assert(seconds >= 0);
  signal(SIGALRM, handleSigAlrm);
  setTimer(seconds);
}

void JoinTreeProcessor::disarmTimer() {
  setTimer(0);
  cout << "c disarmed timer\n";
}

const JoinNonterminal* JoinTreeProcessor::getJoinTreeRoot() const {
  return joinTree->getJoinRoot();
}

void JoinTreeProcessor::processCommentLine(const vector<string>& words) {
  if (words.size() == 3) {
    string key = words.at(1);
    string val = words.at(2);
    if (key == "pid") {
      plannerPid = stoll(val);
    }
    else if (key == "joinTreeWidth") {
      if (joinTree != nullptr) {
        joinTree->width = stoll(val);
      }
    }
    else if (key == "seconds") {
      if (joinTree != nullptr) {
        joinTree->plannerDuration = stold(val);
      }
    }
  }
}

void JoinTreeProcessor::processProblemLine(const vector<string>& words) {
  if (problemLineIndex != MIN_INT) {
    throw MyError("multiple problem lines: ", problemLineIndex, " and ", lineIndex);
  }
  problemLineIndex = lineIndex;

  if (words.size() != 5) {
    throw MyError("problem line ", lineIndex, " has ", words.size(), " words (should be 5)");
  }

  string jtWord = words.at(1);
  if (jtWord != JOIN_TREE_WORD) {
    throw MyError("expected '", JOIN_TREE_WORD, "'; found '", jtWord, "' | line ", lineIndex);
  }

  Int declaredVarCount = stoll(words.at(2));
  Int declaredClauseCount = stoll(words.at(3));
  Int declaredNodeCount = stoll(words.at(4));

  joinTree = new JoinTree(declaredVarCount, declaredClauseCount, declaredNodeCount);

  for (Int terminalIndex = 0; terminalIndex < declaredClauseCount; terminalIndex++) {
    joinTree->joinTerminals[terminalIndex] = new JoinTerminal(cnf);
  }
}

void JoinTreeProcessor::processNonterminalLine(const vector<string>& words) {
  if (problemLineIndex == MIN_INT) {
    string message = "no problem line before internal node | line " + to_string(lineIndex);
    if (joinTreeEndLineIndex != MIN_INT) {
      message += " (previous join tree ends on line " + to_string(joinTreeEndLineIndex) + ")";
    }
    throw MyError(message);
  }

  Int parentIndex = stoll(words.front()) - 1; // 0-indexing
  if (parentIndex < joinTree->declaredClauseCount || parentIndex >= joinTree->declaredNodeCount) {
    throw MyError("wrong internal-node index | line ", lineIndex);
  }

  vector<JoinNode*> children;
  Set<Int> projectionVars;
  bool parsingElimVars = false;
  for (Int i = 1; i < words.size(); i++) {
    string word = words.at(i);
    if (word == ELIM_VARS_WORD) {
      parsingElimVars = true;
    }
    else {
      Int num = stoll(word);
      if (parsingElimVars) {
        Int declaredVarCount = joinTree->declaredVarCount;
        if (num <= 0 || num > declaredVarCount) {
          throw MyError("var '", num, "' inconsistent with declared var count '", declaredVarCount, "' | line ", lineIndex);
        }
        projectionVars.insert(num);
      }
      else {
        Int childIndex = num - 1; // 0-indexing
        if (childIndex < 0 || childIndex >= parentIndex) {
          throw MyError("child '", word, "' wrong | line ", lineIndex);
        }
        children.push_back(joinTree->getJoinNode(childIndex));
      }
    }
  }
  joinTree->joinNonterminals[parentIndex] = new JoinNonterminal(children, projectionVars, parentIndex);
}

void JoinTreeProcessor::finishReadingJoinTree() {
  Int nonterminalCount = joinTree->joinNonterminals.size();
  Int expectedNonterminalCount = joinTree->declaredNodeCount - joinTree->declaredClauseCount;

  if (nonterminalCount < expectedNonterminalCount) {
    cout << WARNING << "missing internal nodes (" << expectedNonterminalCount << " expected, " << nonterminalCount << " found) before current join tree ends on line " << lineIndex << "\n";
  }
  else {
    if (joinTree->width == MIN_INT) {
      joinTree->width = joinTree->getJoinRoot()->getWidth();
    }

    cout << "c processed join tree ending on line " << lineIndex << "\n";
    dpve::io::printRow("joinTreeWidth", joinTree->width);
    dpve::io::printRow("plannerSeconds", joinTree->plannerDuration);

    if (verboseJoinTree >= 1) {
      cout << DASH_LINE;
      joinTree->printTree();
      cout << DASH_LINE;
    }

    joinTreeEndLineIndex = lineIndex;
    if (backupJoinTree==nullptr || joinTree->width < backupJoinTree->width){
      backupJoinTree = joinTree;
      JoinNode::resetStaticFields(); 
    }
  }

  problemLineIndex = MIN_INT;
  joinTree = nullptr;
}

void JoinTreeProcessor::readInputStream() {
  string line;
  while (getline(std::cin, line)) {
    lineIndex++;

    if (verboseJoinTree >= 2) {
      dpve::io::printInputLine(line, lineIndex);
    }

    vector<string> words = dpve::util::splitInputLine(line);
    if (words.empty()) {}
    else if (words.front() == "=") { // LG's tree separator "="
      if (joinTree != nullptr) {
        finishReadingJoinTree();
      }
      if (hasDisarmedTimer()) { // timer expires before first join tree ends
        break;
      }
    }
    else if (words.front() == "c") { // possibly special comment line
      processCommentLine(words);
    }
    else if (words.front() == "p") { // problem line
      processProblemLine(words);
    }
    else { // nonterminal-node line
      processNonterminalLine(words);
    }
  }

  if (joinTree != nullptr) {
    finishReadingJoinTree();
  }

  if (!hasDisarmedTimer()) { // stdin ends before timer expires
    cout << "c stdin ends before timer expires; disarming timer\n";
    disarmTimer();
  }
}

JoinTreeProcessor::JoinTreeProcessor(Float plannerWaitDuration, const Cnf& cnf): cnf(cnf) {
  cout << "c processing join tree...\n";

  armTimer(plannerWaitDuration);
  cout << "c getting join tree from stdin with " << plannerWaitDuration << "s timer (end input with 'enter' then 'ctrl d')\n";

  readInputStream();

  if (joinTree == nullptr) {
    if (backupJoinTree == nullptr) {
      throw MyError("no join tree before line ", lineIndex);
    }
    joinTree = backupJoinTree;
    JoinNode::restoreStaticFields();
    // joinTree->printTree();
  } else{
    if (backupJoinTree != nullptr && backupJoinTree->width < joinTree->width){
      joinTree = backupJoinTree;
      JoinNode::restoreStaticFields();
    }
  }

  cout << "c getting join tree from stdin: done\n";

  if (plannerPid != MIN_INT) { // timer expires before first join tree ends
    killPlanner();
  }
}


/* class JoinNode =========================================================== */

Int JoinNode::nodeCount;
Int JoinNode::terminalCount;
Set<Int> JoinNode::nonterminalIndices;

Int JoinNode::backupNodeCount;
Int JoinNode::backupTerminalCount;
Set<Int> JoinNode::backupNonterminalIndices;

// Cnf JoinNode::cnf;

void JoinNode::resetStaticFields() {
  backupNodeCount = nodeCount;
  backupTerminalCount = terminalCount;
  backupNonterminalIndices = nonterminalIndices;

  nodeCount = 0;
  terminalCount = 0;
  nonterminalIndices.clear();
}

void JoinNode::restoreStaticFields() {
  nodeCount = backupNodeCount;
  terminalCount = backupTerminalCount;
  nonterminalIndices = backupNonterminalIndices;
}

Set<Int> JoinNode::getPostProjectionVars() const {
  return util::getDiff(preProjectionVars, projectionVars);
}

Int JoinNode::chooseClusterIndex(Int clusterIndex, const vector<Set<Int>>& projectableVarSets, string clusteringHeuristic) {
  if (clusterIndex < 0 || clusterIndex >= projectableVarSets.size()) {
    throw MyError("clusterIndex == ", clusterIndex, " whereas projectableVarSets.size() == ", projectableVarSets.size());
  }

  Set<Int> projectableVars = util::getUnion(projectableVarSets); // Z = Z_1 \cup .. \cup Z_m
  Set<Int> postProjectionVars = getPostProjectionVars(); // of this node
  if (util::isDisjoint(projectableVars, postProjectionVars)) {
    return projectableVarSets.size(); // special cluster
  }

  if (clusteringHeuristic == BUCKET_ELIM_LIST || clusteringHeuristic == BOUQUET_METHOD_LIST) {
    return clusterIndex + 1;
  }
  for (Int target = clusterIndex + 1; target < projectableVarSets.size(); target++) {
    if (!util::isDisjoint(postProjectionVars, projectableVarSets.at(target))) {
      return target;
    }
  }
  return projectableVarSets.size();
}

Int JoinNode::getNodeRank(const vector<Int>& restrictedVarOrder, string clusteringHeuristic) {
  Set<Int> postProjectionVars = getPostProjectionVars();

  if (clusteringHeuristic == BUCKET_ELIM_LIST || clusteringHeuristic == BUCKET_ELIM_TREE) { // min var rank
    Int rank = MAX_INT;
    for (Int varRank = 0; varRank < restrictedVarOrder.size(); varRank++) {
      if (postProjectionVars.contains(restrictedVarOrder.at(varRank))) {
        rank = min(rank, varRank);
      }
    }
    return (rank == MAX_INT) ? restrictedVarOrder.size() : rank;
  }

  Int rank = MIN_INT;
  for (Int varRank = 0; varRank < restrictedVarOrder.size(); varRank++) {
    if (postProjectionVars.contains(restrictedVarOrder.at(varRank))) {
      rank = max(rank, varRank);
    }
  }
  return (rank == MIN_INT) ? restrictedVarOrder.size() : rank;
}

bool JoinNode::isTerminal() const {
  return nodeIndex < terminalCount;
}

/* class JoinTerminal ======================================================= */

Int JoinTerminal::getWidth(const Assignment& assignment) const {
  return util::getDiff(preProjectionVars, assignment).size();
}

void JoinTerminal::updateVarSizes(Map<Int, size_t>& varSizes, const Cnf& cnf) const {
  Set<Int> vars = cnf.clauses.at(nodeIndex).getClauseVars();
  for (Int var : vars) {
    varSizes[var] = max(varSizes[var], vars.size());
  }
}

JoinTerminal::JoinTerminal(const Cnf& cnf) {
  nodeIndex = terminalCount;
  terminalCount++;
  nodeCount++;

  preProjectionVars = cnf.clauses.at(nodeIndex).getClauseVars();
}

/* class JoinNonterminal ===================================================== */

void JoinNonterminal::printNode(const string& startWord) const {
  cout << startWord << nodeIndex + 1 << " ";

  for (const JoinNode* child : children) {
    cout << child->nodeIndex + 1 << " ";
  }

  cout << ELIM_VARS_WORD;
  for (Int var : projectionVars) {
    cout << " " << var;
  }

  cout << "\n";
}

void JoinNonterminal::printSubtree(const string& startWord) const {
  for (const JoinNode* child : children) {
    if (!child->isTerminal()) {
      static_cast<const JoinNonterminal*>(child)->printSubtree(startWord);
    }
  }
  printNode(startWord);
}

Int JoinNonterminal::getWidth(const Assignment& assignment) const {
  Int width = util::getDiff(preProjectionVars, assignment).size();
  for (JoinNode* child : children) {
    width = max(width, child->getWidth(assignment));
  }
  return width;
}

void JoinNonterminal::updateVarSizes(Map<Int, size_t>& varSizes, const Cnf& cnf) const {
  for (Int var : preProjectionVars) {
    varSizes[var] = max(varSizes[var], preProjectionVars.size());
  }
  for (JoinNode* child : children) {
    child->updateVarSizes(varSizes, cnf);
  }
}

vector<Int> JoinNonterminal::getBiggestNodeVarOrder(const Cnf& cnf) const {
  Map<Int, size_t> varSizes; // var x |-> size of biggest node containing x
  for (Int var : cnf.apparentVars) {
    varSizes[var] = 0;
  }

  updateVarSizes(varSizes, cnf);

  multimap<size_t, Int, greater<size_t>> sizedVars = dpve::util::flipMap(varSizes); // size |-> var

  size_t prevSize = 0;

  if (verboseSolving >= 2) {
    cout << DASH_LINE;
  }

  vector<Int> varOrder;
  for (const auto& [varSize, var] : sizedVars) {
    varOrder.push_back(var);

    if (verboseSolving >= 2) {
      if (prevSize == varSize) {
        cout << " " << var;
      }
      else {
        if (prevSize > 0) {
          cout << "\n";
        }
        prevSize = varSize;
        cout << "c vars in nodes of size " << std::right << std::setw(5) << varSize << ": " << var;
      }
    }
  }

  if (verboseSolving >= 2) {
    cout << "\n" << DASH_LINE;
  }

  return varOrder;
}

vector<Int> JoinNonterminal::getHighestNodeVarOrder() const {
  vector<Int> varOrder;
  std::queue<const JoinNonterminal*> q;
  q.push(this);
  while (!q.empty()) {
    const JoinNonterminal* n = q.front();
    q.pop();
    for (Int var : n->projectionVars) {
      varOrder.push_back(var);
    }
    for (const JoinNode* child : n->children) {
      if (!child->isTerminal()) {
        q.push(static_cast<const JoinNonterminal*>(child));
      }
    }
  }
  return varOrder;
}

vector<Int> JoinNonterminal::getLexPVarOrder(const Cnf& cnf) const {
  Set<Int> processedVars;
  vector<Int> varOrder;
  Graph fullPrimalGraph = cnf.getPrimalGraph();
  vector<Int> tiebreakerinv = cnf.getMostClausesVarOrder();
  Map<Int, Int> tiebreaker;
  for (Int i = 0; i<tiebreakerinv.size(); i++){
    tiebreaker[tiebreakerinv[i]] = i; 
  }
  std::queue<const JoinNonterminal*> q;
  q.push(this);
  while (!q.empty()) {
    const JoinNonterminal* n = q.front();
    q.pop();
    auto vo = n->getLexPVarRanking(fullPrimalGraph, processedVars, tiebreaker);
    varOrder.insert(varOrder.end(),vo.begin(),vo.end());
    for (const JoinNode* child : n->children) {
      if (!child->isTerminal()) {
        q.push(static_cast<const JoinNonterminal*>(child));
      }
    }
  }
  return varOrder;
}

bool commp(const std::pair<Int,std::tuple<Label,Int>> a, const std::pair<Int,std::tuple<Label,Int>> b){
  return a.second < b.second;
}

vector<Int> JoinNonterminal::getLexPVarRanking(Graph fullPrimalGraph, Set<Int>& processedVars, Map<Int, Int> tiebreaker) const {
  Set<Int> varSet = util::getDiff(preProjectionVars, processedVars);
  Map<Int, std::tuple<Label,Int>> unnumberedVertices;
  for (Int vertex : varSet) {
    unnumberedVertices[vertex] = {Label(),tiebreaker.at(vertex)};
  }
  vector<Int> numberedVertices; // whose alpha numbers are decreasing
  Graph curGraph = fullPrimalGraph.projectOnto(varSet);
  for (Int number = varSet.size(); number > 0; number--) {
    Int vertex = max_element(unnumberedVertices.begin(), unnumberedVertices.end(), commp)->first; // ignores label
    numberedVertices.push_back(vertex);
    unnumberedVertices.erase(vertex);
    for (Int neighbor : curGraph.adjacencyMap.at(vertex)) {
      auto unnumberedNeighborIt = unnumberedVertices.find(neighbor);
      if (unnumberedNeighborIt != unnumberedVertices.end()) {
        Int unnumberedNeighbor = unnumberedNeighborIt->first;
        std::get<0>(unnumberedVertices.at(unnumberedNeighbor)).addNumber(number);
      }
    }
  }
  processedVars.insert(varSet.begin(),varSet.end());
  return numberedVertices;
}

vector<Int> JoinNonterminal::getVarOrder(Int varOrderHeuristic, const Cnf& cnf) const {
  if (CNF_VAR_ORDER_HEURISTICS.contains(abs(varOrderHeuristic))) {
    return cnf.getCnfVarOrder(varOrderHeuristic);
  }

  vector<Int> varOrder;
  if (abs(varOrderHeuristic) == dpve::BIGGEST_NODE_HEURISTIC) { // 8
    varOrder = getBiggestNodeVarOrder(cnf);
  } else if (abs(varOrderHeuristic)== 10){
    varOrder = getLexPVarOrder(cnf);
  } else {
    assert(abs(varOrderHeuristic) == HIGHEST_NODE_HEURISTIC); //9
    varOrder = getHighestNodeVarOrder();
  }

  if (varOrderHeuristic < 0) {
    reverse(varOrder.begin(), varOrder.end());
  }

  return varOrder;
}

vector<Assignment> JoinNonterminal::getOuterAssignments(Int varOrderHeuristic, Int sliceVarCount, const Cnf& cnf) const {
  if (sliceVarCount <= 0) {
    return {Assignment()};
  }

  TimePoint sliceVarOrderStartPoint = util::getTimePoint();
  vector<Int> varOrder = getVarOrder(varOrderHeuristic, cnf);
  if (verboseSolving >= 1) {
    io::printRow("sliceVarSeconds", util::getDuration(sliceVarOrderStartPoint));
  }

  TimePoint assignmentsStartPoint = util::getTimePoint();
  vector<Assignment> assignments;

  if (verboseSolving >= 2) {
    cout << "c slice var order heuristic: {";
  }

  for (Int i = 0, assignedVars = 0; i < varOrder.size() && assignedVars < sliceVarCount; i++) {
    Int var = varOrder.at(i);
    if (cnf.outerVars.contains(var)) {
      assignments = Assignment::getExtendedAssignments(assignments, var);
      assignedVars++;
      if (verboseSolving >= 2) {
        cout << " " << var;
      }
    }
  }

  if (verboseSolving >= 2) {
    cout << " }\n";
  }

  if (verboseSolving >= 1) {
    dpve::io::printRow("sliceAssignmentsSeconds", dpve::util::getDuration(assignmentsStartPoint));
  }

  return assignments;
}

JoinNonterminal::JoinNonterminal(const vector<JoinNode*>& children, const Set<Int>& projectionVars, Int requestedNodeIndex) {
  this->children = children;
  this->projectionVars = projectionVars;

  if (requestedNodeIndex == MIN_INT) {
    requestedNodeIndex = nodeCount;
  }
  else if (requestedNodeIndex < terminalCount) {
    throw MyError("requestedNodeIndex == ", requestedNodeIndex, " < ", terminalCount, " == terminalCount");
  }
  else if (nonterminalIndices.contains(requestedNodeIndex)) {
    throw MyError("requestedNodeIndex ", requestedNodeIndex, " already taken");
  }

  nodeIndex = requestedNodeIndex;
  nonterminalIndices.insert(nodeIndex);
  nodeCount++;

  for (JoinNode* child : children) {
    util::unionize(preProjectionVars, child->getPostProjectionVars());
  }
}


