#include "dmc.hpp"
#include "util.hpp"
#include "util/util.h"

#include <queue>
#include <tuple>
/* class SatFilter =========================================================== */

using std::tuple;

using dpve::io::printRow;
using dpve::io::printLine;
using dpve::Dd;
using dpve::Number;
using dpve::SatFilter;
using dpve::Executor;
using dpve::Dpve;
using dpve::Assignment;

using std::max;
using std::to_string;

Dd SatFilter::getClauseBdd(const Clause& clause){
  Dd clauseDd = Dd::getZeroBdd();
  for (Int literal : clause) {
    bool val = literal > 0;
    Int cnfVar = abs(literal);
    Int ddVar = cnfVarToDdVarMap.at(cnfVar);
    Dd literalDd = Dd::getVarBdd(ddVar, val);
    clauseDd = clauseDd.getBddOr(literalDd);
  }
  return clauseDd;
}

Dd SatFilter::solveSubtree(const JoinNode* joinNode) {
  if (joinNode->isTerminal()) {
    TimePoint terminalStartPoint = dpve::util::getTimePoint();
    // cout << "c 1\n";
    Dd d = getClauseBdd(cnf.clauses.at(joinNode->nodeIndex));
    ((JoinNode*) joinNode)->dd = (void*) new Dd(Dd::getOneBdd()); //cast away const-ness of joinNode to modify dd. Should be fine since object itself is not const
    joinNodesProcessed ++;
    if (((joinNodesProcessed-1)%(std::max((JoinNode::nodeCount/10),1LL)))==1) printLine(to_string(joinNodesProcessed)+"/"+to_string(JoinNode::nodeCount)+":"+to_string(util::getDuration(satFilterStartPoint))+" ");
    return d;
  }

  // cout << "c 3\n";
  Dd prod = Dd::getOneBdd();
  
  if (joinPriority == FCFS){
    for (JoinNode* child : joinNode->children) {
      prod = prod.getBddAnd(solveSubtree(child));
    }  
  } else {
    vector<Dd> childDdList;
    for (JoinNode* child : joinNode->children) {
      // cout << "c before recurse\n";
      childDdList.push_back(solveSubtree(child));
    }

    if (joinPriority == ARBITRARY_PAIR) { // arbitrarily multiplies child decision diagrams
      for (Dd childDd : childDdList) {
        prod = prod.getBddAnd(childDd);
      }
    }
    else { //define lambda to be used for priority queue comparison operation
      auto compare = [this](Dd leftDd, Dd rightDd) {
        if (joinPriority == SMALLEST_PAIR) { // top = rightmost = smallest
          return leftDd.getNodeCount() > rightDd.getNodeCount();
        }
        return leftDd.getNodeCount() < rightDd.getNodeCount(); 
      };
      std::priority_queue<Dd,vector<Dd>, decltype(compare)> childDdQueue(compare);
      childDdQueue.push(prod);
      for (Dd childDd : childDdList) {
        childDdQueue.push(childDd);
      }
      assert(!childDdQueue.empty());
      // Dd::manualReorder();
      while (childDdQueue.size() >= 2) {
        Dd dd1 = childDdQueue.top();
        childDdQueue.pop();
        Dd dd2 = childDdQueue.top();
        childDdQueue.pop();
        Dd dd3 = dd1.getBddAnd(dd2);
        childDdQueue.push(dd3);
      }
      prod = childDdQueue.top();
    }
  }
  Dd retDD = Dd::getOneDd();
  if (joinNode->projectionVars.size()>0){
    ((JoinNode*) joinNode)->dd = (void*) new Dd(prod); //to make sure it does not go out of scope since we are storing a pointer
    //cast away const-ness of joinNode to modify dd. Should be fine since object itself is not const
    
    vector<Int> ddVars;
    for (Int cnfVar : joinNode->projectionVars) {
      ddVars.push_back(cnfVarToDdVarMap.at(cnfVar));
    }
    // cout << "c 7\n";
    retDD = prod.getBddExists(ddVars, ddVarToCnfVarMap);
  } else{
    ((JoinNode*) joinNode)->dd = (void*) new Dd(Dd::getOneBdd()); //cast away const-ness of joinNode to modify dd. Should be fine since object itself is not const
    retDD = prod;
  }
  joinNodesProcessed ++;
  if (((joinNodesProcessed-1)%(std::max((JoinNode::nodeCount/10),1LL)))==1) printLine(to_string(joinNodesProcessed)+"/"+to_string(JoinNode::nodeCount)+":"+to_string(util::getDuration(satFilterStartPoint))+" ");
  // cout << "c 8\n";
  return retDD;
}

bool SatFilter::filterBdds(const JoinNode* joinNode, Dd parentBdd){
  Dd* d = (Dd*)joinNode->dd;
  bool hasNewClsDescendents = false;
  bool bottomMost = true; //joinnode with projvars.size>0 and at least one jointerminal descendant. 
  //joinnode->dd should not be set to one if true, as it will be used by executor.
  if (joinNode->isTerminal()){
    hasNewClsDescendents = true;
  }
  if (joinNode->projectionVars.size()==0){
    d = &parentBdd;   // okay since *d gets passed to filterBdds below
    bottomMost = false;
  } else{
    *(Dd*)(((JoinNode*) joinNode)->dd) = d->getFilteredBdd(parentBdd);
    assert(d == (Dd*)joinNode->dd); 
    bottomMost = true;
    hasNewClsDescendents = false;
  }
  // if bottommost is false it should remain false no matter what
  // if it is true, it should become false if even one other value is true
  // so take OR of return values of children, take negation of that, and take and of that
  // by demorgans law it'll be bottomost &= !filterbdds
  for (JoinNode* child : joinNode->children) {
    hasNewClsDescendents |= (filterBdds(child, *d)); // bitwise & so that no short-circuiting happens
  }
  if (hasNewClsDescendents && bottomMost){
    //do nothing
  } else{
    *(Dd*)(((JoinNode*) joinNode)->dd) = Dd::getOneBdd();
  }
  joinNodesProcessed ++;
  if(joinNodesProcessed>JoinNode::nodeCount){
    joinNodesProcessed=1;
  }
  if (((joinNodesProcessed-1)%(std::max((JoinNode::nodeCount/10),1LL)))==1) printLine(to_string(joinNodesProcessed)+"/"+to_string(JoinNode::nodeCount)+":"+to_string(util::getDuration(satFilterStartPoint))+" ");
  return hasNewClsDescendents;
}

SatFilter::SatFilter(const Cnf& cnf, const Map<Int, Int>& cnfVarToDdVarMap,const vector<Int>& ddVarToCnfVarMap):
cnf(cnf), cnfVarToDdVarMap(cnfVarToDdVarMap), ddVarToCnfVarMap(ddVarToCnfVarMap)
{
  joinNodesProcessed = 0;
  satFilterStartPoint = util::getTimePoint();
}

/* class Executor =========================================================== */

Dd Executor::solveSubtree(const JoinNode* joinNode, const PruneMaxParams& pmParams, const Assignment& assignment ) {
  // cout<<"Starting visit of joinNode number "<<joinNode->nodeIndex+1<<"\n";
  if (joinNode->isTerminal()) {
    TimePoint terminalStartPoint = util::getTimePoint();

    joinNodesProcessed ++;
    if (((joinNodesProcessed-1)%(std::max((JoinNode::nodeCount/10),1LL)))==1) printLine(to_string(joinNodesProcessed)+"/"+to_string(JoinNode::nodeCount)+":"+to_string(util::getDuration(executorStartPoint))+" ");
    if (satFilter>0) {
      return (((Dd*)joinNode->dd)->getAdd());
    } else{
      Map<Int, pair<Int,Int>> ddVarSignAndAsnmts;
      const Clause& c = cnf.clauses.at(joinNode->nodeIndex);
      for (auto cnfLit:c){
        Int ddVar = cnfVarToDdVarMap.at(abs(cnfLit));
        //NOTE: ddVars are 0-indexed so can't use sign to indicate polarity
        auto it = assignment.find(abs(cnfLit));
        if (it != assignment.end()) { // variable has assigned value in assignment
          ddVarSignAndAsnmts[ddVar] = {cnfLit>0,it->second?1:-1}; //if asnmt is true then positive else negative
        } else{
          ddVarSignAndAsnmts[ddVar] = {cnfLit>0,0}; // unassigned then zero
        }
      }
      return Dd::getClauseDd(ddVarSignAndAsnmts,c.xorFlag);
    }  
  }

  Dd dd = satFilter>0? ((Dd*)joinNode->dd)->getAdd() : Dd::getOneDd();
  if (satFilter>0){
    *(Dd*)(((JoinNode*) joinNode)->dd) = Dd::getOneBdd(); //once you get the ADD no need for the BDD
  }

  if (joinPriority == FCFS){
    for (JoinNode* child : joinNode->children) {
      dd = dd.getProduct(solveSubtree(child, pmParams, assignment));
    }
  } else{
    vector<Dd> childDdList;
    for (JoinNode* child : joinNode->children) {
      childDdList.push_back(solveSubtree(child, pmParams, assignment));
    }
    
    //Following call considers reordering if enabled.
    //Has internal checks to decide when to reorder
    // Dd::manualReorder(levelMaps);
    /*
    if (joinNodesProcessed>JoinNode::nodeCount*reOrdThresh){
      reOrdThresh += 0.05;
      
      if(voInd == 0){
        assert (ratio == 0);
        printLine("No better order found. Order is unchanged!");
      } else{
        printLine("Better order found: "+CNF_VAR_ORDER_HEURISTICS.at(abs(voInd))+(voInd<0?" reversed": " unreverseed")+
          " with size reduction ratio: "+to_string(ratio));
      }
    }
    */

    if (joinPriority == ARBITRARY_PAIR) { // arbitrarily multiplies child decision diagrams
      for (Dd childDd : childDdList) {
        dd = dd.getProduct(childDd);
      }
    }
    else { 
      //define lambda to be used for priority queue comparison operation
      auto compare = [this](Dd leftDd, Dd rightDd) {
        if (joinPriority == SMALLEST_PAIR) { // top = rightmost = smallest
          return leftDd.getNodeCount() > rightDd.getNodeCount();
        }
        return leftDd.getNodeCount() < rightDd.getNodeCount(); 
      };
      std::priority_queue<Dd,vector<Dd>, decltype(compare)> childDdQueue(compare);
      childDdQueue.push(dd);
      for (Dd childDd : childDdList) {
        childDdQueue.push(childDd);
      }
      assert(!childDdQueue.empty());
      while (childDdQueue.size() >= 2) {
        Dd dd1 = childDdQueue.top();
        childDdQueue.pop();
        Dd dd2 = childDdQueue.top();
        childDdQueue.pop();
        Dd dd3 = dd1.getProduct(dd2);
        childDdQueue.push(dd3);
      }
      dd = childDdQueue.top();
    }
  }
  // Map<Int,tuple<Float,Float,bool, Int>> ddVarWts;
  Map<Int,tuple<Number,Number,bool, Int>> ddVarWts;
  for (auto& pVar: joinNode->projectionVars){
    Int ddVar = cnfVarToDdVarMap.at(pVar);
    Number posWt = cnf.literalWeights.at(pVar);//.fraction;
    Number negWt = cnf.literalWeights.at(-pVar);//.fraction;
    bool additiveFlag = cnf.outerVars.contains(pVar);
    if (existRandom) {
      additiveFlag = !additiveFlag;
    }
    Int asmt;
    auto it = assignment.find(pVar);
    if (it != assignment.end()) { // literal has assigned value
      asmt = it->second?1:-1;
    }else{ //unassigned
      asmt = 0;
    }
    ddVarWts[ddVar]={posWt,negWt,additiveFlag,asmt};
  }
  dd = dd.getAbstraction(ddVarWts,pmParams.logBound,maximizationStack,pmParams.maximizerFormat,pmParams.substitutionMaximization,verboseSolving);
  if (dd.isZero()){
    printLine("WARNING: Returned Dd after abstraction is zero at joinNode number "+to_string(joinNodesProcessed));
  }
  const Set<Int> sup = dd.getSupport();
  for (auto& cnfVar:joinNode->projectionVars){
    Int ddVar = cnfVarToDdVarMap.at(cnfVar);
    assert (!sup.contains(ddVar));
  }

  joinNodesProcessed ++;
  if (((joinNodesProcessed-1)%(std::max((JoinNode::nodeCount/10),1LL)))==1) printLine(to_string(joinNodesProcessed)+"/"+to_string(JoinNode::nodeCount)+":"+to_string(util::getDuration(executorStartPoint))+" ");
  return dd;
}

Assignment Executor::getMaximizer(Int declaredVarCount) {
  vector<int> ddVarAssignment(ddVarToCnfVarMap.size(), -1); // uses init value -1 (neither 0 nor 1) to test assertion in function Cudd_Eval
  Assignment cnfVarAssignment;

  while (!maximizationStack.empty()) {
    pair<Int, Dd> ddVarAndDsgn = maximizationStack.back();
    Int ddVar = ddVarAndDsgn.first;
    Dd dsgn = ddVarAndDsgn.second;

    bool val = dsgn.evalAssignment(ddVarAssignment);
    ddVarAssignment[ddVar] = val;
    cnfVarAssignment.insert({ddVarToCnfVarMap.at(ddVar), val});

    maximizationStack.pop_back();
  }
  return cnfVarAssignment;
}

void Dpve::setLogBound() {
  if (p.pmParams.logBound > -INF) {} // LOG_BOUND_OPTION
  else if (!p.pmParams.thresholdModel.empty()) { // THRESHOLD_MODEL_OPTION
    p.pmParams.logBound = e->solveSubtree(joinRoot,p.pmParams, Assignment(p.pmParams.thresholdModel)).extractConst().fraction;
    printRow("logBound", p.pmParams.logBound);
  }
  else if (p.pmParams.satSolverPruning) { // SAT_SOLVER_PRUNING
    SatSolver satSolver(p.cnf);
    satSolver.checkSat(true);
    Assignment model = satSolver.getModel();
    p.pmParams.logBound = e->solveSubtree(joinRoot,p.pmParams, model).extractConst().fraction;
    printRow("logBound", p.pmParams.logBound);
    printLine(model.getShortFormat(p.cnf.declaredVarCount));
  }
}

Number Dpve::adjustSolutionToHiddenVar(const Number &apparentSolution, Int cnfVar, const bool additiveFlag) {
  if (p.cnf.apparentVars.contains(cnfVar)) {
    return apparentSolution;
  }

  const Number positiveWeight = p.cnf.literalWeights.at(cnfVar);
  const Number negativeWeight = p.cnf.literalWeights.at(-cnfVar);
  if (additiveFlag) {
    Number s = positiveWeight+negativeWeight;
    return p.logCounting ? (apparentSolution + (positiveWeight + negativeWeight).getLog10()) : (apparentSolution * (positiveWeight + negativeWeight));
  }
  else {
    return p.logCounting ? (apparentSolution + max(positiveWeight, negativeWeight).getLog10()) : (apparentSolution * max(positiveWeight, negativeWeight)); // weights are positive
  }
}

Number Dpve::getAdjustedSolution(const Number &apparentSolution) {
  Number n = apparentSolution;
  for (Int var = 1; var <= p.cnf.declaredVarCount; var++) { // processes inner vars
    if (!p.cnf.outerVars.contains(var)) {
      n = adjustSolutionToHiddenVar(n, var, p.existRandom);
    }
  }
  for (Int var : p.cnf.outerVars) {
    n = adjustSolutionToHiddenVar(n, var, !p.existRandom);
  }
  if (p.scalingFactor == 0){
    //do nothing
  } else{
    if(p.logCounting){
      //convert to log first
      Number sf = Number(p.scalingFactor/log2(10));
      n += sf; //adding since we are in logscale
    } else{
      n = Number::mul_exp2(n,p.scalingFactor); //scalingFactor is implicitly cast to int here
    } 
  }
  return n;
}

Number Dpve::getMaximizerValue(const Assignment& maximizer) {
  Dd dd = e->solveSubtree(joinRoot, p.pmParams, maximizer);
  Number solution = dd.extractConst();
  return getAdjustedSolution(solution);
}

Executor::Executor(const Cnf& cnf, const Map<Int, Int>& cnfVarToDdVarMap, const vector<Int>& ddVarToCnfVarMap,
      const bool existRandom, const string joinPriority, const Int satFilter, const Int verboseSolving, const Int verboseProfiling,
      const Map<Int, vector<Int>> levelMaps_): 
    cnf(cnf),
    cnfVarToDdVarMap(cnfVarToDdVarMap),
    ddVarToCnfVarMap(ddVarToCnfVarMap),
    existRandom(existRandom),
    joinPriority(joinPriority),
    satFilter(satFilter),
    verboseSolving(verboseSolving),
    verboseProfiling(verboseProfiling),
    levelMaps(levelMaps_)
    {
      joinNodesProcessed = 0;
      executorStartPoint = util::getTimePoint();
    }

void Dpve::reorder(){ 
  // condition used by cudds internal sifting : mgr->keys - mgr->isolated; 
}

Dpve::Dpve(const io::InputParams& p_):p(p_){
  //construct join tree
  //compute var order
  Dd::init(p.ddPackage,p.cnf.apparentVars.size(),p.logCounting,p.atomicAbstract, p.weightedCounting, p.multiplePrecision,p.tableRatio,p.initRatio,
    p.threadCount,p.maxMem,p.dynVarOrdering,0);
}

Dpve::~Dpve(){
  Dd::stop();
  if(p.satFilter>0){
    delete s;
  }
  if(p.satFilter!=1){
    delete e;
  }
}

pair<Number, Assignment> Dpve::computeSolution(){
  JoinTreeProcessor::toolStartPoint = p.toolStartPoint;
  JoinTreeProcessor::verboseJoinTree = p.verboseJoinTree;
  JoinTreeProcessor joinTreeProcessor(p.plannerWaitDuration, p.cnf);
  
  Map<Int, Number> unprunableWeights = p.cnf.getUnprunableWeights();
  if (!unprunableWeights.empty() && (p.pmParams.logBound > -INF || !p.pmParams.thresholdModel.empty() || p.pmParams.satSolverPruning)) {
    JoinTreeProcessor::killPlanner();
    printLine();
    printLine("unprunable literal weights:");
    for (const auto& [literal, weight] : unprunableWeights) {
      p.cnf.printLiteralWeight(literal, weight);
    }
    throw util::MyError("must not prune if there are unprunable weights");
  }
  
  TimePoint ddVarOrderStartPoint = util::getTimePoint();
  const JoinNonterminal* joinRoot = joinTreeProcessor.getJoinTreeRoot();
  vector<Int> ddVarToCnfVarMap ;
  ddVarToCnfVarMap = joinRoot->getVarOrder(p.ddVarOrderHeuristic, p.cnf); // e.g. [42, 13], i.e. ddVarOrder
  if (p.verboseSolving >= 1) {
    io::printRow("diagramVarSeconds", util::getDuration(ddVarOrderStartPoint));
  }
  Map<Int, Int> cnfVarToDdVarMap; // e.g. {42: 0, 13: 1}
  for (Int ddVar = 0; ddVar < ddVarToCnfVarMap.size(); ddVar++) {
    Int cnfVar = ddVarToCnfVarMap.at(ddVar);
    cnfVarToDdVarMap[cnfVar] = ddVar;
  }
  if (p.dynVarOrdering == 1){
    printLine("Computing all reordering var orders.."," ");
    TimePoint levelOrdersStartPoint = util::getTimePoint();
    for (auto [v,s] : CNF_VAR_ORDER_HEURISTICS){
      if (v == COLAMD_HEURISTIC || v == RANDOM_HEURISTIC || v == LEX_M_HEURISTIC){
        continue;
      }
      printLine(s," ");
      auto vo = joinRoot->getVarOrder(v, p.cnf); // returns d2cMap. i.e. vo[ddVarIndex] = cnfVarIndex
      //we will interpret vo as a reordering map. i.e. 
      // ddVarIndex will now be ddVarLevel since ddVarIndex for all variables will always be fixed
      // i.e. vo[ddVarLevel] = cnfVarIndex
      // i.e. cnfVarIndex must be put in the level specified by ddVarLevel
      for (auto cnfVar : vo){
        levelMaps[v].push_back(cnfVarToDdVarMap.at(cnfVar)); // so we find ddVarIndex of cnfVar using c2dMap and insert into levelmap.
      }
      //so levelMap[ddVarLevel] = ddVarIndex
      
      //Simiarly for -v
      for (auto it = vo.end(); it--!=vo.begin();){
        Int cnfVar = *it;
        levelMaps[-v].push_back(cnfVarToDdVarMap.at(cnfVar));
      }
      //levelMap will be used as permutation for Cudd_ShuffleHeap(perm)
      //"The i-th entry of the permutation array contains the index of the variable that should be brought to the i-th level."
    }
    io::printRow("allDiagramOrdersSeconds", util::getDuration(levelOrdersStartPoint));
  }
  
  if (p.satFilter>0){
    s = new SatFilter(p.cnf,cnfVarToDdVarMap,ddVarToCnfVarMap);
    printLine("Computing SatFilter ...");
    bool solution = s->solveSubtree(static_cast<const JoinNode*>(joinRoot)).isTrue();
    if (!solution){
      throw util::UnsatException();
    }
    printLine("Done constructing SatFilter. Applying SatFilter...");
    s->filterBdds(joinRoot,Dd::getOneBdd());
    printLine("Done Applying SatFilter!");
    printLine();
  }
  if(p.satFilter!=1){
    printLine("Starting executor...");
    e = new Executor(p.cnf,cnfVarToDdVarMap,ddVarToCnfVarMap,p.existRandom,p.joinPriority,p.satFilter,p.verboseSolving,p.verboseProfiling, levelMaps);
    setLogBound();

    Dd res = e->solveSubtree(static_cast<const JoinNode*>(joinRoot), p.pmParams);
    Number apparentSolution = res.extractConst();

    if (p.pmParams.logBound > -INF) {
      printRow("prunedDiagrams", Dd::prunedDdCount);
      printRow("pruningSeconds", Dd::pruningDuration);
    }

    if (p.verboseSolving >= 1) {
      printRow("apparentSolution", apparentSolution);
    }
    const Number adjustedSolution = getAdjustedSolution(apparentSolution);
    Assignment maximizer;
    if (p.pmParams.maximizerFormat) {
      maximizer = e->getMaximizer(p.cnf.declaredVarCount);
    }
    return(pair<Number,Assignment>(adjustedSolution,maximizer));
  }
  return pair<Number,Assignment>(Number("0"),Assignment());
}