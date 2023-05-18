#pragma once

/* inclusions =============================================================== */
#include "common.hpp"
#include "decision_diagrams.hpp"
#include "formula.hpp"
#include "io.hpp"
#include "jointrees.hpp"
#include "sat_solver.hpp"

using dpve::io::PruneMaxParams;

namespace dpve{
class Executor {
  public:
    Dd solveSubtree(const JoinNode* joinNode, const PruneMaxParams& pmParams, const Assignment& assignment = Assignment());
    Assignment getMaximizer(Int declaredVarCount);
    Executor(const Cnf& cnf, const Map<Int, Int>& cnfVarToDdVarMap, const vector<Int>& ddVarToCnfVarMap, const bool existRandom, 
      const string joinPriority, const Int satFilter, const Int verboseSolving, const Int verboseProfiling,
      const Map<Int, vector<Int>> levelMaps_ = Map<Int, vector<Int>>());
    
    Float reOrdThresh = 0.7;

  private:
    const Cnf& cnf;
    const Map<Int, Int>& cnfVarToDdVarMap;
    const vector<Int>& ddVarToCnfVarMap;
    vector<pair<Int, Dd>> maximizationStack; // pair<DD var, derivative sign>
    
    const Map<Int, vector<Int>>& levelMaps;

    const bool existRandom;
    const string joinPriority; 
    const Int satFilter; 

    const Int verboseSolving;
    const Int verboseProfiling;
    
    TimePoint executorStartPoint;
    Int joinNodesProcessed=0;
    Map<Int, Float> varDurations; // CNF var |-> total execution time in seconds
    Map<Int, size_t> varDdSizes; // CNF var |-> max DD size
};

class SatFilter{
  public:
    Dd solveSubtree(const JoinNode* joinNode);
    bool filterBdds(const JoinNode* joinNode,const Dd parentBdd);
    SatFilter(const Cnf& cnf, const Map<Int, Int>& cnfVarToDdVarMap,const vector<Int>& ddVarToCnfVarMap);
  private:
  const Map<Int, Int>& cnfVarToDdVarMap;
  const vector<Int>& ddVarToCnfVarMap;
  const Cnf& cnf;

  TimePoint satFilterStartPoint;
  Int joinNodesProcessed=0;

  string joinPriority;

  Dd getClauseBdd(const Clause& clause);
  // recursively computes valuation of join tree node
};

class Dpve{
  private:
    Executor *e;
    SatFilter *s;
    const JoinNonterminal* joinRoot;
    const io::InputParams& p;    
    Map<Int, Int> cnfVarToDdVarMap;
    const vector<Int> ddVarToCnfVarMap;
    Map<Int,vector<Int>> levelMaps;

    void setLogBound();

    Number adjustSolutionToHiddenVar(const Number &apparentSolution, Int cnfVar, const bool additiveFlag);
    Number getAdjustedSolution(const Number &apparentSolution);
    void reorder();
  public:
    Dpve(const io::InputParams& p);
    pair<Number, Assignment> computeSolution();
    Number getMaximizerValue(const Assignment& maximizer);
    ~Dpve();
};

} //end namepace dpve
