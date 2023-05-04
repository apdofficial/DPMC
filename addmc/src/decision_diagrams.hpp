#pragma once

#include "common.hpp"
#include "formula.hpp"

#include "../libraries/cudd/cplusplus/cuddObj.hh"
#include "../libraries/cudd/cudd/cuddInt.h"
#include "../libraries/cudd/util/util.h"

#include "../libraries/sylvan/src/sylvan.h"
#include "../libraries/sylvan/src/sylvan_gmp.h"
#include "../libraries/sylvan/src/sylvan_obj.hpp"

#include <gmpxx.h>

#include <tuple>

using sylvan::gmp_op_max_CALL;
using sylvan::gmp_op_plus_CALL;
using sylvan::gmp_op_times_CALL;
using sylvan::mtbdd_apply_CALL;
using sylvan::mtbdd_fprintdot_nc;
using sylvan::mtbdd_getdouble;
using sylvan::mtbdd_getvalue;
using sylvan::mtbdd_gmp;
using sylvan::mtbdd_leafcount_more;
using sylvan::mtbdd_makenode;

using sylvan::sylvan_set_limits;
using sylvan::sylvan_init_package;
using sylvan::sylvan_init_mtbdd;

using sylvan::Mtbdd;
using sylvan::MTBDD;

using sylvan::Bdd;
using sylvan::mtbdd_apply_RUN;

using std::tuple;

namespace dpve{

class Dd { // wrapper for CUDD and Sylvan
public:
  ADD cuadd; // CUDD
  Mtbdd mtbdd; // Sylvan

  BDD cubdd; 
  Bdd sybdd;

  static size_t maxDdLeafCount;
  static size_t maxDdNodeCount;

  static size_t prunedDdCount;
  static Float pruningDuration;
 
  static bool dynOrderEnabled;

  size_t getLeafCount() const;
  size_t getNodeCount() const;

  Dd(const ADD& cuadd);
  Dd(const Mtbdd& mtbdd);
  Dd(const Dd& dd);
  Dd(const BDD& cubdd);
  Dd(const Bdd& sybdd);


  Number extractConst() const; // does not read logCounting
  static Dd getConstDd(const Number& n); // reads logCounting
  static Dd getZeroDd(); // returns minus infinity if logCounting
  static Dd getOneDd(); // returns zero if logCounting
  static Dd getVarDd(Int ddVar, bool val);
  bool operator!=(const Dd& rightDd) const;
  Dd getComposition(Int ddVar, bool val) const; // restricts *this to ddVar=val
  Dd getProduct(const Dd& dd) const; // reads logCounting
  Dd getSum(const Dd& dd) const; // reads logCounting
  Dd getMax(const Dd& dd) const; // real max (not 0-1 max)
  Dd getXor(const Dd& dd) const; // must be 0-1 DDs
  Set<Int> getSupport() const;
  Dd getBoolDiff(const Dd& rightDd) const; // returns 0-1 DD for *this >= rightDd
  bool evalAssignment(vector<int>& ddVarAssignment) const;
  //getAbstraction is not a const method
  Dd getAbstraction(Map<Int,tuple<Float,Float,bool,Int>> ddVarWts, Float logBound, vector<pair<Int, Dd>>& maximizationStack, bool maximizerFormat, bool substitutionMaximization, Int verboseSolving);
  Dd getPrunedDd(Float lowerBound) const;
  void writeDotFile(const string& dotFileDir = "./") const;
  static void writeInfoFile(const string& filePath);
  
  static bool enableDynamicOrdering();
  static bool disableDynamicOrdering();
  static int postReorderHook(DdManager* dd, const char *str, void *data);
  static int preReorderHook(DdManager* dd, const char *str, void *data);
  static int postGCHook(DdManager *dd, const char *str, void *data);
  static int preGCHook(DdManager *dd, const char *str, void *data);
  
  static Dd getZeroBdd(); // returns minus infinity if logCounting
  static Dd getOneBdd(); // returns zero if logCounting
  static Dd getVarBdd(Int ddVar, bool val);

  Dd getBddExists(
    vector<Int> ddVars,
    const vector<Int>& ddVarToCnfVarMap
  ) const;
  Dd getBddAnd(const Dd& dd) const;
  Dd getBddOr(const Dd& dd) const; // must be Bdd
  bool isTrue() const; // must be Bdd
  Set<Int> getBddSupport() const;
  Dd getFilteredBdd(const Dd);
  Dd getAdd();

  static Dd getClauseDd(Map<Int, pair<Int,Int>> clauseDDVarSignAndAsmts, bool xorFlag);

  static void init(string ddPackage_, bool logCounting_, bool atomicAbstract_=1, bool weightedCounting_=0, bool multiplePrecision_=0, Int tableRatio=0, Int initRatio=0, Int threadCount=1, Float maxMem=0, bool dynVarOrdering=0, Int dotFileIndex_=0);
  static Float getNegWt(Int ddVar);
  static void stop();
  private:
    static string ddPackage;
    static Cudd* mgr;
    static bool logCounting;
    static bool atomicAbstract;
    static bool weightedCounting;
    static bool multiplePrecision;
    static Int dotFileIndex;
    static Map<Int, pair<Float,Float>> wtMap; 
};
} //end namespace dpve