#include "decision_diagrams.hpp"
#include "io.hpp"
#include "util.hpp"

using dpve::Dd;
using dpve::Float;
using dpve::Set;
using dpve::Int;
using dpve::Number;
using dpve::io::printLine;
using std::to_string;
using dpve::Map;
/* class Dd ================================================================= */

size_t Dd::maxDdLeafCount;
size_t Dd::maxDdNodeCount;

size_t Dd::prunedDdCount;
Float Dd::pruningDuration;

bool Dd::dynOrderEnabled = false;

string Dd::ddPackage;
Cudd* Dd::mgr = 0;
bool Dd::logCounting = 0;
bool Dd::atomicAbstract = 0;
bool Dd::weightedCounting = 0;
bool Dd::multiplePrecision = 0;
Int Dd::dotFileIndex = 0;
Map<Int, pair<Float,Float>> Dd::wtMap;

bool Dd::enableDynamicOrdering() {
  assert(ddPackage == CUDD_PACKAGE);
  mgr->AutodynEnable();
  Dd::dynOrderEnabled = true;
  mgr->AddHook(postReorderHook, CUDD_POST_REORDERING_HOOK);
  mgr->AddHook(preReorderHook, CUDD_PRE_REORDERING_HOOK);
  return (Dd::dynOrderEnabled);
}

bool Dd::disableDynamicOrdering() {
  assert(ddPackage == CUDD_PACKAGE);
  mgr->AutodynDisable();
  Dd::dynOrderEnabled = false;
  return (Dd::dynOrderEnabled);
}

int Dd::postGCHook(DdManager *dd, const char *str, void *data){
  printLine("GC done!");
  return 1;
}

int Dd::preGCHook(DdManager *dd, const char *str, void *data){
  printLine("Starting GC..","   ");
  return 1;
}

int Dd::postReorderHook(DdManager *dd, const char *str, void *data){
    unsigned long initialTime = (unsigned long) (ptruint) data;
    int retval;
    unsigned long finalTime = util_cpu_time();
    double totalTimeSec = (double)(finalTime - initialTime) / 1000.0;

    retval = fprintf(dd->out,"%ld nodes in %g sec\n", strcmp(str, "BDD") == 0 ? Cudd_ReadNodeCount(dd) : Cudd_zddReadNodeCount(dd), totalTimeSec);
    if (retval == EOF) return(0);
    retval = fflush(dd->out);
    if (retval == EOF) return(0);
    return(1);

}

int Dd::preReorderHook(DdManager *dd, const char *str, void *data){
    Cudd_ReorderingType method = (Cudd_ReorderingType) (ptruint) data;
    int retval;

    retval = fprintf(dd->out,"%s reordering with ", str);
    if (retval == EOF) return(0);
    switch (method) {
    case CUDD_REORDER_SIFT_CONVERGE:
    case CUDD_REORDER_SYMM_SIFT_CONV:
    case CUDD_REORDER_GROUP_SIFT_CONV:
    case CUDD_REORDER_WINDOW2_CONV:
    case CUDD_REORDER_WINDOW3_CONV:
    case CUDD_REORDER_WINDOW4_CONV:
    case CUDD_REORDER_LINEAR_CONVERGE:
        retval = fprintf(dd->out,"converging ");
        if (retval == EOF) return(0);
        break;
    default:
        break;
    }
    switch (method) {
    case CUDD_REORDER_RANDOM:
    case CUDD_REORDER_RANDOM_PIVOT:
        retval = fprintf(dd->out,"random");
        break;
    case CUDD_REORDER_SIFT:
    case CUDD_REORDER_SIFT_CONVERGE:
        retval = fprintf(dd->out,"sifting");
        break;
    case CUDD_REORDER_SYMM_SIFT:
    case CUDD_REORDER_SYMM_SIFT_CONV:
        retval = fprintf(dd->out,"symmetric sifting");
        break;
    case CUDD_REORDER_LAZY_SIFT:
        retval = fprintf(dd->out,"lazy sifting");
        break;
    case CUDD_REORDER_GROUP_SIFT:
    case CUDD_REORDER_GROUP_SIFT_CONV:
        retval = fprintf(dd->out,"group sifting");
        break;
    case CUDD_REORDER_WINDOW2:
    case CUDD_REORDER_WINDOW3:
    case CUDD_REORDER_WINDOW4:
    case CUDD_REORDER_WINDOW2_CONV:
    case CUDD_REORDER_WINDOW3_CONV:
    case CUDD_REORDER_WINDOW4_CONV:
        retval = fprintf(dd->out,"window");
        break;
    case CUDD_REORDER_ANNEALING:
        retval = fprintf(dd->out,"annealing");
        break;
    case CUDD_REORDER_GENETIC:
        retval = fprintf(dd->out,"genetic");
        break;
    case CUDD_REORDER_LINEAR:
    case CUDD_REORDER_LINEAR_CONVERGE:
        retval = fprintf(dd->out,"linear sifting");
        break;
    case CUDD_REORDER_EXACT:
        retval = fprintf(dd->out,"exact");
        break;
    default:
        return(0);
    }
    if (retval == EOF) return(0);

    retval = fprintf(dd->out,": from %ld to ... ", strcmp(str, "BDD") == 0 ? Cudd_ReadNodeCount(dd) : Cudd_zddReadNodeCount(dd));
    if (retval == EOF) return(0);
    fflush(dd->out);
    return(1);

} /* end of Cudd_StdPreReordHook */


size_t Dd::getLeafCount() const {
  if (ddPackage == CUDD_PACKAGE) {
    return cuadd.CountLeaves();
  }
  MTBDD d = mtbdd.GetMTBDD();
  return mtbdd_leafcount(d);
}

size_t Dd::getNodeCount() const {
  if (ddPackage == CUDD_PACKAGE) {
    if (cuadd.getNode()==0){
      return cubdd.nodeCount();
    }
    return cuadd.nodeCount();
  }
  return mtbdd.NodeCount();
}

Dd::Dd(const ADD& cuadd) {
  assert(ddPackage == CUDD_PACKAGE);
  this->cuadd = cuadd;
}

Dd::Dd(const Mtbdd& mtbdd) {
  assert(ddPackage == SYLVAN_PACKAGE);
  this->mtbdd = mtbdd;
}

Dd::Dd(const Dd& dd) {
  if (ddPackage == CUDD_PACKAGE) {
    // *this = Dd(dd.cuadd);
    this->cuadd = dd.cuadd;
    this->cubdd = dd.cubdd;
  }
  else {
    // *this = Dd(dd.mtbdd);
    this->mtbdd = dd.mtbdd;
    this->sybdd = dd.sybdd;
  }

  // maxDdLeafCount = max(maxDdLeafCount, getLeafCount());
  // maxDdNodeCount = max(maxDdNodeCount, getNodeCount());
}

Number Dd::extractConst() const {
  if (ddPackage == CUDD_PACKAGE) {
    ADD minTerminal = cuadd.FindMin();
    assert(minTerminal == cuadd.FindMax());
    return Number(cuddV(minTerminal.getNode()));
  }
  assert(mtbdd.isLeaf());
  if (multiplePrecision) {
    uint64_t val = mtbdd_getvalue(mtbdd.GetMTBDD());
    return Number(mpq_class(reinterpret_cast<mpq_ptr>(val)));
  }
  return Number(mtbdd_getdouble(mtbdd.GetMTBDD()));
}

Dd Dd::getConstDd(const Number& n) {
  if (ddPackage == CUDD_PACKAGE) {
    return logCounting ? Dd(mgr->constant(n.getLog10())) : Dd(mgr->constant(n.fraction));
  }
  if (multiplePrecision) {
    mpq_t q; // C interface
    mpq_init(q);
    mpq_set(q, n.quotient.get_mpq_t());
    Dd dd(Mtbdd(mtbdd_gmp(q)));
    mpq_clear(q);
    return dd;
  }
  return logCounting ? Dd(Mtbdd::doubleTerminal(n.getLog10())) : Dd(Mtbdd::doubleTerminal(n.fraction));
}

Dd Dd::getZeroDd() {
  return getConstDd(Number());
}

Dd Dd::getOneDd() {
  return getConstDd(Number("1"));
}

Dd Dd::getVarDd(Int ddVar, bool val) {
  if (ddPackage == CUDD_PACKAGE) {
    if (logCounting) {
      return Dd(mgr->addLogVar(ddVar, val));
    }
    ADD d = mgr->addVar(ddVar);
    return val ? Dd(d) : Dd(d.Cmpl());
  }
  MTBDD d0 = getZeroDd().mtbdd.GetMTBDD(); // handles logcounting case
  MTBDD d1= getOneDd().mtbdd.GetMTBDD(); // handles logcounting case
  // if (logCounting){
  //   return val? Dd(Mtbdd(sylvan::mtbdd_makenode(ddVar, d0, d1)).Log()) : Dd(Mtbdd(sylvan::mtbdd_makenode(ddVar, d1, d0)).Log());
  // }
  return val ? Dd(Mtbdd(sylvan::mtbdd_makenode(ddVar, d0, d1))) : Dd(Mtbdd(sylvan::mtbdd_makenode(ddVar, d1, d0)));
}

Dd::Dd(const BDD& cubdd){
  assert(ddPackage == CUDD_PACKAGE);
  this->cubdd = cubdd;
}

Dd::Dd(const Bdd& sybdd) {
  assert(ddPackage == SYLVAN_PACKAGE);
  this->sybdd = sybdd;
}

Dd Dd::getZeroBdd() {
  if (ddPackage == CUDD_PACKAGE){
    return mgr->bddZero();
  } else{
    return Bdd::bddZero();
  }
}

Dd Dd::getOneBdd() {
  if (ddPackage == CUDD_PACKAGE){
    return mgr->bddOne();
  } else{
    return Bdd::bddOne();
  }
}

Dd Dd::getVarBdd(Int ddVar, bool val) {
  if (ddPackage == CUDD_PACKAGE) {
    BDD d = mgr->bddVar(ddVar);
    return val ? Dd(d) : Dd(!d);
  } else {
    return val? Bdd::bddVar(ddVar) : !Bdd::bddVar(ddVar); 
  }
}

Dd Dd::getBddAnd(const Dd& dd) const {
  if (ddPackage == CUDD_PACKAGE) {
    return Dd(cubdd * dd.cubdd);
  } else{
    return Dd(sybdd * dd.sybdd);
  }
}

Dd Dd::getBddExists(vector<Int> ddVars, const vector<Int>& ddVarToCnfVarMap) const {
  Dd cube = getOneBdd();
  for (auto& var : ddVars){
    cube = cube.getBddAnd(getVarBdd(var,true));
  }
  if (ddPackage == CUDD_PACKAGE){
    return cubdd.ExistAbstract(cube.cubdd);
  } else{
    return sybdd.ExistAbstract(sylvan::BddSet(cube.sybdd));
  }
}

Dd Dd::getBddOr(const Dd& dd) const {
  if(ddPackage == CUDD_PACKAGE){
    return Dd(cubdd.Or(dd.cubdd));
  } else{
    return Dd(sybdd.Or(dd.sybdd));
  }
}

bool Dd::isTrue() const{
  if(ddPackage == CUDD_PACKAGE){
     return cubdd.IsOne();
  } else{
    return sybdd.isOne();
  }
}

Set<Int> Dd::getBddSupport() const{
  Set<Int> support;
  if (ddPackage == CUDD_PACKAGE) {
    for (Int ddVar : cubdd.SupportIndices()) {
      support.insert(ddVar);
    }
   
  }
  else {
    Mtbdd cube = sybdd.Support(); // conjunction of all vars appearing in mtbdd
    while (!cube.isOne()) {
      support.insert(cube.TopVar());
      cube = cube.Then();
    }
  }
  return support;
}

Dd Dd::getFilteredBdd(const Dd other){
  Set<Int> otherSupport = other.getBddSupport();
  for (const auto& elem : getBddSupport()) {
    otherSupport.erase(elem);
  }

  Dd cube = getOneBdd();
  for (auto& var : otherSupport){
    cube = cube.getBddAnd(getVarBdd(var,true));
  }
  if (ddPackage == CUDD_PACKAGE){
    return Dd(cubdd.AndAbstract(other.cubdd, cube.cubdd));
  } else{
    return Dd(sybdd.AndAbstract(other.sybdd,cube.sybdd));
  }
}

Dd Dd::getAdd(){
  if (ddPackage == CUDD_PACKAGE){
    return logCounting? Dd(cubdd.Add().Log()) : Dd(cubdd.Add());
  } else{
    return logCounting? Dd(Mtbdd(sybdd).Log()): Dd(Mtbdd(sybdd));
  }
}

/*
First parameter has 3 elements: 
  the ddVar for each cnfvar in the clause.
  the sign/polarity (+1 or 0) for the ddVar
  the assignment to the variable (0: unassigned, +1: positive asnmt -1: negative asnmt)
*/
Dd Dd::getClauseDd(Map<Int, pair<Int,Int>> clauseDDVarSignAndAsmts, bool xorFlag) {
  Dd clauseDd = Dd::getZeroDd();
  for (auto [ddVar , valuePair] : clauseDDVarSignAndAsmts) {
    assert(ddVar>=0);
    auto [sign, asmt] = valuePair;
    assert(sign==1 || sign==0);
    assert(asmt==1 || asmt==0 || asmt==-1);
    if (asmt == 0){ //variable is unassigned
      Dd literalDd = Dd::getVarDd(ddVar, sign);
      clauseDd = xorFlag ? clauseDd.getXor(literalDd) : clauseDd.getMax(literalDd);
    }else if ((sign>0) == (asmt>0)){ //literal is satisfied by assignment
      if (xorFlag) { // flips polarity
        clauseDd = clauseDd.getXor(Dd::getOneDd());
      }
      else { // returns satisfied disjunctive clause
        return Dd::getOneDd();
      }
    } // excludes unsatisfied literal from clause otherwise
  }
  return clauseDd;
}

template<typename TReal>
static bool isApproximatelyEqual(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
{
    TReal diff = std::fabs(a - b);
    if (diff <= tolerance)
        return true;

    if (diff < std::fmax(std::fabs(a), std::fabs(b)) * tolerance)
        return true;

    return false;
}

/*
Does not convert weights to logscale regardless of status of logCounting flag
Callers responsibility to convert as necessary.
*/
Float Dd::getNegWt(Int ddVar){
  assert (ddVar>=0);
  auto [posWt, negWt] = wtMap[ddVar];
  if (negWt == 1){
    assert (posWt == 1);
  } else{
    if (!isApproximatelyEqual(negWt + posWt,1.0l,0.001l)){
      printLine("literal weights not summing to 1!! :"+to_string(ddVar)+"  "+to_string(posWt)+":"+to_string(negWt));
      exit(1);
    }
  }
  // return logCounting ? log10l(negWt) : negWt;
  return negWt;
}

Dd Dd::getAbstraction(Map<Int,tuple<Float,Float,bool,Int>> ddVarWts, Float logBound, vector<pair<Int, Dd>>& maximizationStack, bool maximizerFormat, bool substitutionMaximization, Int verboseSolving){
  if (atomicAbstract){
    if (ddPackage == CUDD_PACKAGE){
      ADD wCube = mgr->addOne();
      wtMap.clear();
      for (auto ddVarWt : ddVarWts) {
        Int ddVar = ddVarWt.first;
        const auto [posWt, negWt, additiveFlag, asmt] = ddVarWt.second;
        wtMap[ddVar] = {posWt,negWt};
        assert(asmt==0); //cases with assignments dont currently work with atomic abstract
        ADD v = mgr->addVar(ddVar);
        wCube *= v;
      }
      if (weightedCounting){
        return logCounting? cuadd.WeightedLogSumExistAbstract(wCube, getNegWt) : cuadd.WeightedExistAbstract(wCube, getNegWt);
      } else {
        return logCounting? cuadd.LogSumExistAbstract(wCube) : cuadd.ExistAbstract(wCube);
      }
    } else{
      assert(!weightedCounting);
      sylvan::BddSet b;
      for (auto ddVarWt : ddVarWts) {
        Int ddVar = ddVarWt.first;
        b.add((uint32_t)ddVar);
      }
      return logCounting? mtbdd.AbstractLogSumExp(b) : mtbdd.AbstractPlus(b);
    }
  } else{
    Dd dd = *this;
    for (auto ddVarWt : ddVarWts) {
      Int ddVar = ddVarWt.first;
      const auto [posWt, negWt, additiveFlag, asmt] = ddVarWt.second;
      if (asmt != 0) { //variable has an assignment
        dd = dd.getProduct(asmt>0 ? getConstDd(posWt) : getConstDd(negWt));
      } else{
        Dd highTerm = dd.getComposition(ddVar, true).getProduct(getConstDd(posWt));
        Dd lowTerm = dd.getComposition(ddVar, false).getProduct(getConstDd(negWt));

        if (maximizerFormat && !additiveFlag) {
          Dd dsgn = highTerm.getBoolDiff(lowTerm); // derivative sign
          maximizationStack.push_back({ddVar, dsgn});
          if (substitutionMaximization) {
            dd = Dd(cuadd.Compose(dsgn.cuadd, ddVar));
          }
        }
        if (!substitutionMaximization){
          dd = additiveFlag ? highTerm.getSum(lowTerm) : highTerm.getMax(lowTerm);
        }
      }
      if (logBound > -INF) {
        if (posWt != 1 || negWt != 1) {
          Dd prunedDd = dd.getPrunedDd(logBound);
          if (prunedDd != dd) {
            if (verboseSolving >= 3) {
              printLine("writing pre-pruning decision diagram...");
              dd.writeDotFile();

              printLine("writing post-pruning decision diagram...");
              prunedDd.writeDotFile();
            }
            dd = prunedDd;
          }
        }
      }
      //set *this cuadd/mtbdd to One so that the underlying dd can potentially be garbage collected.
      //if this is not done, this cuadd/mtbdd wont be gc'ed until return from getAbstraction
      if (ddPackage==CUDD_PACKAGE){
        this->cuadd = mgr->addOne();
      } else{
        mtbdd = Mtbdd::mtbddOne();
      }
    }
    return dd;
  }
}

bool Dd::operator!=(const Dd& rightDd) const {
  if (ddPackage == CUDD_PACKAGE) {
    return cuadd != rightDd.cuadd;
  }
  return mtbdd != rightDd.mtbdd;
}

Dd Dd::getComposition(Int ddVar, bool val) const {
  if (ddPackage == CUDD_PACKAGE) {
    if (dpve::util::isFound(ddVar, cuadd.SupportIndices())) {
      return Dd(cuadd.Compose(val ? mgr->addOne() : mgr->addZero(), ddVar));
    }
    return *this;
  }
  sylvan::MtbddMap m;
  m.put(ddVar, val ? Mtbdd::mtbddOne() : Mtbdd::mtbddZero());
  return Dd(mtbdd.Compose(m));
}

Dd Dd::getProduct(const Dd& dd) const {
  if (ddPackage == CUDD_PACKAGE) {
    return logCounting ? Dd(cuadd + dd.cuadd) : Dd(cuadd * dd.cuadd);
  }
  if (multiplePrecision) {
    // LACE_ME;
    return Dd(Mtbdd(gmp_times(mtbdd.GetMTBDD(), dd.mtbdd.GetMTBDD())));
  }
  return logCounting? Dd(mtbdd + dd.mtbdd) : Dd(mtbdd * dd.mtbdd);
}

Dd Dd::getSum(const Dd& dd) const {
  if (ddPackage == CUDD_PACKAGE) {
    return logCounting ? Dd(cuadd.LogSumExp(dd.cuadd)) : Dd(cuadd + dd.cuadd);
  }
  if (multiplePrecision) {
    // LACE_ME;
    return Dd(Mtbdd(gmp_plus(mtbdd.GetMTBDD(), dd.mtbdd.GetMTBDD())));
  }
  return logCounting? Dd(mtbdd.LogSumExp(dd.mtbdd)) : Dd(mtbdd + dd.mtbdd);
}

Dd Dd::getMax(const Dd& dd) const {
  if (ddPackage == CUDD_PACKAGE) {
    return Dd(cuadd.Maximum(dd.cuadd));
  }
  if (multiplePrecision) {
    // LACE_ME;
    return Dd(Mtbdd(gmp_max(mtbdd.GetMTBDD(), dd.mtbdd.GetMTBDD())));
  }
  return Dd(mtbdd.Max(dd.mtbdd));
}

Dd Dd::getXor(const Dd& dd) const {
  assert(ddPackage == CUDD_PACKAGE);
  return logCounting ? Dd(cuadd.LogXor(dd.cuadd)) : Dd(cuadd.Xor(dd.cuadd));
}

Set<Int> Dd::getSupport() const {
  Set<Int> support;
  if (ddPackage == CUDD_PACKAGE) {
    for (Int ddVar : cuadd.SupportIndices()) {
      support.insert(ddVar);
    }
  }
  else {
    Mtbdd cube = mtbdd.Support(); // conjunction of all vars appearing in mtbdd
    while (!cube.isOne()) {
      support.insert(cube.TopVar());
      cube = cube.Then();
    }
  }
  return support;
}

Dd Dd::getBoolDiff(const Dd& rightDd) const {
  assert(ddPackage == CUDD_PACKAGE);
  return Dd((cuadd - rightDd.cuadd).BddThreshold(0).Add());
}

bool Dd::evalAssignment(vector<int>& ddVarAssignment) const {
  assert(ddPackage == CUDD_PACKAGE);
  Number n = Dd(cuadd.Eval(&ddVarAssignment.front())).extractConst();
  return n == Number(1);
}


Dd Dd::getPrunedDd(Float lowerBound) const {
  assert(logCounting);

  TimePoint pruningStartPoint = util::getTimePoint();

  ADD bound = mgr->constant(lowerBound);
  ADD prunedDd = cuadd.LogThreshold(bound);

  pruningDuration += util::getDuration(pruningStartPoint);

  if (prunedDd != cuadd) {
    prunedDdCount++;
  }

  return Dd(prunedDd);
}

void Dd::writeDotFile(const string& dotFileDir) const {
  string filePath = dotFileDir + "dd" + to_string(dotFileIndex++) + ".dot";
  FILE* file = fopen(filePath.c_str(), "wb"); // writes to binary file

  if (ddPackage == CUDD_PACKAGE) { // davidkebo.com/cudd#cudd6
    DdNode** ddNodeArray = static_cast<DdNode**>(malloc(sizeof(DdNode*)));
    ddNodeArray[0] = cuadd.getNode();
    Cudd_DumpDot(mgr->getManager(), 1, ddNodeArray, NULL, NULL, file);
    free(ddNodeArray);
  }
  else {
    mtbdd_fprintdot_nc(file, mtbdd.GetMTBDD());
  }

  fclose(file);
  printLine("wrote decision diagram to file "+filePath);
}

void Dd::writeInfoFile(const string& filePath) {
  assert(ddPackage == CUDD_PACKAGE);
  FILE* file = fopen(filePath.c_str(), "w");
  Cudd_PrintInfo(mgr->getManager(), file);
  fclose(file);
  printLine("wrote CUDD info to file "+filePath);
}


void Dd::init(string ddPackage_, bool logCounting_, bool atomicAbstract_, bool weightedCounting_, bool multiplePrecision_, Int tableRatio, Int initRatio, Int threadCount, Float maxMem, bool dynVarOrdering, Int dotFileIndex_){
  ddPackage = ddPackage_;
  logCounting = logCounting_;
  atomicAbstract = atomicAbstract_;
  weightedCounting = weightedCounting_;
  multiplePrecision = multiplePrecision_;
  dotFileIndex = dotFileIndex_;
  if (ddPackage == CUDD_PACKAGE){
    mgr = new Cudd(
    0, // init num of BDD vars
    0, // init num of ZDD vars
    CUDD_UNIQUE_SLOTS, // init num of unique-table slots; cudd.h: #define CUDD_UNIQUE_SLOTS 256
    CUDD_CACHE_SLOTS, // init num of cache-table slots; cudd.h: #define CUDD_CACHE_SLOTS 262144
    maxMem * MEGA // maxMemory
    );
    if (dynVarOrdering){
      enableDynamicOrdering();
    }
    mgr->AddHook(postGCHook, CUDD_POST_GC_HOOK);
    mgr->AddHook(preGCHook, CUDD_PRE_GC_HOOK);
    printLine("CUDD Max Mem: "+to_string(mgr->ReadMaxMemory()));
  } else{
    lace_start(threadCount, 1000000); // auto-detect number of workers, use a 1,000,000 size task queue
    // Init Sylvan  
    sylvan_set_limits(maxMem*MEGA, tableRatio, initRatio);
    sylvan_init_package();
    sylvan_init_mtbdd();
    if (multiplePrecision) {
      sylvan::gmp_init();
    }
  }
}

void Dd::stop(){
  if (ddPackage == SYLVAN_PACKAGE) { // quits Sylvan
    sylvan::sylvan_quit();
    lace_stop();
  }else{
    printLine("Current CUDD used Mem: "+to_string(mgr->ReadMemoryInUse()));
    printLine("Total GC time: "+to_string(mgr->ReadGarbageCollectionTime()));
    mgr->info();
    Cudd_Quit(mgr->getManager());
  }
}