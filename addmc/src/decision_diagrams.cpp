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
Float Dd::reordThresh, Dd::reordThreshInc;
Int Dd::maxSwaps, Dd::maxSwapsInc, Dd::swapTime, Dd::dynVarOrdering = 0, Dd::lut;
bool Dd::didReordering, Dd::noReordSinceGC; 

string Dd::ddPackage;
Cudd* Dd::mgr = 0;
bool Dd::logCounting = 0;
bool Dd::atomicAbstract = 0;
bool Dd::weightedCounting = 0;
bool Dd::multiplePrecision = 0;
Int Dd::dotFileIndex = 0;
Map<Int, pair<Number,Number>> Dd::wtMap;

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
  Int memused = mgr->ReadMemoryInUse();
  Int memMax = mgr->ReadMaxMemory();
  printLine("GC done! memused: "+to_string(memused)+" / "+to_string(memMax)+" = "+to_string((memused+0.0)/memMax));
  noReordSinceGC = true;
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
    mpq_ptr op = (mpq_ptr)val;
    mpq_t mres;
    mpq_init(mres);
    mpq_set(mres,op);
    // printLine(" "," ");
    // mpq_out_str(stdout,10,mres);
    // printLine();
    Number n = Number(mpq_class(mres));
    mpq_clear(mres);
    return n;
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
  manualReorder();
  if (ddPackage == CUDD_PACKAGE) {
    return Dd(cubdd * dd.cubdd);
  } else{
    return Dd(sybdd * dd.sybdd);
  }
}

Dd Dd::getBddExists(vector<Int> ddVars, const vector<Int>& ddVarToCnfVarMap) const {
  manualReorder();
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
    if (multiplePrecision){
      Mtbdd temp = Mtbdd(gmp_convertToGMP(sybdd.GetBDD()));
      return Dd(temp);
    }
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
  assert(!multiplePrecision);
  auto [posWtN, negWtN] = wtMap[ddVar];
  Float posWt = posWtN.fraction, negWt = negWtN.fraction;
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

Dd Dd::getAbstraction(Map<Int,tuple<Number,Number,bool,Int>> ddVarWts, Float logBound, vector<pair<Int, Dd>>& maximizationStack, bool maximizerFormat, bool substitutionMaximization, Int verboseSolving){
  manualReorder();
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
      if (multiplePrecision){
        Mtbdd temp = Mtbdd::mtbddOne();
        for (auto ddVarWt : ddVarWts) {
          Int ddVar = ddVarWt.first;
          assert(ddVar>=0);
          temp = temp.Times(Mtbdd::mtbddVar(ddVar));
        }
        return Dd(Mtbdd(gmp_abstract_plus(mtbdd.GetMTBDD(),temp.GetMTBDD())));
      } else{
        sylvan::BddSet b;
        for (auto ddVarWt : ddVarWts) {
          Int ddVar = ddVarWt.first;
          b.add((uint32_t)ddVar);
        }
        return logCounting? mtbdd.AbstractLogSumExp(b) : mtbdd.AbstractPlus(b);
      }
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
  manualReorder();
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

bool Dd::isZero() const{
  if (ddPackage == CUDD_PACKAGE){
    return logCounting? cuadd.getNode() == mgr->minusInfinity().getNode() : cuadd.IsZero();
  } else{
    if (multiplePrecision){
      mpq_t z;
      mpq_init(z);
      mpq_set_d(z,0);
      MTBDD mz = mtbdd_gmp(z);
      bool isZero = (mz == mtbdd.GetMTBDD());
      mpq_clear(z);
      return isZero;
    } else{
      // return logCounting? getConstDd(Number()).mtbdd.GetMTBDD() == mtbdd.GetMTBDD() : mtbdd.isZero();
      // if (mtbdd.isLeaf()){
      //   printLine("Found leaf value "+to_string(mtbdd_getdouble(mtbdd.GetMTBDD())));
      // }
      return mtbdd_getdouble(mtbdd.GetMTBDD()) == (logCounting? Number().getLog10() : 0);
    }
  }
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

#ifdef SYLDEV
static size_t prev_size = 0;
static int terminate_reordering = 0;

VOID_TASK_0(reordering_start)
{
    sylvan::sylvan_gc_RUN();
    size_t size = llmsset_count_marked(sylvan::nodes);
    // std::cout << "\nSylvan:RE: start:    " << size << " size\n";
    printLine("Sylvan size before reordering: "+to_string(size));
}

VOID_TASK_0(reordering_progress)
{
    size_t size = llmsset_count_marked(sylvan::nodes);
    // we need at least 40% reduction in size to continue
    if (size >= prev_size * 0.60) terminate_reordering = 1;
    else prev_size = size;
    // std::cout << "Sylvan:RE: progress: " << size << " size\n";
    printLine("Sylvan reordering progress size: "+to_string(size));
}

VOID_TASK_0(reordering_end)
{
    sylvan::sylvan_gc_RUN();
    size_t size = llmsset_count_marked(sylvan::nodes);
    // std::cout << "Sylvan:RE: end:      " << size << " size\n";
    printLine("Sylvan size after reordering: "+to_string(size));
}

int should_reordering_terminate()
{
    return terminate_reordering;
}
#endif

VOID_TASK_0(gc_start)
{
   size_t used, total;
   sylvan::sylvan_table_usage_RUN(&used, &total);
  //  cout << "\nSylvan:GC: start:" << used << "/" << total << "\n";
  printLine("Sylvan before GC. Used : "+to_string(used)+" out of "+to_string(total),"    ");
}

VOID_TASK_0(gc_end)
{
   size_t used, total;
   sylvan::sylvan_table_usage_RUN(&used, &total);
  //  cout << "Sylvan:GC: end:" << used << "/" << total << "\n";
  Dd::noReordSinceGC = true;
  printLine("Sylvan after GC. Used : "+to_string(used)+" out of "+to_string(total));
}

bool Dd::beforeReorder(){
  //check if growth has finished
  //if cant grow further, check ratio and time available
  //if ratio greater than 90 and time available then reorder
  // 
  if (ddPackage == SYLVAN_PACKAGE){
    if (!noReordSinceGC){
      return false;
    }
    size_t preUsed = llmsset_count_marked(sylvan::nodes);
    // if the table is filled with more than x% data, start reordering
    if (preUsed > sylvan::nodes->max_size * reordThresh) { // comparing to max_size ensures table has grown fully
      printLine("Used fraction: "+to_string(preUsed)+"/ "+to_string(sylvan::nodes->max_size)+"*"+to_string(reordThresh)+" = "+to_string(preUsed/(sylvan::nodes->max_size*reordThresh)));
      return true;
    } else{
      return false;
    }
  } else{ //CUDD
    if (!noReordSinceGC){
      return false;
    }
    if(mgr->ReadMemoryInUse()<0.75*mgr->ReadMaxMemory()){
      return false;
    }
    if (mgr->ReadSlots() < lut){
      return false; // dont do reordering until fast growth is complete
    } else{
      Float usedFrac = Cudd_ReadUsedSlots(mgr->getManager());
      //printLine("Used fraction of slots: "+to_string(usedFrac)+"  ReordThresh:"+to_string(reordThresh)+" Expected usage: "+to_string(Cudd_ExpectedUsedSlots(mgr->getManager())));
      if (usedFrac > reordThresh){
        printLine("Used fraction of slots: "+to_string(usedFrac)+"  ReordThresh:"+to_string(reordThresh)+" Expected usage: "+to_string(Cudd_ExpectedUsedSlots(mgr->getManager())));
        return true;
      } else{
        return false;
      }
    }
  /*
    set loose upto to enable fast growth of unique table say to 80% of maxmem
      fraction can be 2 or 1.5
      it is given as number of slots
    then can compare cudd_readslots to check size of unique table
      this can give an idea of whether fast growth has terminated
    then use cudd_readusedslots to see what fraction of slots are used.
    or directly use Cudd_ReadMemoryInUse
    how to set timelimit, or maxvartoswaplimit or maxswapslimit?
      need to take into account time left, mem used 
      need a parameter like minImprovement, which says terminate reordering after at least minImp improvment
  */
  }
}

void Dd::afterReorder(){
  // if (ddPackage==SYLVAN_PACKAGE){
  //   Float ratio;
  //   sylvan::sylvan_table_usage_RUN(&postUsed, &postTotal);
  //   printLine("Post Used:"+to_string(postUsed)+"Post total:"+to_string(postTotal));
  //   ratio = (preUsed+0.0-postUsed)/preUsed;
  //   printLine("Ratio is:"+to_string(ratio) );
  // }
  if (didReordering){
    if (ddPackage == CUDD_PACKAGE){
      if (mgr->getManager()->ddTotalNumberSwapping >= maxSwaps){
        maxSwaps = mgr->getManager()->ddTotalNumberSwapping + maxSwapsInc;
      }
      reordThresh += reordThreshInc;
      reordThreshInc /= 2.5; //need to be slightly more aggressive than 2
    } else{ //Sylvan
      reordThresh += reordThreshInc;
      reordThreshInc /= 2;
      swapTime *= 2;
      maxSwaps += maxSwapsInc;
      #ifdef SYLDEV
      sylvan_set_reorder_timelimit(1 * swapTime * 1000);
      #endif
    }
    noReordSinceGC = false;
  }
  didReordering = false;
}

void Dd::manualReorderCUDD1(Map<Int, vector<Int>> levelMaps){
  DdManager* m = mgr->getManager();
  assert(ddReorderPreprocess(m));
  Int initialSize = m->keys - m->isolated;
  int* perm = new int[mgr->ReadSize()];
  Int minSize = initialSize;
  Int minVOInd = 0; //unchanged
  vector<Int> initialVO(mgr->ReadSize(),-1);
  for (int i = 0; i<mgr->ReadSize();i++){
    initialVO[mgr->ReadPerm(i)] = i; //inverted so as to allow easy restoring later
  }
  for (auto [i,vo] : levelMaps){
    assert (i!=0);
    assert(vo.size() == mgr->ReadSize());
    for (int j = 0; j<vo.size();j++){
      perm[j] = vo[j];
    }
    // cout<<"0\n";
    mgr->ShuffleHeap(perm);
    // cout<<"1\n";
    Int newSize = m->keys - m->isolated;
    if (newSize < minSize){
      minSize = newSize;
      minVOInd = i;
    }
  }
  if (minVOInd != 0){
    auto minVO = levelMaps[minVOInd];
    for (int j = 0; j<minVO.size();j++){
      perm[j] = minVO[j];
    }
    mgr->ShuffleHeap(perm);
  } else{   // restore initialVO
    assert(minSize == initialSize);
    for (int j = 0; j<initialVO.size();j++){
      perm[j] = initialVO[j];
    }
    mgr->ShuffleHeap(perm);
  }
  delete perm;
}

void Dd::manualReorderCUDD2(){
  mgr->ReduceHeap(CUDD_REORDER_SYMM_SIFT);
  // mgr->ReduceHeap(CUDD_REORDER_SIFT); //in one test this took 294 secs compared to 299 secs for the symmetric version. So not much difference.
}

void Dd::manualReorder(Map<Int, vector<Int>> levelMaps){
  if (dynVarOrdering != 1 && dynVarOrdering != 2){
    return;
  }
  if (beforeReorder()){
    didReordering = true;
    TimePoint reordStart = util::getTimePoint();
    printLine("Starting reordering..");
    if (ddPackage == SYLVAN_PACKAGE){
      #ifdef SYLDEV
      sylvan::Sylvan::reduceHeap();
      #endif
    } else{ //CUDD_PACKAGE
      printLine("NodeCount before reordering: "+to_string(mgr->ReadNodeCount()));
      if (dynVarOrdering == 1){
        manualReorderCUDD1(levelMaps);
      }else { //dynVarOrdering == 2
        manualReorderCUDD2();
      }
      printLine("NodeCount after reordering: "+to_string(mgr->ReadNodeCount())); 
    }
    afterReorder();
    printLine("Reordering done! Time taken: "+to_string(util::getDuration(reordStart)));
  }
}

void Dd::init(string ddPackage_, Int numVars, bool logCounting_, bool atomicAbstract_, bool weightedCounting_, bool multiplePrecision_, Int tableRatio, Int initRatio, Int threadCount, Float maxMem, Int dynVarOrdering_, Int dotFileIndex_){
  ddPackage = ddPackage_;
  logCounting = logCounting_;
  atomicAbstract = atomicAbstract_;
  weightedCounting = weightedCounting_;
  multiplePrecision = multiplePrecision_;
  dotFileIndex = dotFileIndex_;
  dynVarOrdering = dynVarOrdering_;
  
  reordThresh = 0.65;  reordThreshInc = 0.1;
  maxSwaps = 250; maxSwapsInc = 250; swapTime = 15;
  noReordSinceGC = false;
  
  if (ddPackage == CUDD_PACKAGE){
    mgr = new Cudd(
    numVars, // init num of BDD vars
    0, // init num of ZDD vars
    CUDD_UNIQUE_SLOTS, // init num of unique-table slots; cudd.h: #define CUDD_UNIQUE_SLOTS 256
    CUDD_CACHE_SLOTS, // init num of cache-table slots; cudd.h: #define CUDD_CACHE_SLOTS 262144
    maxMem * MEGA // maxMemory
    );
    mgr->SetMaxMemory(maxMem*MEGA);
    if (dynVarOrdering == 3){
      enableDynamicOrdering();
    }
    mgr->makeTerse(); // to prevent printing spurious lines whcih wont have "c o " prefix required by MCC
    Int lut = (maxMem*MEGA/sizeof(DdNode))/3;  //the default inside CUDD is 5. We want lut to be larger.
    mgr->SetLooseUpTo(lut);
    mgr->SetSiftMaxVar(10);
    mgr->SetSiftMaxSwap(maxSwaps);
    mgr->AddHook(postGCHook, CUDD_POST_GC_HOOK);
    mgr->AddHook(preGCHook, CUDD_PRE_GC_HOOK);
    printLine("CUDD Max Mem: "+to_string(mgr->ReadMaxMemory()));
    printLine("CUDD Max Cache Hard: "+to_string(mgr->ReadMaxCacheHard()));
  } else{
    lace_start(threadCount, 1000000); // auto-detect number of workers, use a 1,000,000 size task queue
    // Init Sylvan  
    sylvan_set_limits(maxMem*MEGA, tableRatio, initRatio);
    sylvan_init_package();
    sylvan_init_mtbdd();
    if (multiplePrecision) {
      sylvan::gmp_init();
    }
    sylvan::sylvan_gc_hook_pregc(TASK(gc_start));
    sylvan::sylvan_gc_hook_postgc(TASK(gc_end));
    #ifdef SYLDEV
    sylvan::mtbdd_newlevels(numVars);

    if(dynVarOrdering > 0){
      sylvan_init_reorder();

      sylvan_set_reorder_maxswap(maxSwaps);
      sylvan_set_reorder_maxvar(10);
      sylvan_set_reorder_threshold(256);
      sylvan_set_reorder_maxgrowth(1.2f);
      sylvan_set_reorder_timelimit(1 * swapTime * 1000);

      sylvan_re_hook_prere(TASK(reordering_start));
      sylvan_re_hook_progre(TASK(reordering_progress));
      sylvan_re_hook_postre(TASK(reordering_end));
      //gc hooks set above since they will be used even without reordering.
    }
    #endif
  }
}

void Dd::stop(){
  if (ddPackage == SYLVAN_PACKAGE) { // quits Sylvan
    sylvan::sylvan_quit();
    lace_stop();
  }else{
    printLine("Current CUDD used Mem: "+to_string(mgr->ReadMemoryInUse()));
    printLine("Total GC time: "+to_string(mgr->ReadGarbageCollectionTime()));
    //mgr->info();
    Cudd_Quit(mgr->getManager());
  }
}