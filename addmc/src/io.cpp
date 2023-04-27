#include "io.hpp"


#include "../../addmc/libraries/cxxopts/include/cxxopts.hpp"
#include "formula.hpp"
#include "util.hpp"

#include <cmath>
#include <iomanip>
#include <thread>

using dpve::io::printRow;
using dpve::io::InputParams;
using dpve::io::PruneMaxParams;
using dpve::Int;
using std::cout;
using std::left;
using std::max;
using std::min;
using std::right;
using std::setw;
using std::to_string;

using std::thread;

namespace {  // anonymous namespace. Local to this file

  const string ATOMIC_ABSTRACT_FLAG = "aa";
  const string CNF_FILE_FLAG = "cf";
  const string DD_PACKAGE_FLAG = "dp";
  const string DD_VAR_FLAG = "dv";
  const string DYN_ORDER_FLAG = "dy";
  const string EXIST_RANDOM_FLAG = "er";
  const string INIT_RATIO_FLAG = "ir";
  const string HELP_FLAG = "h";
  const string JOIN_PRIORITY_FLAG = "jp";
  const string LOG_BOUND_FLAG = "lb";
  const string LOG_COUNTING_FLAG = "lc";
  const string MAXIMIZER_FORMAT_FLAG = "mf";
  const string MAX_MEM_FLAG = "mm";
  const string MULTIPLE_PRECISION_FLAG = "mp";
  const string MEM_SENSITIVITY_FLAG = "ms";
  const string MAXIMIZER_VERIFICATION_FLAG = "mv";
  const string PROJECTED_COUNTING_FLAG = "pc";
  const string PLANNER_WAIT_FLAG = "pw";
  const string RANDOM_SEED_FLAG = "rs";
  const string SAT_FILTER_FLAG = "sa";
  const string SUBSTITUTION_MAXIMIZATION_FLAG = "sm";
  const string SAT_SOLVER_PRUNING = "sp";
  const string SLICE_VAR_FLAG = "sv";
  const string THREAD_COUNT_FLAG = "tc";
  const string THRESHOLD_MODEL_FLAG = "tm";
  const string TABLE_RATIO_FLAG = "tr";
  const string THREAD_SLICE_COUNT_FLAG = "ts";
  const string VERBOSE_CNF_FLAG = "vc";
  const string VERBOSE_JOIN_TREE_FLAG = "vj";
  const string VERBOSE_PROFILING_FLAG = "vp";
  const string VERBOSE_SOLVING_FLAG = "vs";
  const string WEIGHTED_COUNTING_FLAG = "wc";
  
  class OptionRequirement {
    public:
      string name;
      string value;
      string comparator;

      OptionRequirement(const string& name, const string& value, const string& comparator = "=");
      string getRequirement() const;
  };

  OptionRequirement::OptionRequirement(const string& name, const string& value, const string& comparator) {
    this->name = name;
    this->value = value;
    this->comparator = comparator;
  }

  string OptionRequirement::getRequirement() const {
    return name + "_arg " + comparator + " " + value;
  }

  string requireOptions(const vector<OptionRequirement>& requirements) {
    string s = " [needs ";
    for (auto it = requirements.begin(); it != requirements.end(); it++) {
      s += it->getRequirement();
      if (next(it) == requirements.end()) {
        s += "]";
      }
      else {
        s += ", ";
      }
    }
    return s;
  }

  string requireOption(const string& name, const string& value, const string& comparator = "=") {
    return requireOptions({OptionRequirement(name, value, comparator)});
  }

  string requireDdPackage(const string& ddPackageArg) {
    assert(dpve::DD_PACKAGES.contains(ddPackageArg));
    return requireOption(DD_PACKAGE_FLAG, ddPackageArg);
  }

  string helpDdPackage() {
    string s = "diagram package: ";
    for (auto it = dpve::DD_PACKAGES.begin(); it != dpve::DD_PACKAGES.end(); it++) {
      s += it->first + "/" + it->second;
      if (next(it) != dpve::DD_PACKAGES.end()) {
        s += ", ";
      }
    }
    return s + "; string";
  }

  map<Int, string> getVarOrderHeuristics() {
    map<Int, string> m = dpve::CNF_VAR_ORDER_HEURISTICS;
    m.insert(dpve::JOIN_TREE_VAR_ORDER_HEURISTICS.begin(), dpve::JOIN_TREE_VAR_ORDER_HEURISTICS.end());
    return m;
  }

  string helpLogBound() {
    string s = "log10(bound) for pruning";
    s += requireOptions({
      OptionRequirement(PROJECTED_COUNTING_FLAG, "0"),
      OptionRequirement(EXIST_RANDOM_FLAG, "1"),
      OptionRequirement(LOG_COUNTING_FLAG, "1")
    });
    return s + "; float";
  }

  string helpThresholdModel() {
    string s = "threshold model for pruning";
    s += requireOptions({
      OptionRequirement(PROJECTED_COUNTING_FLAG, "0"),
      OptionRequirement(EXIST_RANDOM_FLAG, "1"),
      OptionRequirement(LOG_COUNTING_FLAG, "1"),
      OptionRequirement(LOG_BOUND_FLAG, "-inf")
    });
    return s + "; string";
  }

  string helpSatSolverPruning() {
    string s = "SAT pruning with CryptoMiniSat";
    s += requireOptions({
      OptionRequirement(PROJECTED_COUNTING_FLAG, "0"),
      OptionRequirement(EXIST_RANDOM_FLAG, "1"),
      OptionRequirement(LOG_COUNTING_FLAG, "1"),
      OptionRequirement(LOG_BOUND_FLAG, "-inf"),
      OptionRequirement(THRESHOLD_MODEL_FLAG, "\"\""),
    });
    return s + ": 0, 1; int";
  }

  string helpMaximizerFormat() {
    string s = "maximizer format";
    s += requireOptions({
      OptionRequirement(EXIST_RANDOM_FLAG, "1"),
      OptionRequirement(DD_PACKAGE_FLAG, dpve::CUDD_PACKAGE)
    });
    s += ": ";
    for (auto it = dpve::MAXIMIZER_FORMATS.begin(); it != dpve::MAXIMIZER_FORMATS.end(); it++) {
      s += to_string(it->first) + "/" + it->second;
      if (next(it) != dpve::MAXIMIZER_FORMATS.end()) {
        s += ", ";
      }
    }
    return s + "; int";
  }

  string helpSubstitutionMaximization() {
    string s = "substitution-based maximization";
    s += requireOptions({
      OptionRequirement(WEIGHTED_COUNTING_FLAG, "0"),
      OptionRequirement(MAXIMIZER_FORMAT_FLAG, to_string(dpve::NEITHER_FORMAT), ">")
    });
    return s + ": 0, 1; int";
  }

  string helpVarOrderHeuristic(const map<Int, string>& heuristics) {
    string s = ": ";
    for (auto it = heuristics.begin(); it != heuristics.end(); it++) {
      s += to_string(it->first) + "/" + it->second;
      if (next(it) != heuristics.end()) {
        s += ", ";
      }
    }
    return s + " (negatives for inverse orders); int";
  }

  string helpVerboseCnfProcessing() {
    return "verbose CNF processing: 0, 1, 2, 3; int";
  }

  string helpVerboseSolving() {
    return "verbose solving: 0, 1, 2; int";
  }

  string helpDiagramVarOrderHeuristic() {
    return "diagram var order" + helpVarOrderHeuristic(dpve::CNF_VAR_ORDER_HEURISTICS);
  }

  string helpSliceVarOrderHeuristic() {
    string s = "slice var order";
    s += requireOption(THREAD_SLICE_COUNT_FLAG, "1", ">");
    s += helpVarOrderHeuristic(getVarOrderHeuristics());
    return s;
  }

  string helpJoinPriority() {
    string s = "join priority: ";
    for (auto it = dpve::JOIN_PRIORITIES.begin(); it != dpve::JOIN_PRIORITIES.end(); it++) {
      s += it->first + "/" + it->second;
      if (next(it) != dpve::JOIN_PRIORITIES.end()) {
        s += ", ";
      }
    }
    return s + "; string";
  }

  string helpDynamicVarOrdering() {
    return "dynamic variable ordering. DD_PACKAGE must be CUDD. 0/1. Default 0.";
  }

  string helpSatFilter() {
    return "0 - Disable SatFilter (Only Executor) / 1 - Only SatFilter / 2 - SatFilter + Executor. Default 0.";
  }

  string helpAtomicAbstract() {
    return "0/1 - Disable/Enable single step abstraction operation. (to enable, er_arg=0 pc_arg=0 dp_arg=c or if dp_arg=s then additionally wc_arg=0 required) Default 0.";
  }
}

void dpve::io::printPreamble(int argc, char** argv){
   std::cout << "name of program: " << argv[0] << '\n' ;
   std::cout << "there are " << argc-1 << " (more) arguments, they are:\n" ;
   std::copy( argv+1, argv+argc, std::ostream_iterator<const char*>( std::cout, " " ) ) ;
   cout<<"\n";
}

void dpve::io::printInputLine(const string& line, Int lineIndex) {
  cout << "c line " << right << setw(5) << lineIndex << ":" << (line.empty() ? "" : " " + line) << "\n";
}

void dpve::io::printRowKey(const string& key, size_t keyWidth) {
  string prefix = key;
  if (key != "s") {
    prefix = "c " + prefix;
  }

  keyWidth = max(keyWidth, prefix.size() + 1);
  cout << left << setw(keyWidth) << prefix;
}

void dpve::io::printTypeRow(size_t keyWidth, bool existRandom, bool weightedCounting, bool projectedCounting) {
  string type = "maximum";
  if (!existRandom) {
    type = "mc";
    if (weightedCounting) {
      type = "w" + type;
    }
    if (projectedCounting) {
      type = "p" + type;
    }
  }
  printRow("s type", type, keyWidth);
}

void dpve::io::printEstRow(const Number& solution, size_t keyWidth, bool logCounting) {
  printRow("s log10-estimate", logCounting ? solution.fraction : solution.getLog10(), keyWidth);
}

void dpve::io::printSatRow(const Number& solution, bool unsatFlag, size_t keyWidth, bool satSolverPruning, bool logCounting, 
  bool weightedCounting, bool multiplePrecision) {
  
  const string SAT_WORD = "SATISFIABLE";
  const string UNSAT_WORD = "UN" + SAT_WORD;

  string satisfiability = "UNKNOWN";

  if (unsatFlag) {
    satisfiability = UNSAT_WORD;
  }
  else if (satSolverPruning) {
    satisfiability = SAT_WORD; // otherwise, UnsatSolverException would have been thrown earlier
  }
  else if (logCounting) {
    if (solution == Number(-INF)) {
      if (!weightedCounting) {
        satisfiability = UNSAT_WORD;
      }
    }
    else {
      satisfiability = SAT_WORD;
    }
  }
  else {
    if (solution == Number()) {
      if (!weightedCounting || multiplePrecision) {
        satisfiability = UNSAT_WORD;
      }
    }
    else {
      satisfiability = SAT_WORD;
    }
  }

  printRow("s", satisfiability, keyWidth);
}

void dpve::io::printArbRow(const Number& solution, bool frac, size_t keyWidth, bool weightedCounting) {
  string key = "s exact arb ";

  if (weightedCounting) {
    if (frac) {
      printRow(key + "frac", solution, keyWidth);
    }
    else {
      printRow(key + "float", mpf_class(solution.quotient), keyWidth);
    }
  }
  else {
    printRow(key + "int", solution, keyWidth);
  }
}

void dpve::io::printDoubleRow(const Number& solution, size_t keyWidth, bool logCounting) {
  Float f = solution.fraction;
  printRow("s exact double prec-sci", logCounting ? exp10l(f) : f, keyWidth);
}

void dpve::io::printAdjustedSolutionRows(const Number& adjustedSolution, bool satSolverPruning, bool logCounting, 
  bool weightedCounting, bool multiplePrecision, bool existRandom, bool projectedCounting, bool unsatFlag, size_t keyWidth) {
  cout << DASH_LINE;
  
  printSatRow(adjustedSolution, unsatFlag, keyWidth, satSolverPruning, logCounting, weightedCounting, multiplePrecision);
  printTypeRow(keyWidth,existRandom, weightedCounting, projectedCounting);
  printEstRow(adjustedSolution, keyWidth, logCounting);

  if (multiplePrecision) {
    printArbRow(adjustedSolution, false, keyWidth, weightedCounting); // notation = weighted ? int : float
    if (weightedCounting) {
      printArbRow(adjustedSolution, true, keyWidth, weightedCounting); // notation = frac
    }
  }
  else {
    printDoubleRow(adjustedSolution, keyWidth, logCounting);
  }
  cout << DASH_LINE;
}

void dpve::io::printAssignmentString(string s){
  cout << "v";
  cout << s;
  cout << "\n";
}

void dpve::io::printLine(){
  cout << "\n";
}

void dpve::io::printLine(string s){
  cout << "c c "<<s<<std::endl<<std::flush;
}

void dpve::io::printLine(char *s){
  cout << "c c "<<s<<std::endl<<std::flush;
}

PruneMaxParams::PruneMaxParams(const Float logBound, const Int maximizerFormat, const bool maximizerVerification, 
      const bool satSolverPruning, const bool substitutionMaximization, const string thresholdModel):
      logBound(logBound), maximizerFormat(maximizerFormat), maximizerVerification(maximizerVerification), satSolverPruning(satSolverPruning),
      substitutionMaximization(substitutionMaximization), thresholdModel(thresholdModel)
        {}

InputParams::InputParams(const bool atomicAbstract, const Cnf cnf, const string ddPackage, 
    const Int ddVarOrderHeuristic, const Int dynVarOrdering, const bool existRandom, const Int initRatio, const string joinPriority, const bool logCounting,
    const bool multiplePrecision, const Float maxMem, const Float plannerWaitDuration, 
    const bool projectedCounting, const PruneMaxParams pmParams, const Int randomSeed, const Int satFilter, 
    const Int tableRatio, const Int threadCount, const TimePoint toolStartPoint, 
    const Int verboseCnf, const Int verboseJoinTree, const Int verboseProfiling, const Int verboseSolving, const bool weightedCounting):
    
    atomicAbstract(atomicAbstract),
    cnf(cnf),
    ddPackage(ddPackage),
    ddVarOrderHeuristic(ddVarOrderHeuristic),
    dynVarOrdering(dynVarOrdering),
    initRatio(initRatio),
    joinPriority(joinPriority),
    multiplePrecision(multiplePrecision),
    maxMem(maxMem),
    plannerWaitDuration(plannerWaitDuration),
    existRandom(existRandom),
    logCounting(logCounting),
    projectedCounting(projectedCounting),
    pmParams(pmParams),
    randomSeed(randomSeed),
    satFilter(satFilter),
    tableRatio(tableRatio),
    threadCount(threadCount),
    toolStartPoint(toolStartPoint),
    verboseCnf(verboseCnf),
    verboseJoinTree(verboseJoinTree),
    verboseProfiling(verboseProfiling),
    verboseSolving(verboseSolving),
    weightedCounting(weightedCounting) 
    {}

const InputParams dpve::io::parseOptions(int argc, char** argv) {
  // print all command line arguments
   
  cxxopts::Options options("dmc", "Diagram Model Counter (reads join tree from stdin)");
  options.set_width(118);

  using cxxopts::value;
  options.add_options()
    (CNF_FILE_FLAG, "CNF file path; string (required)", value<string>())
    (WEIGHTED_COUNTING_FLAG, "weighted counting: 0, 1; int", value<Int>()->default_value("1"))
    (PROJECTED_COUNTING_FLAG, "projected counting (graded join tree): 0, 1; int", value<Int>()->default_value("0"))
    (EXIST_RANDOM_FLAG, "exist-random SAT (max-sum instead of sum-max): 0, 1; int", value<Int>()->default_value("0"))
    (DD_PACKAGE_FLAG, helpDdPackage(), value<string>()->default_value(CUDD_PACKAGE))
    (LOG_COUNTING_FLAG, "logarithmic counting: 0, 1; int", value<Int>()->default_value("0"))
    (LOG_BOUND_FLAG, helpLogBound(), value<string>()->default_value(to_string(-INF))) // cxxopts fails to parse "-inf" as Float
    (THRESHOLD_MODEL_FLAG, helpThresholdModel(), value<string>()->default_value(""))
    (SAT_SOLVER_PRUNING, helpSatSolverPruning(), value<Int>()->default_value("0"))
    (MAXIMIZER_FORMAT_FLAG, helpMaximizerFormat(), value<Int>()->default_value(to_string(NEITHER_FORMAT)))
    (MAXIMIZER_VERIFICATION_FLAG, "maximizer verification" + requireOption(MAXIMIZER_FORMAT_FLAG, to_string(NEITHER_FORMAT), ">") + ": 0, 1; int", value<Int>()->default_value("0"))
    (SUBSTITUTION_MAXIMIZATION_FLAG, helpSubstitutionMaximization(), value<Int>()->default_value("0"))
    (PLANNER_WAIT_FLAG, "planner wait duration minimum (in seconds); float", value<Float>()->default_value("0.0"))
    (THREAD_COUNT_FLAG, "thread count [or 0 for hardware_concurrency value]; int", value<Int>()->default_value("1"))
    (RANDOM_SEED_FLAG, "random seed; int", value<Int>()->default_value("0"))
    (DYN_ORDER_FLAG, helpDynamicVarOrdering(), value<Int>()->default_value("0"))
    (SAT_FILTER_FLAG, helpSatFilter(), value<Int>()->default_value("0"))
    (ATOMIC_ABSTRACT_FLAG, helpAtomicAbstract(), value<Int>()->default_value("0"))
    (DD_VAR_FLAG, helpDiagramVarOrderHeuristic(), value<Int>()->default_value(to_string(MCS_HEURISTIC)))
    (MAX_MEM_FLAG, "maximum memory (in MB) for unique table and cache table combined [or 0 for unlimited memory with CUDD]; float", value<Float>()->default_value("4e3"))
    (TABLE_RATIO_FLAG, "table ratio" + requireDdPackage(SYLVAN_PACKAGE) + ": log2(unique_size/cache_size); int", value<Int>()->default_value("1"))
    (INIT_RATIO_FLAG, "init ratio for tables" + requireDdPackage(SYLVAN_PACKAGE) + ": log2(max_size/init_size); int", value<Int>()->default_value("10"))
    (MULTIPLE_PRECISION_FLAG, "multiple precision" + requireDdPackage(SYLVAN_PACKAGE) + ": 0, 1; int", value<Int>()->default_value("0"))
    (JOIN_PRIORITY_FLAG, helpJoinPriority(), value<string>()->default_value(SMALLEST_PAIR))
    (VERBOSE_CNF_FLAG, helpVerboseCnfProcessing(), value<Int>()->default_value("0"))
    (VERBOSE_JOIN_TREE_FLAG, "verbose join-tree processing: 0, 1, 2", value<Int>()->default_value("0"))
    (VERBOSE_PROFILING_FLAG, "verbose profiling: 0, 1, 2; int", value<Int>()->default_value("0"))
    (VERBOSE_SOLVING_FLAG, helpVerboseSolving(), value<Int>()->default_value("0"))
    (HELP_FLAG, "help")
  ;

  cxxopts::ParseResult result = options.parse(argc, argv);
  if (result.count(HELP_FLAG) || !result.count(CNF_FILE_FLAG)) {
    cout << options.help();
    exit(0);
  }
  auto cnfFilePath = result[CNF_FILE_FLAG].as<string>();
  auto weightedCounting = result[WEIGHTED_COUNTING_FLAG].as<Int>(); 
  auto projectedCounting = result[PROJECTED_COUNTING_FLAG].as<Int>(); 
  auto existRandom = result[EXIST_RANDOM_FLAG].as<Int>(); 
  auto ddPackage = result[DD_PACKAGE_FLAG].as<string>(); 
  auto logCounting = result[LOG_COUNTING_FLAG].as<Int>();
  auto logBound = stold(result[LOG_BOUND_FLAG].as<string>());
  auto thresholdModel = result[THRESHOLD_MODEL_FLAG].as<string>(); // global var
  auto satSolverPruning = result[SAT_SOLVER_PRUNING].as<Int>(); // global var
  auto maximizerFormat = result[MAXIMIZER_FORMAT_FLAG].as<Int>(); // global var
  auto maximizerVerification = result[MAXIMIZER_VERIFICATION_FLAG].as<Int>(); // global var
  auto substitutionMaximization = result[SUBSTITUTION_MAXIMIZATION_FLAG].as<Int>(); // global var
  auto plannerWaitDuration = result[PLANNER_WAIT_FLAG].as<Float>();
    plannerWaitDuration = max(plannerWaitDuration, 0.0l);
  auto threadCount = result[THREAD_COUNT_FLAG].as<Int>(); // global var
  if (threadCount <= 0) {
    threadCount = thread::hardware_concurrency();
  }
  auto randomSeed = result[RANDOM_SEED_FLAG].as<Int>(); // global var
  auto dynVarOrdering = result[DYN_ORDER_FLAG].as<Int>();
  auto satFilter = result[SAT_FILTER_FLAG].as<Int>();
  auto atomicAbstract = result[ATOMIC_ABSTRACT_FLAG].as<Int>();
  auto ddVarOrderHeuristic = result[DD_VAR_FLAG].as<Int>();
  auto maxMem = result[MAX_MEM_FLAG].as<Float>(); // global var
    maxMem = max(maxMem, 0.0l);
  auto tableRatio = result[TABLE_RATIO_FLAG].as<Int>();
  auto initRatio = result[INIT_RATIO_FLAG].as<Int>();
  auto multiplePrecision = result[MULTIPLE_PRECISION_FLAG].as<Int>(); // global var
  assert(!result.count(TABLE_RATIO_FLAG) || ddPackage == SYLVAN_PACKAGE);
  assert(!result.count(INIT_RATIO_FLAG) || ddPackage == SYLVAN_PACKAGE);
  auto joinPriority = result[JOIN_PRIORITY_FLAG].as<string>(); //global var
  auto verboseCnf = result[VERBOSE_CNF_FLAG].as<Int>(); // global var
  auto verboseJoinTree = result[VERBOSE_JOIN_TREE_FLAG].as<Int>(); // global var
  auto verboseProfiling = result[VERBOSE_PROFILING_FLAG].as<Int>(); // global var
  auto verboseSolving = result[VERBOSE_SOLVING_FLAG].as<Int>(); // global var
  auto toolStartPoint = util::getTimePoint(); // global var

  PruneMaxParams pmParams(logBound,maximizerFormat,maximizerVerification,satSolverPruning,substitutionMaximization,thresholdModel);

  Cnf cnf(verboseCnf,randomSeed,weightedCounting,projectedCounting);
  cnf.readCnfFile(cnfFilePath);
  return InputParams(atomicAbstract, cnf, ddPackage, ddVarOrderHeuristic, dynVarOrdering, existRandom, initRatio, joinPriority, logCounting, multiplePrecision, maxMem, plannerWaitDuration, projectedCounting, pmParams, randomSeed, satFilter, tableRatio, threadCount, toolStartPoint, verboseCnf, verboseJoinTree, verboseProfiling, verboseSolving, weightedCounting);
}

bool dpve::io::validateOptions(InputParams& p){
  assert(DD_PACKAGES.contains(p.ddPackage));
  assert(p.pmParams.logBound == -INF || !p.projectedCounting);
  assert(p.pmParams.logBound == -INF || p.existRandom);
  assert(p.pmParams.logBound == -INF || p.logCounting);
  assert(p.pmParams.thresholdModel.empty() || !p.projectedCounting);
  assert(p.pmParams.thresholdModel.empty() || p.existRandom);
  assert(p.pmParams.thresholdModel.empty() || p.logCounting);
  assert(p.pmParams.thresholdModel.empty() || p.pmParams.logBound == -INF);
  assert(!p.pmParams.satSolverPruning || !p.projectedCounting);
  assert(!p.pmParams.satSolverPruning || p.existRandom);
  assert(!p.pmParams.satSolverPruning || p.logCounting);
  assert(!p.pmParams.satSolverPruning || p.pmParams.logBound == -INF);
  assert(!p.pmParams.satSolverPruning || p.pmParams.thresholdModel.empty());
  assert(MAXIMIZER_FORMATS.contains(p.pmParams.maximizerFormat));
  assert(!p.pmParams.maximizerFormat || p.existRandom);
  assert(!p.pmParams.maximizerFormat || p.ddPackage == CUDD_PACKAGE);
  assert(!p.pmParams.maximizerVerification || p.pmParams.maximizerFormat);
  assert(!p.pmParams.substitutionMaximization || !p.weightedCounting);
  assert(!p.pmParams.substitutionMaximization || p.pmParams.maximizerFormat);
  assert(p.threadCount > 0);
  assert(p.dynVarOrdering == 0 || p.ddPackage == CUDD_PACKAGE);
  assert(p.satFilter >= 0 && p.satFilter <=2);
  assert((p.atomicAbstract == false) || (p.projectedCounting == false && p.existRandom == false && p.ddPackage == CUDD_PACKAGE) || 
        (p.projectedCounting == false && p.existRandom == false && p.weightedCounting == false && p.ddPackage == SYLVAN_PACKAGE));
  //assert(CNF_VAR_ORDER_HEURISTICS.contains(abs(ddVarOrderHeuristic)));
  // assert(!result.count(SLICE_VAR_FLAG) || p.threadSliceCount > 1);
  // assert(util:SLICE_VAR:getVarOrderHeuristics().contains(abs(p.sliceVarOrderHeuristic)));
  assert(!p.multiplePrecision || p.ddPackage == SYLVAN_PACKAGE);
  assert(JOIN_PRIORITIES.contains(p.joinPriority));
  assert(p.verboseProfiling <= 0 || p.threadCount == 1);
  return true;
}

void dpve::io::InputParams::printParsed(){
  if (verboseSolving >= 1) {
    cout << "c processing command-line options...\n";
    // printRow("cnfFile", cnfFilePath);
    printRow("weightedCounting", weightedCounting);
    printRow("projectedCounting", projectedCounting);
    printRow("existRandom", existRandom);
    printRow("diagramPackage", DD_PACKAGES.at(ddPackage));
    if (ddPackage == CUDD_PACKAGE) {
      printRow("logCounting", logCounting);
      printRow("dynamic var ordering", dynVarOrdering);
      printRow("atomic abstract",atomicAbstract);
    }
    if (!projectedCounting && existRandom && logCounting) {
      if (pmParams.logBound > -INF) {
        printRow("logBound", pmParams.logBound);
      }
      else if (!pmParams.thresholdModel.empty()) {
        printRow("thresholdModel", pmParams.thresholdModel);
      }
      else if (pmParams.satSolverPruning) {
        printRow("satSolverPruning", pmParams.satSolverPruning);
      }
    }
    if (existRandom && ddPackage == CUDD_PACKAGE) {
      printRow("maximizerFormat", MAXIMIZER_FORMATS.at(pmParams.maximizerFormat));
    }
    if (pmParams.maximizerFormat) {
      printRow("maximizerVerification", pmParams.maximizerVerification);
    }
    if (!weightedCounting && pmParams.maximizerFormat) {
      printRow("substitutionMaximization", pmParams.substitutionMaximization);
    }
    printRow("plannerWaitSeconds", plannerWaitDuration);
    printRow("threadCount", threadCount);
    printRow("randomSeed", randomSeed);
    printRow("diagramVarOrderHeuristic", (ddVarOrderHeuristic < 0 ? "INVERSE_" : "TODO!!"));// + CNF_VAR_ORDER_HEURISTICS.at(abs(ddVarOrderHeuristic)));
    printRow("maxMemMegabytes", maxMem);
    if (ddPackage == SYLVAN_PACKAGE) {
      printRow("tableRatio", tableRatio);
      printRow("initRatio", initRatio);
      printRow("multiplePrecision", multiplePrecision);
    }
    printRow("joinPriority", JOIN_PRIORITIES.at(joinPriority));
    cout << "\n";
  }
}