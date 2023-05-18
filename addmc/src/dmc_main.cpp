#include "dmc.hpp"
#include "util.hpp"

using dpve::io::printRow;
using dpve::Number;
using dpve::io::InputParams;
using dpve::Dpve;
using dpve::io::printAdjustedSolutionRows;
using dpve::io::printAssignmentString;
using dpve::util::getTimePoint;
using dpve::util::getDuration;
using dpve::io::printAdjustedSolutionRows;
using dpve::TimePoint;

int main(int argc, char** argv) {
  std::cout << std::unitbuf; // enables automatic flushing
  dpve::io::printPreamble(argc,argv);
  InputParams p = dpve::io::parseOptions(argc,argv);
  assert(dpve::io::validateOptions(p));
  p.printParsed();
  try{
    Dpve d(p);
    auto [adjustedSolution, maximizer] = d.computeSolution();
    printAdjustedSolutionRows(adjustedSolution,p.pmParams.satSolverPruning,p.logCounting,p.weightedCounting,p.multiplePrecision,p.existRandom,p.projectedCounting);
    if (p.pmParams.maximizerFormat) {
      switch (p.pmParams.maximizerFormat) {
        case dpve::NEITHER_FORMAT:
          break;
        case dpve::SHORT_FORMAT:
          printAssignmentString(maximizer.getShortFormat(p.cnf.declaredVarCount));
          break;
        case dpve::LONG_FORMAT:
          printAssignmentString(maximizer.getLongFormat(p.cnf.declaredVarCount));
          break;
        default:
          printAssignmentString(maximizer.getShortFormat(p.cnf.declaredVarCount));
          printAssignmentString(maximizer.getLongFormat(p.cnf.declaredVarCount));
      }
      if (p.pmParams.maximizerVerification) {
        TimePoint maximizerVerificationStartPoint = getTimePoint();
        Number maximizerSolution = d.getMaximizerValue(maximizer);
        printRow("adjustedSolution", adjustedSolution);
        printRow("maximizerSolution", maximizerSolution);
        printRow("solutionMatch", (adjustedSolution - maximizerSolution).getAbsolute() < Number("1/1000000")); // 1e-6 is tolerance in ProCount paper
        if (p.verboseSolving >= 1) {
          printRow("maximizerVerificationSeconds", getDuration(maximizerVerificationStartPoint));
        }
      }
    }
  }
  catch (dpve::util::UnsatException) {
    dpve::io::printAdjustedSolutionRows(p.logCounting ? Number(-dpve::INF) : Number(),p.pmParams.satSolverPruning,p.logCounting,p.weightedCounting,p.multiplePrecision,p.existRandom,p.projectedCounting, true);
  }
  printRow("seconds", getDuration(p.toolStartPoint));
}