#pragma once

#include "common.hpp"
#include "formula.hpp"
#include <iostream>
/* consts =================================================================== */

namespace dpve{
  extern ostream& operator<<(ostream& stream, const Number& n);
namespace io{
  class PruneMaxParams{
    public:
      mutable Float logBound; //can change
      const Int maximizerFormat;
      const bool maximizerVerification;
      const bool satSolverPruning;
      const bool substitutionMaximization;
      const string thresholdModel;
      PruneMaxParams(const Float logBound, const Int maximizerFormat, const bool maximizerVerification, const bool satSolverPruning, 
        const bool substitutionMaximization, const string thresholdModel);
    private:
      PruneMaxParams();
  };
  class InputParams{
    public:
      const bool atomicAbstract;
      // const string cnfFilePath;
      const Cnf cnf;
      const string ddPackage;
      const Int ddVarOrderHeuristic;
      const Int dynVarOrdering;
      const bool existRandom;
      const Int initRatio; // log2(max_size / init_size)
      const string joinPriority;
      const bool logCounting;
      const bool multiplePrecision;
      const Float maxMem;
      const Float plannerWaitDuration;
      const bool projectedCounting;
      const PruneMaxParams pmParams;
      const Int randomSeed;
      const Int satFilter;
      const Int tableRatio; // log2(unique_table / cache_table)
      const Int threadCount;
      const TimePoint toolStartPoint;
      const Int verboseCnf;
      const Int verboseJoinTree;
      const Int verboseProfiling;
      const Int verboseSolving;
      const bool weightedCounting;
   
      void printParsed();
      InputParams(const bool atomicAbstract, const Cnf cnf, const string ddPackage, 
        const Int ddVarOrderHeuristic, const Int dynVarOrdering, const bool existRandom, const Int initRatio, const string joinPriority, 
        const bool logCounting, const bool multiplePrecision, const Float maxMem, const Float plannerWaitDuration, 
        const bool projectedCounting, const PruneMaxParams pmParams, const Int randomSeed, const Int satFilter, 
        const Int tableRatio, const Int threadCount, const TimePoint toolStartPoint, 
        const Int verboseCnf, const Int verboseJoinTree, const Int verboseProfiling, const Int verboseSolving, const bool weightedCounting);
    private:
      InputParams();
  };
  
  inline const string WARNING = "c MY_WARNING: ";
  inline const string DASH_LINE = "c ------------------------------------------------------------------\n";
  
  const InputParams parseOptions(int argc, char** argv);
  bool validateOptions(InputParams&);

  void printRowKey(const string& key, size_t keyWidth);

  template<typename T> void printRow(const string& key, const T& val, size_t keyWidth = 32){
    printRowKey(key, keyWidth);

    Int p = std::cout.precision();
    // cout.precision(std::numeric_limits<Float>::digits10); // default for Float: 6 digits
    std::cout << val << "\n";
    std::cout.precision(p);
  }

  void printPreamble(int argc, char** argv);
  void printInputLine(const string& line, Int lineIndex);
  
  void printTypeRow(size_t keyWidth, bool existRandom, bool weightedCounting, bool projectedCounting);
  void printEstRow(const Number& solution, size_t keyWidth, bool logCounting);
  void printSatRow(const Number& solution, bool unsatFlag, size_t keyWidth, bool satSolverPruning, bool logCounting, 
      bool weightedCounting, bool multiplePrecision);
  void printArbRow(const Number& solution, bool frac, size_t keyWidth, bool weightedCounting);
  void printDoubleRow(const Number& solution, size_t keyWidth, bool logCounting);
  void printAdjustedSolutionRows(const Number& adjustedSolution, bool satSolverPruning, bool logCounting, 
      bool weightedCounting, bool multiplePrecision, bool existRandom, bool projectedCounting, bool unsatFlag=false, size_t keyWidth=0);
  void printAssignmentString(string s);
  void printLine(string s);
  void printLine(char* s);
  void printLine();
} //namespace io
} //namespace dpve