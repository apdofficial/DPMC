#pragma once

#include "../libraries/cryptominisat/build/include/cryptominisat5/cryptominisat.h"
#include "formula.hpp"

using CMSat::Lit;
using CMSat::lbool; // generally uint8_t; typically {l_True, l_False, l_Undef}
using CMSat::l_True;
using CMSat::l_False;

namespace dpve{

class SatSolver {
public:
  CMSat::SATSolver cmsSolver;

  bool checkSat(bool exceptionThrowing); // may throw UnsatSolverException
  Assignment getModel(); // also bans returned model for future solving
  Lit getLit(Int cnfVar, bool val);
  SatSolver(const Cnf& cnf);
};
} //end namespace dpve