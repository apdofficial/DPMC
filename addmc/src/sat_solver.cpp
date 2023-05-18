#include "util.hpp"
#include "sat_solver.hpp"

namespace dpve{
bool SatSolver::checkSat(bool exceptionThrowing) {
  lbool satisfiability = cmsSolver.solve();
  if (satisfiability == l_False) {
    if (exceptionThrowing) {
      throw dpve::util::UnsatSolverException();
    }
    return false;
  }
  assert(satisfiability == l_True);
  return true;
}

Assignment SatSolver::getModel() {
  vector<Lit> banLits;
  Assignment model;
  vector<lbool> lbools = cmsSolver.get_model();
  for (Int i = 0; i < lbools.size(); i++) {
    Int cnfVar = i + 1;
    bool val = true;
    lbool b = lbools.at(i);
    if (b != CMSat::l_Undef) {
      assert(b == l_True || b == l_False);
      val = b == l_True;
      banLits.push_back(getLit(cnfVar, !val));
    }
    cmsSolver.add_clause(banLits);
    model[cnfVar] = val;
  }
  return model;
}

Lit SatSolver::getLit(Int cnfVar, bool val) {
  assert(cnfVar > 0);
  return Lit(cnfVar - 1, !val);
}

SatSolver::SatSolver(const Cnf& cnf) {
  cmsSolver.new_vars(cnf.declaredVarCount);
  for (const Clause& clause : cnf.clauses) {
    if (clause.xorFlag) {
      vector<unsigned> vars;
      bool rhs = true;
      for (Int literal : clause) {
        vars.push_back(abs(literal) - 1);
        if (literal < 0) {
          rhs = !rhs;
        }
      }
      cmsSolver.add_xor_clause(vars, rhs);
    }
    else {
      vector<Lit> lits;
      for (Int literal : clause) {
        lits.push_back(getLit(abs(literal), literal > 0));
      }
      cmsSolver.add_clause(lits);
    }
  }
}
}