#include <tuple>
#include "util.hpp"
#include "common.hpp"

using dpve::Int;
using std::ostream;

/* global functions ========================================================= */

ostream& dpve::operator<<(ostream& stream, const dpve::Number& n) {
  if (dpve::Number::multiplePrecision) {
    stream << n.quotient;
  }
  else {
    stream << n.fraction;
  }

  return stream;
}
