#include "util.hpp"
#include "io.hpp"
#include <iterator>
#include <sstream>

using dpve::Int;
using dpve::TimePoint;
using dpve::Float;

using std::istream_iterator;

vector<Int> dpve::util::getSortedNums(const Set<Int>& nums) {
  vector<Int> v(nums.begin(), nums.end());
  sort(v.begin(), v.end());
  return v;
}

TimePoint dpve::util::getTimePoint() {
  return std::chrono::steady_clock::now();
}

Float dpve::util::getDuration(TimePoint start) {
  return std::chrono::duration_cast<std::chrono::milliseconds>(getTimePoint() - start).count() / 1e3l;
}

vector<string> dpve::util::splitInputLine(const string& line) {
  std::istringstream inStringStream(line);
  vector<string> words;
  copy(istream_iterator<string>(inStringStream), istream_iterator<string>(), back_inserter(words));
  return words;
}

/* classes for exceptions =================================================== */

/* class EmptyClauseException =============================================== */

dpve::util::EmptyClauseException::EmptyClauseException(Int lineIndex, const string& line) {
  cout << dpve::io::WARNING << "empty clause | line " << lineIndex << ": " << line << "\n";
}

/* class EmptyClauseException =============================================== */

dpve::util::UnsatSolverException::UnsatSolverException() {
  cout << dpve::io::WARNING << "unsatisfiable CNF, according to SAT solver\n";
}
