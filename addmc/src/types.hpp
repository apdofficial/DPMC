#pragma once

#include <chrono>
#include <gmpxx.h>
#include <unordered_set>
#include <unordered_map>

/* uses ===================================================================== */

// using std::cout;
// using std::greater;
// using std::istream_iterator;
// using std::left;
// using std::map;
// using std::max;
// using std::min;
// using std::multimap;
// using std::next;
// using std::ostream;
// using std::pair;
using std::string;
// using std::to_string;
// using std::vector;

/* types ==================================================================== */
namespace dpve{
using Float = long double;
using Int = long long;
using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;

template<typename K, typename V> using Map = std::unordered_map<K, V>;
template<typename T> using Set = std::unordered_set<T>;


const Float INF = std::numeric_limits<Float>::infinity();

const Int MIN_INT = std::numeric_limits<Int>::min();
const Int MAX_INT = std::numeric_limits<Int>::max();

const Float MEGA = 1e6l; // same as countAntom (1 MB = 1e6 B)

class Number {
public:
  static bool multiplePrecision;
  
  mpq_class quotient;
  Float fraction;

  Number(const mpq_class& q); // multiplePrecision
  Number(Float f); // !multiplePrecision
  Number(const Number& n);
  Number(const string& repr = "0"); // `repr` is `<int>/<int>` or `<float>`

  Number getAbsolute() const;
  Float getLog10() const;
  Float getLogSumExp(const Number& n) const;
  bool operator==(const Number& n) const;
  bool operator!=(const Number& n) const;
  bool operator<(const Number& n) const;
  bool operator<=(const Number& n) const;
  bool operator>(const Number& n) const;
  bool operator>=(const Number& n) const;
  Number operator*(const Number& n) const;
  Number& operator*=(const Number& n);
  Number operator+(const Number& n) const;
  Number& operator+=(const Number& n);
  Number operator-(const Number& n) const;

  static Number mul_exp2(const Number n, const Int exp);
};
} //end namespace dpve