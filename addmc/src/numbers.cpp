#include "types.hpp"

#include <cassert>
#include <cmath>

using dpve::Number;
using dpve::Float;
using std::max;

bool Number::multiplePrecision = 0;

Number::Number(const mpq_class& q) {
  assert(multiplePrecision);
  quotient = q;
}

Number::Number(Float f) {
  assert(!multiplePrecision);
  fraction = f;
}

Number::Number(const Number& n) {
  if (multiplePrecision) {
    *this = Number(n.quotient);
  }
  else {
    *this = Number(n.fraction);
  }
}

Number::Number(const string& repr) {
  Int divPos = repr.find('/');
  if (multiplePrecision) {
    if (divPos != string::npos) { // repr is <int>/<int>
      *this = Number(mpq_class(repr));
    }
    else { // repr is <float>
      *this = Number(mpq_class(mpf_class(repr)));
    }
  }
  else {
    if (divPos != string::npos) { // repr is <int>/<int>
      Float numerator = stold(repr.substr(0, divPos));
      Float denominator = stold(repr.substr(divPos + 1));
      *this = Number(numerator / denominator);
    }
    else { // repr is <float>
      *this = Number(stold(repr));
    }
  }
}

Number Number::getAbsolute() const {
  if (multiplePrecision) {
    return Number(abs(quotient));
  }
  return Number(fabsl(fraction));
}

Float Number::getLog10() const {
  if (multiplePrecision) {
    mpf_t f; // C interface
    mpf_init(f);
    mpf_set_q(f, quotient.get_mpq_t());
    long int exponent;
    Float d = mpf_get_d_2exp(&exponent, f); // f == d * 2^exponent
    Float lgF = log10l(d) + exponent * log10l(2);
    mpf_clear(f);
    return lgF;
  }
  return log10l(fraction);
}

Float Number::getLogSumExp(const Number& n) const {
  assert(!multiplePrecision);
  if (fraction == -INF) {
    return n.fraction;
  }
  if (n.fraction == -INF) {
    return fraction;
  }
  Float m = max(fraction, n.fraction);
  return log10l(exp10l(fraction - m) + exp10l(n.fraction - m)) + m; // base-10 Cudd_addLogSumExp
}

Number Number::mul_exp2(const Number n, const Int exp){
  if (multiplePrecision){
    mpz_t pow;
    mpz_init(pow);
    mpz_ui_pow_ui(pow,2,exp);
    mpz_class factor(pow);
    mpz_clear(pow);
    mpq_class res(n.quotient);
    res *= factor;
    res /= n.quotient.get_den();
    res.canonicalize();
    return Number(res);
  }
  return Number(n*(pow(2,exp)));
}

bool Number::operator==(const Number& n) const {
  if (multiplePrecision) {
    return quotient == n.quotient;
  }
  return fraction == n.fraction;
}

bool Number::operator!=(const Number& n) const {
  return !(*this == n);
}

bool Number::operator<(const Number& n) const {
  if (multiplePrecision) {
    return quotient < n.quotient;
  }
  return fraction < n.fraction;
}

bool Number::operator<=(const Number& n) const {
  return *this < n || *this == n;
}

bool Number::operator>(const Number& n) const {
  if (multiplePrecision) {
    return quotient > n.quotient;
  }
  return fraction > n.fraction;
}

bool Number::operator>=(const Number& n) const {
  return *this > n || *this == n;
}

Number Number::operator*(const Number& n) const {
  if (multiplePrecision) {
    return Number(quotient * n.quotient);
  }
  return Number(fraction * n.fraction);
}

Number& Number::operator*=(const Number& n) {
  *this = *this * n;
  return *this;
}

Number Number::operator+(const Number& n) const {
  if (multiplePrecision) {
    return Number(quotient + n.quotient);
  }
  return Number(fraction + n.fraction);
}

Number& Number::operator+=(const Number& n) {
  *this = *this + n;
  return *this;
}

Number Number::operator-(const Number& n) const {
  if (multiplePrecision) {
    return Number(quotient - n.quotient);
  }
  return Number(fraction - n.fraction);
}
