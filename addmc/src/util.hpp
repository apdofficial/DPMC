#pragma once

#include "types.hpp"

#include <iostream>
#include <map>
#include <vector>
#include <string>

/* namespaces =============================================================== */
using std::cout;
using std::vector;
using std::multimap;
using std::pair;
using std::greater;
namespace dpve::util {
  class UnsatException : public std::exception {};

  class EmptyClauseException : public UnsatException {
  public:
    EmptyClauseException(Int lineIndex, const string& line);
  };

  class UnsatSolverException : public UnsatException {
  public:
    UnsatSolverException();
  };

  class MyError : public std::exception {
  public:
    template<typename ... Ts> MyError(const Ts& ... args) { // en.cppreference.com/w/cpp/language/fold
      cout << "\n";
      cout << "c ******************************************************************\n";
      cout << "c MY_ERROR: ";
      (cout << ... << args); // fold expression
      cout << "\n";
      cout << "c ******************************************************************\n";
    }
  };

  vector<Int> getSortedNums(const Set<Int>& nums);

  TimePoint getTimePoint();
  Float getDuration(TimePoint start); // in seconds

  vector<string> splitInputLine(const string& line);
  

  template<typename T, typename U> pair<U, T> flipPair(const pair<T, U>& p) {
    return pair<U, T>(p.second, p.first);
  }

  template<typename T, typename U> multimap<U, T, greater<U>> flipMap(const Map<T, U>& inMap) { // decreasing
    multimap<U, T, greater<U>> outMap;
    transform(inMap.begin(), inMap.end(), inserter(outMap, outMap.begin()), flipPair<T, U>);
    return outMap;
  }

  template<typename T, typename U> bool isFound(const T& element, const vector<U>& container) {
    return std::find(begin(container), end(container), element) != end(container);
  }

  template<typename T> Set<T> getIntersection(const Set<T>& container1, const Set<T>& container2) {
    Set<T> intersection;
    for (const T& member : container1) {
      if (container2.contains(member)) {
        intersection.insert(member);
      }
    }
    return intersection;
  }

  template<typename T, typename U> Set<T> getDiff(const Set<T>& members, const U& nonMembers) {
    Set<T> diff;
    for (const T& member : members) {
      if (!nonMembers.contains(member)) {
        diff.insert(member);
      }
    }
    return diff;
  }

  template<typename T, typename U> void unionize(Set<T>& unionSet, const U& container) {
    for (const auto& member : container) {
      unionSet.insert(member);
    }
  }

  template<typename T> Set<T> getUnion(const vector<Set<T>>& containers) {
    Set<T> s;
    for (const Set<T>& container : containers) {
      unionize(s, container);
    }
    return s;
  }

  template<typename T> bool isDisjoint(const Set<T>& container1, const Set<T>& container2) {
    for (const T& member : container1) {
      if (container2.contains(member)) {
        return false;
      }
    }
    return true;
  }
}
