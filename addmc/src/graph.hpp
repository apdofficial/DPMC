#pragma once
#include "types.hpp"
#include <vector>
using std::pair;
using std::vector;
namespace dpve{
class Graph { // undirected
public:
  Set<Int> vertices;
  Map<Int, Set<Int>> adjacencyMap;

  Graph(const Set<Int>& vs);
  void addEdge(Int v1, Int v2);

  bool isNeighbor(Int v1, Int v2) const;
  bool hasPath(Int from, Int to, Set<Int>& visitedVertices) const; // path length >= 0
  bool hasPath(Int from, Int to) const;
  void removeVertex(Int v); // also removes edges from and to `v`
  void fillInEdges(Int v); // does not remove `v`
  Int getFillInEdgeCount(Int v) const;
  Int getMinFillVertex() const;
  Graph projectOnto(Set<Int> vars) const;
};

class Label : public vector<Int> { // for lexicographic search
public:
  void addNumber(Int i); // retains descending order
  static bool hasSmallerLabel(const pair<Int, Label>& a, const pair <Int, Label>& b);
};
} //end namespace dpve