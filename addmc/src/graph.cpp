#include "graph.hpp"
#include "util.hpp"

using dpve::Graph;
using dpve::Label;
using dpve::Int;

using std::next;
/* class Graph ============================================================== */

Graph::Graph(const Set<Int>& vs) {
  for (Int v : vs) {
    vertices.insert(v);
    adjacencyMap[v] = Set<Int>();
  }
}

void Graph::addEdge(Int v1, Int v2) {
  adjacencyMap.at(v1).insert(v2);
  adjacencyMap.at(v2).insert(v1);
}

bool Graph::isNeighbor(Int v1, Int v2) const {
  return adjacencyMap.at(v1).contains(v2);
}

bool Graph::hasPath(Int from, Int to, Set<Int>& visitedVertices) const {
  if (from == to) {
    return true;
  }

  visitedVertices.insert(from);

  Set<Int> unvisitedNeighbors = dpve::util::getDiff(adjacencyMap.at(from), visitedVertices);

  for (Int v : unvisitedNeighbors) {
    if (hasPath(v, to, visitedVertices)) {
      return true;
    }
  }
  return false;
}

bool Graph::hasPath(Int from, Int to) const {
  Set<Int> visitedVertices;
  return hasPath(from, to, visitedVertices);
}

void Graph::removeVertex(Int v) {
  vertices.erase(v);

  adjacencyMap.erase(v); // edges from v

  for (pair<const Int, Set<Int>>& vertexAndNeighbors : adjacencyMap) {
    vertexAndNeighbors.second.erase(v); // edge to v
  }
}

void Graph::fillInEdges(Int v) {
  const Set<Int>& neighbors = adjacencyMap.at(v);
  for (auto neighbor1 = neighbors.begin(); neighbor1 != neighbors.end(); neighbor1++) {
    for (auto neighbor2 = next(neighbor1); neighbor2 != neighbors.end(); neighbor2++) {
      addEdge(*neighbor1, *neighbor2);
    }
  }
}

Int Graph::getFillInEdgeCount(Int v) const {
  Int edgeCount = 0;
  const Set<Int>& neighbors = adjacencyMap.at(v);
  for (auto neighbor1 = neighbors.begin(); neighbor1 != neighbors.end(); neighbor1++) {
    for (auto neighbor2 = next(neighbor1); neighbor2 != neighbors.end(); neighbor2++) {
      if (!isNeighbor(*neighbor1, *neighbor2)) {
        edgeCount++;
      }
    }
  }
  return edgeCount;
}

Int Graph::getMinFillVertex() const {
  Int vertex = dpve::MIN_INT;
  Int fillInEdgeCount = dpve::MAX_INT;

  for (Int v : vertices) {
    Int edgeCount = getFillInEdgeCount(v);
    if (edgeCount < fillInEdgeCount) {
      fillInEdgeCount = edgeCount;
      vertex = v;
    }
  }

  if (vertex == dpve::MIN_INT) {
    throw dpve::util::MyError("graph has no vertex");
  }

  return vertex;
}

Graph Graph::projectOnto(Set<Int> vars) const{
  Graph graph(vars);
  for (auto& var1 : vars){
    for (auto& var2: vars){
      if (isNeighbor(var1,var2)){
        graph.addEdge(var1,var2);
      }
    }
  }
  return graph;
}

/* class Label ============================================================== */

void Label::addNumber(Int i) {
  push_back(i);
  sort(begin(), end(), greater<Int>());
}

bool Label::hasSmallerLabel(const pair<Int, Label>& a, const pair <Int, Label>& b) {
  return a.second < b.second;
}

