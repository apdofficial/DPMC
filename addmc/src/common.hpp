#pragma once

/* inclusions =============================================================== */

#include "types.hpp"

#include <cassert>
#include <iostream>
#include <map>

using std::map;
using std::ostream;

namespace dpve{
   /* diagram packages: */
  inline const string CUDD_PACKAGE = "c";   //inline after c++17 allows constant definition in header without violating odr
  inline const string SYLVAN_PACKAGE = "s";
  inline const map<string, string> DD_PACKAGES = {
    {CUDD_PACKAGE, "CUDD"},
    {SYLVAN_PACKAGE, "SYLVAN"}
  };

  /* CNF var order heuristics: */
  inline const Int RANDOM_HEURISTIC = 0;
  inline const Int DECLARATION_HEURISTIC = 1;
  inline const Int MOST_CLAUSES_HEURISTIC = 2;
  inline const Int MIN_FILL_HEURISTIC = 3;
  inline const Int MCS_HEURISTIC = 4;
  inline const Int LEX_P_HEURISTIC = 5;
  inline const Int LEX_M_HEURISTIC = 6;
  inline const Int COLAMD_HEURISTIC = 7;

  inline const map<Int, string> CNF_VAR_ORDER_HEURISTICS = {
    {RANDOM_HEURISTIC, "RANDOM"},
    {DECLARATION_HEURISTIC, "DECLARATION"},
    {MOST_CLAUSES_HEURISTIC, "MOST_CLAUSES"},
    {MIN_FILL_HEURISTIC, "MIN_FILL"},
    {MCS_HEURISTIC, "MCS"},
    {LEX_P_HEURISTIC, "LEX_P"},
    {LEX_M_HEURISTIC, "LEX_M"},
    {COLAMD_HEURISTIC, "COLAMD"}
  };

  /* JT var order heuristics: */
  inline const Int BIGGEST_NODE_HEURISTIC = 8;
  inline const Int HIGHEST_NODE_HEURISTIC = 9;
  inline const map<Int, string> JOIN_TREE_VAR_ORDER_HEURISTICS = {
    {BIGGEST_NODE_HEURISTIC, "BIGGEST_NODE"},
    {HIGHEST_NODE_HEURISTIC, "HIGHEST_NODE"}
  };

  /* clustering heuristics: */
  inline const string BUCKET_ELIM_LIST = "bel";
  inline const string BUCKET_ELIM_TREE = "bet";
  inline const string BOUQUET_METHOD_LIST = "bml";
  inline const string BOUQUET_METHOD_TREE = "bmt";
  inline const map<string, string> CLUSTERING_HEURISTICS = {
    {BUCKET_ELIM_LIST, "BUCKET_ELIM_LIST"},
    {BUCKET_ELIM_TREE, "BUCKET_ELIM_TREE"},
    {BOUQUET_METHOD_LIST, "BOUQUET_METHOD_LIST"},
    {BOUQUET_METHOD_TREE, "BOUQUET_METHOD_TREE"}
  };

  /* maximizer formats: */
  inline const Int NEITHER_FORMAT = 0;
  inline const Int SHORT_FORMAT = 1;
  inline const Int LONG_FORMAT = 2;
  inline const Int DUAL_FORMAT = 3;
  inline const map<Int, string> MAXIMIZER_FORMATS = {
    {NEITHER_FORMAT, "NEITHER"},
    {SHORT_FORMAT, "SHORT"},
    {LONG_FORMAT, "LONG"},
    {DUAL_FORMAT, "DUAL"}
  };

  /* join priorities: */
  inline const string ARBITRARY_PAIR = "a";
  inline const string BIGGEST_PAIR = "b";
  inline const string SMALLEST_PAIR = "s";
  inline const string FCFS = "f";
  inline const map<string, string> JOIN_PRIORITIES = {
    {ARBITRARY_PAIR, "ARBITRARY_PAIR"},
    {BIGGEST_PAIR, "BIGGEST_PAIR"},
    {SMALLEST_PAIR, "SMALLEST_PAIR"},
    {FCFS, "FCFS"}
  };
/* global functions ========================================================= */

ostream& operator<<(ostream& stream, const dpve::Number& n);
} //end namespace dpve