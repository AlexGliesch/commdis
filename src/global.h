/*
* A multistart alternating tabu search for commercial districting
* Copyright (c) 2018 Alex Gliesch, Marcus Ritt, Mayron C. O. Moreira
*
* Permission is hereby granted, free of charge, to any person (the "Person")
* obtaining a copy of this software and associated documentation files (the
* "Software"), to deal in the Software, including the rights to use, copy, modify,
* merge, publish, distribute the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* 1. The above copyright notice and this permission notice shall be included in
*    all copies or substantial portions of the Software.
* 2. Under no circumstances shall the Person be permitted, allowed or authorized
*    to commercially exploit the Software.
* 3. Changes made to the original Software shall be labeled, demarcated or
*    otherwise identified and attributed to the Person.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
* FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
* IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
* CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#pragma once
#include "pre.h"
#define for3(a) for (int a = 0; a < 3; ++a)
// #define NDEBUG
// #define DEBUG_DIAMETER     // dynamic diameter
// #define DEBUG_BRIDGES      // dynamic connectivity
// #define DEBUG_BALANCE      // dynamic balance computation
// #define DEBUG_SOLUTION_DS  // data structures of "solution" class
// #define DEBUG_CONSTRUCTIVE // constructive algorithm
// #define DEBUG_GEOMETRY     // geometry
// #define DEBUG_BAL_NBH      // balance search
// #define DEBUG_MIXED_NBH    // weigthed obj. fun. search

extern vector<int> E; // List of edges

struct node {
  array<double, 3> a;
  double value;
  double x, y;

  struct {
    auto begin() const { return E.begin() + i0; }
    auto end() const { return E.begin() + i1; }
    int i0, i1;
  } nbs;

  size_t hash() const {
    size_t seed = 0;
    boost::hash_combine(seed, x);
    boost::hash_combine(seed, y);
    return seed;
  }
};

extern int time_limit_seconds;
extern int max_iter;
extern uint64_t random_seed;
extern string instance_name;
extern string out_file;
extern int tabu_tenure;
extern double tenure_mult;
extern int I_max;
extern int A_max;
extern int pi_max;
extern double alpha;
extern bool do_only_cons;
extern bool do_cmp_search;
extern bool do_bal_search;
extern bool do_pi;
extern bool do_ls;
extern bool do_weighted;
extern string seeds_alg;
extern string cons_alg;
extern bool use_bridges;
extern bool use_calipers;
extern bool bootstrap;
extern bool irace;

// Global variables
extern double tau[3]; 
extern array<double, 3> mu;
extern array<double, 3> mu_tau;
extern double d_max;           
extern int nnodes;             
extern int ndistr;              
extern vector<node> nodes;
extern int nedges;            
extern size_t execution_hash; 

inline double dist2(int a, int b) {
  double x1 = nodes[a].x, y1 = nodes[a].y;
  double x2 = nodes[b].x, y2 = nodes[b].y;
  return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

inline double dist(int a, int b) { return sqrt(dist2(a, b)); }

inline bool lexi_less(double a1, double a2, double b1, double b2) {
  return f_eq(a1, b1) ? (f_eq(a2, b2) ? false : a2 < b2) : a1 < b1;
}

inline bool lexi_leq(double a1, double a2, double b1, double b2) {
  return f_eq(a1, b1) ? (f_eq(a2, b2) ? true : a2 < b2) : a1 < b1;
}

inline bool timeout() { return elapsed_time(start_time) > time_limit_seconds; }

void statistic(const string& s, double value);

void print_statistics();

void read_instance();

void read_cmd_line(int argc, char** argv);
