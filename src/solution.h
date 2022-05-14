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
#include "dynamic_diameter.h"
#include "global.h"
#include <fmt/ostream.h>

struct solution {
  solution() { init_empty(); }

  inline bool operator<(const solution& s) const {
    if (is_incomplete() or s.is_incomplete()) return num_assigned > s.num_assigned;
    return lexi_less(balance(), cmp(), s.balance(), s.cmp());
  }

  bool is_incomplete() const { return num_assigned != nnodes; }

  void init_empty();

  void assign_initial_nodes(const vector<int>& initial_nodes);

  double balance_cost(int u, int k) const;

  inline double balance() const { return bal; }

  bool move_keeps_connectivity(int u, int k) const;

  void update_bridges(int k);

  vector<int>& get_neighbors(int k, bool nbs_ls) const;

  bool lot_is_connected_bfs(int k, int ignore = -1) const;

  void recompute_all(bool nbs_ls, bool recompute_nbh);

  void recompute_balance(bool populate_activities);

  dynamic_diameter diam;

  double cmp() const { return diam.diameter() + FAREPS * diam.count(); }

  double diameter() const { return diam.diameter(); }

  double cmp_cost(int u, int k) const {
    auto p = diam.move_cost(*this, u, k);
    return p.first + FAREPS * p.second;
  }

  void recompute_neighborhoods(bool nbs_ls);

  void assign_node(int u, int k, bool nbs_ls, bool upd_objs, bool upd_bridges);

  bool test_nb(int u, int k, bool nbs_ls) const {
    if (assigned[u] == k) return false;
    if (not nbs_ls and assigned[u] != -1) return false;
    assert(nbs_ls or assigned[u] == -1);
    for (int n : nodes[u].nbs)
      if (assigned[n] == k) return true;
    return false;
  }

  inline void add_nb(int u, int k) {
    if (find(begin(lot_nbs[k]), end(lot_nbs[k]), u) == end(lot_nbs[k]))
      lot_nbs[k].push_back(u);
  }

  inline int area(int k) const {
    return k >= 0 and k < ndistr ? lot_assigned[k].size() : -1;
  }

  uint64_t hash() const;

  double bal = -1;
  int num_assigned = 0;
  array<vector<double>, 3> w;
  vector<int> assigned;
  vector<vector<int>> lot_nbs;
  vector<vector<int>> lot_assigned;
  vector<bool> is_bridge;
};

inline ostream& operator<<(ostream& o, const solution& s) {
  o << setprecision(4) << fixed;
  o << "(cmp: " << (s.is_incomplete() ? -1.0 : s.diam.diameter());
  o << ", bal: " << (s.is_incomplete() ? -1.0 : s.balance()) << ")";
  return o;
}

void assert_correctness(const solution& s, bool ls);
