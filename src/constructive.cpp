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
#include "constructive.h"
#include "balance_search.h"
#include "solution.h"

void random_initial_nodes(solution& s) {
  assert(s.num_assigned == 0);
  static vector<int> indices, v;
  indices.resize(nnodes);
  iota(begin(indices), end(indices), 0);
  choice(indices, v, ndistr);
  assert((int)v.size() == ndistr);
  s.assign_initial_nodes(v);
}

void max_disp_initial_nodes(solution& s) {
  assert(s.num_assigned == 0);
  vector<int> centers(ndistr);

  centers[0] = rand_int(0, nnodes - 1);

  s.assign_node(centers[0], 0, false, true, false);

  static vector<int> indices, v;
  indices.resize(nnodes);
  iota(begin(indices), end(indices), 0);
  remove_from_sorted(indices, centers[0]);

  for (int k = 1; k < ndistr; ++k) {
    int best_cand = -1;
    double max_min_dist = -1;
    choice(indices, v, max(2, (int)sqrt(nnodes)));
    for (int j : v) {
      assert(s.assigned[j] == -1);
      double min_dist = DBL_MAX;
      for (int l = 0; l < k; ++l)
        min_dist = min(min_dist, dist2(j, centers[l]));
      if (min_dist > max_min_dist) {
        max_min_dist = min_dist;
        best_cand = j;
      }
    }

    assert(best_cand != -1);
    centers[k] = best_cand;

    s.assign_node(centers[k], k, false, true, false);
    remove_from_sorted(indices, centers[k]);
  }
}

void fill_greedy(solution& s) {
  vector<int> cnd(ndistr, -1);
  vector<double> cnd_cmp(ndistr, -1);

  auto recompute_cand = [&](int k) {
    int& u = cnd[k];
    double& cmp = cnd_cmp[k];
    u = cmp = -1;
    reservoir_sampling rs;
    auto& nbs = s.get_neighbors(k, false);

    for (int n : nbs) {
      assert(s.assigned[n] == -1);
      double n_cmp;
      if (cons_alg.back() == 'd') { // "greedy-d"
        n_cmp = max(s.diam.max_dist_u_in_k(s, n, k).first, s.diam.lot_diameter(k));
      } else /*if (opt == B)*/ {
        n_cmp = 0;
        for3(j) n_cmp +=
            max(0.0, abs(s.w[j][k] + nodes[n].a[j] - mu[j]) - mu_tau[j]) / mu[j];
      }
      if (u != -1 and f_eq(n_cmp, cmp)) {
        if (rs.consider()) u = n;
      } else if (u == -1 or n_cmp < cmp) {
        rs.reset(), u = n, cmp = n_cmp;
      }
    }
    assert(u == -1 or s.assigned[u] == -1);
    if (cons_alg.back() == 'b') cmp = -cmp;
  };

  for (int k = 0; k < ndistr; ++k) {
    assert(s.area(k) == 1);
    recompute_cand(k);
  }

  while (s.num_assigned != nnodes) {
    int k = -1;

    double sum_k = -1;
    for (int i = 0; i < ndistr; ++i) {
      assert(cnd[i] == -1 or s.assigned[cnd[i]] == -1);
      if (cnd[i] != -1) {
        double sum_i = 0;
        for3(j) sum_i += s.w[j][i] / mu[j];
        if (k == -1 or lexi_less(sum_i, cnd_cmp[i], sum_k, cnd_cmp[k]))
          k = i, sum_k = sum_i;
      }
    }

    assert(k != -1);
    assert(cnd[k] != -1);
    assert(s.assigned[cnd[k]] == -1);

    s.assign_node(cnd[k], k, false, true, false);

    for (int i = 0; i < ndistr; ++i)
      if (s.assigned[cnd[i]] != -1) recompute_cand(i);
  }
}

void fill_bfs(solution& s) {
  queue<int> q;
  for (int k = 0; k < ndistr; ++k) {
    assert(s.lot_assigned[k].size() == 1);
    for (int n : nodes[s.lot_assigned[k][0]].nbs)
      q.push(n);
  }
  while (q.size()) {
    int u = q.front();
    q.pop();
    if (s.assigned[u] != -1) continue;
    int k = -1;
    for (int n : nodes[u].nbs) {
      if (s.assigned[n] == -1)
        q.push(n);
      else
        k = s.assigned[n];
    }
    assert(k != -1);
    s.assign_node(u, k, false, false, false);
  }
  assert(s.num_assigned == nnodes);
}

solution constructive() {
  solution s;
  if (seeds_alg == "random")
    random_initial_nodes(s);
  else if (seeds_alg == "maxdisp")
    max_disp_initial_nodes(s);
  else
    assert(false);

  if (cons_alg == "bfs")
    fill_bfs(s);
  else if (cons_alg == "greedy-d")
    fill_greedy(s);
  else if (cons_alg == "greedy-b")
    fill_greedy(s);
  else
    assert(false);

  boost::hash_combine(execution_hash, s.balance() + s.diameter());
  return (s);
}
