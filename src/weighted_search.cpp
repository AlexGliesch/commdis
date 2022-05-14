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
#include "weighted_search.h"
#include "local_improvement.h"

weighted_search::weighted_search(const solution& s) {
  cnd.resize(ndistr);
  for (int k = 0; k < ndistr; ++k)
    compute_move_of_lot(s, {}, k);
}

void weighted_search::choose_next_move(const solution& s, const tabu_list& tabu, int& u,
                                       int& k) const {
  (void)tabu;
  u = k = -1;
  double k_val = DBL_MAX;

  reservoir_sampling rs;
  for (int i = 0; i < ndistr; ++i) {
    if (cnd[i] == -1) continue;
    double val = move_value(s, cnd[i], i);
    if (k != -1 and f_eq(val, k_val)) {
      if (rs.consider()) k = i, k_val = val;
    } else if (k == -1 or val < k_val) {
      k = i, k_val = val, rs.reset();
    }
  }

  if (k == -1) {
    u = -1;
    return;
  }

  u = cnd[k];

  if (do_ls and f_geq(k_val, value(s))) k = u = -1;
}

void weighted_search::update(const solution& s, const tabu_list& tabu, int u, int p,
                             int k) {
  static unordered_set<int> must_upd;
  must_upd = {p, k};

  if (not f_eq(s.cmp(), s.cmp_cost(u, p))) {
    for (int i = 0; i < ndistr; ++i)
      compute_move_of_lot(s, tabu, i);
    statistic("mix. cand. recomp.", ndistr);
    return;
  }

  tabu_q.push(u);
  while (tabu_q.size() and not tabu.is_tabu(tabu_q.front())) {
    int v = tabu_q.front();
    tabu_q.pop();
    for (int w : nodes[v].nbs)
      must_upd.insert(s.assigned[w]);
  }

  for (int i : {k, p})
    for (int v : s.lot_nbs[i])
      if (s.assigned[v] != -1) must_upd.insert(s.assigned[v]);

  for (int i : must_upd)
    if (i != -1) compute_move_of_lot(s, tabu, i);
  statistic("mix. cand. recomp.", must_upd.size());
}

void weighted_search::compute_move_of_lot(const solution& s, const tabu_list& tabu,
                                          int k) {
  assert(k >= 0 and k < ndistr);
  int u = -1;
  double u_val = -1;
  reservoir_sampling rs;

  auto& nbs = s.get_neighbors(k, true);
  for (int v : nbs) {
    if (tabu.is_tabu(v)) continue;
    if (not s.move_keeps_connectivity(v, k)) continue;

    double v_val = move_value(s, v, k);
    if (u != -1 and f_eq(v_val, u_val)) {
      if (rs.consider()) u = v, u_val = v_val;
    } else if (u == -1 or v_val < u_val) {
      u = v, u_val = v_val, rs.reset();
    }
  }
  cnd[k] = u;
}
