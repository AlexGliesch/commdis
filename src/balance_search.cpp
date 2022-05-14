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
#include "balance_search.h"
#include "global.h"
#include "local_improvement.h"

balance_search::balance_search(solution& s, double max_diam) {
  cnd.resize(ndistr);
  this->max_diam = max_diam;
  recompute_moves(s, {});
}

void balance_search::choose_next_move(const solution& s, const tabu_list& tabu, int& u,
                                      int& k) const {
  (void)tabu;
  u = k = -1;
  if (f_eq(s.balance(), 0.0)) return;

  reservoir_sampling rs;
  for (int i = 0; i < ndistr; ++i) {
    if (cnd[i].u == -1) continue;

    if (k != -1 and f_eq(cnd[i].bal, cnd[k].bal)) {
      if (rs.consider()) k = i;
    } else if (k == -1 or cnd[i].bal < cnd[k].bal)
      rs.reset(), k = i;
  }

  if (k == -1) {
    u = -1;
    return;
  }

  u = cnd[k].u;

  if (do_ls and f_geq(cnd[k].bal, 0)) u = k = -1;

  assert(k >= -1 and k < ndistr);
  assert(u >= -1 and u < nnodes);
  if (u != -1 and k != -1) {
    assert(k != s.assigned[u]);
    assert(not tabu.is_tabu(u));
    assert(f_eq(s.balance_cost(u, k) - s.balance(), cnd[k].bal));
  }
  statistic("lot index chosen", k);
}

void balance_search::update(const solution& s, const tabu_list& tabu, int u, int p,
                            int k) {
  assert(s.assigned[u] == k);
  assert(k != p);
  assert(k != -1 and p != -1);

  static unordered_set<int> must_upd;
  must_upd = {p, k};

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

#ifdef DEBUG_BAL_NBH
  for (int i = 0; i < ndistr; ++i) {
    auto& c = cnd[i];
    if (c.u == -1) continue;
    if (must_upd.count(i) == 0) {
      assert(s.test_nb(c.u, i, true));
      assert(s.assigned[c.u] != i);
      auto cc = c;
      compute_move_of_lot(s, tabu, i);
      assert(f_eq(s.balance_cost(cc.u, i), s.balance_cost(c.u, i)));
    }
    if (c.u == u and i != k) {
      auto& nbck = s.lot_nbs[i];
      assert(must_upd.count(i));
      assert(find(begin(nbck), end(nbck), u) != end(nbck));
    }
    if (not s.move_keeps_connectivity(c.u, i)) assert(must_upd.count(i));
  }
#endif
  for (int i : must_upd)
    if (i != -1) compute_move_of_lot(s, tabu, i);
  statistic("bal. cand. recomp.", must_upd.size());
}

void balance_search::recompute_moves(const solution& s, const tabu_list& tabu) {
  for (int k = 0; k < ndistr; ++k)
    compute_move_of_lot(s, tabu, k);
}

void balance_search::compute_move_of_lot(const solution& s, const tabu_list& tabu,
                                         int k) {
  assert(k >= 0 and k < ndistr);
  auto& b = cnd[k];
  b.u = b.bal = -1;
  reservoir_sampling rs;

  auto& nbs = s.get_neighbors(k, true);
  for (int v : nbs) {
    if (tabu.is_tabu(v)) continue;
    if (not s.move_keeps_connectivity(v, k)) continue;

    cand_move cc;
    cc.u = v;
    double cmp = s.diam.max_dist_u_in_k(s, v, k).first;
    statistic("discard cmp", f_greater(cmp, max_diam));
    if (f_greater(cmp, max_diam)) continue;

    cc.bal = s.balance_cost(v, k) - s.balance();
    if (cc.u != -1 and f_eq(cc.bal, b.bal)) {
      if (rs.consider()) b = cc;
    } else if (b.u == -1 or cc.bal < b.bal) {
      rs.reset(), b = cc;
    } else {
    }
  }

  if (f_eq(b.bal, 0)) b.bal = 0;

#ifdef DEBUG_BAL_NBH
  if (b.u != -1) {
    assert(not tabu.is_tabu(b.u));
    assert(s.assigned[b.u] != k);
    assert(s.test_nb(b.u, k, true));
    assert(s.assigned[b.u] != -1);
    assert(s.move_keeps_connectivity(b.u, k));
  }
#endif
}
