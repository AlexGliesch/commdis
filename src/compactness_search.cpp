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
#include "compactness_search.h"
#include "local_improvement.h"

void compactness_search::choose_next_move(solution& s, const tabu_list& tabu,
                                          int& u_ret, int& k_ret,
                                          int already_did_pi) const {
  double best_cost = do_ls ? s.cmp() : DBL_MAX;
  u_ret = k_ret = -1;
  for (int i = 0; i < (int)s.diam.far.size(); ++i) {
    int u = s.diam.far[i];
    assert(u >= 0 and u < nnodes);

    if (i != (find(begin(s.diam.far), end(s.diam.far), u) - begin(s.diam.far)))
      continue;

    if (tabu.is_tabu(u)) continue;
    auto& nbs = nodes[u].nbs;
    for (auto j = begin(nbs); j != end(nbs); ++j) {
      int n = *j, k = s.assigned[n];
      if (j != (find_if(begin(nbs), end(nbs),
                        [&](int x) { return s.assigned[x] == k; })))
        continue;

      if (s.assigned[u] != k) {
        if (not s.move_keeps_connectivity(u, k)) continue;
        auto cost = s.cmp_cost(u, k);
        if (f_less(cost, best_cost)) {
          best_cost = cost;
          u_ret = u;
          k_ret = k;
        }
      }
    }
  }

  if (do_pi)
    if (u_ret == -1 and k_ret == -1 and already_did_pi < pi_max) {
      double val_before = s.cmp();
      compactness_path_improve(s);
      statistic("pi improves", f_less(s.cmp(), val_before));
      choose_next_move(s, tabu, u_ret, k_ret, already_did_pi + 1);
    }
}

void compactness_search::compactness_path_improve(solution& s) {
  using namespace boost::heap;
  static vector<double> dist;
  static vector<int> parent;
  static fibonacci_heap<pair<double, int>> pq;
  static vector<decltype(pq)::handle_type> handle;
  static vector<int> path;
  double cost = DBL_MAX;

  for (int i = 0; i < (int)s.diam.far.size(); ++i) {
    int u = s.diam.far[i];
    if (find(begin(s.diam.far), end(s.diam.far), u) - begin(s.diam.far) != i)
      continue;

    dist.assign(nnodes, -1);
    parent.assign(nnodes, -1);
    handle.resize(nnodes);
    pq.clear();
    dist[u] = 0;
    handle[u] = pq.emplace(0, u);
    while (pq.size()) {
      double d = -pq.top().first;
      int v = pq.top().second;
      pq.pop();

      if (s.assigned[v] != s.assigned[u]) {
        if (cost > d) {
          path.clear();
          cost = d;
          while (v != -1) {
            path.push_back(v);
            v = parent[v];
          }
        }
        break;
      }
      for (int n : nodes[v].nbs) {
        double cd = d + ::dist(n, v);
        if (dist[n] == -1) {
          handle[n] = pq.emplace(-cd, n);
          dist[n] = cd, parent[n] = v;
        } else if (dist[n] > cd) { 
          pq.update(handle[n], make_pair(-cd, n));
          dist[n] = cd, parent[n] = v;
        }
      }
    }
  }
  assert(path.size() >= 2);
  assert(count(begin(s.diam.far), end(s.diam.far), path.back()) > 0);
  double cmp_before = s.cmp();

  auto k1 = s.assigned[path.back()];
  assert(s.assigned[path[0]] != k1);
  for (int i = (int)path.size() - 1; i >= 1; --i) {
    assert(s.assigned[path[i]] == k1);
    s.assign_node(path[i], s.assigned[path[0]], true, false, false);
  }
  bool broke_connectivity = not s.lot_is_connected_bfs(k1);
  statistic("pi breaks connectivity", broke_connectivity);
  if (broke_connectivity or (do_ls and f_geq(s.cmp(), cmp_before))) {
    for (int i = 1; i < (int)path.size(); ++i)
      s.assign_node(path[i], k1, true, false, false);
  } else {
    s.recompute_all(true, false);
  }
  assert_correctness(s, true);
}
