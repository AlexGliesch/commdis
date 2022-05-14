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
#include "solution.h"
#include "bridges.h"

void solution::init_empty() {
  assigned.assign(nnodes, -1);
  lot_nbs.assign(ndistr, {});
  lot_assigned.assign(ndistr, {});
  is_bridge.assign(nnodes, false);
  num_assigned = 0;
  diam.init(*this);
  recompute_balance(true);
}

void solution::assign_initial_nodes(const vector<int>& initial_nodes) {
  init_empty();
  assert((int)initial_nodes.size() == ndistr);

  for (int i = 0; i < ndistr; ++i)
    assign_node(initial_nodes[i], i, false, true, false);
#ifndef NDEBUG
  for (auto& i : lot_assigned) assert(i.size() == 1);
  assert(num_assigned == ndistr);
#endif
}

double solution::balance_cost(int u, int k) const {
  int p = assigned[u];
  assert(p != k);
  double val = bal;
  for3(j) {
    // TODO simplify here (see paper)
    val -= max(0.0, abs(w[j][k] - mu[j]) - mu_tau[j]) / mu[j];
    val += max(0.0, abs(w[j][k] + nodes[u].a[j] - mu[j]) - mu_tau[j]) / mu[j];
    if (p != -1) {
      val -= max(0.0, abs(w[j][p] - mu[j]) - mu_tau[j]) / mu[j];
      val += max(0.0, abs(w[j][p] - nodes[u].a[j] - mu[j]) - mu_tau[j]) / mu[j];
    }
  }
#ifdef DEBUG_BALANCE
  auto rw = w;
  if (p != -1) for3(i) rw[i][p] -= nodes[u].a[i];
  for3(i) rw[i][k] += nodes[u].a[i];

  double rbal = 0;
  for (int k = 0; k < ndistr; ++k)
    for3(j) rbal += max(0.0, abs(rw[j][k] - mu[j]) - mu_tau[j]) / mu[j];

  assert(f_eq(rbal, val));
#endif
  return val;
}

bool solution::move_keeps_connectivity(int u, int k) const {
  if (not test_nb(u, k, true) or area(assigned[u]) <= 1) return false;
#ifdef DEBUG_BRIDGES
  assert(not use_bridges or
         lot_is_connected_bfs(assigned[u], u) == (not is_bridge[u]));
#endif
  return use_bridges ? not is_bridge[u] : lot_is_connected_bfs(assigned[u], u);
}

void solution::update_bridges(int k) {
  if (use_bridges) {
    find_bridges(*this, k);
#ifdef DEBUG_BRIDGES
    for (int u : lot_assigned[k])
      assert(is_bridge[u] - lot_is_connected_bfs(k, u));
#endif
  }
}

vector<int>& solution::get_neighbors(int k, bool nbs_ls) const {
  auto& s = const_cast<solution&>(*this);
  auto& nb = s.lot_nbs[k];
  static vector<int> ans;
  ans.resize(nb.size());
  int j = 0;
  for (int i = 0; i < (int)nb.size();) {
    if ((nbs_ls and s.test_nb(nb[i], k, true)) or s.assigned[nb[i]] == -1) {
      ans[j++] = nb[i]; // .push_back(nb[i]);
      ++i;
    } else {
      swap(nb[i], nb.back());
      s.lot_nbs[k].pop_back();
    }
  }
  ans.resize(j);
#ifdef DEBUG_SOLUTION_DS
  assert(count(begin(ans), end(ans), -1) == 0);
#endif
  return ans;
}

bool solution::lot_is_connected_bfs(int k, int ignore) const {
  for (int v : lot_assigned[k])
    if (v != ignore) {
      static queue<int> q;
      static vector<bool> visited;
      assert(q.empty());
      visited.assign(nnodes, false);
      visited[v] = true;
      q.push(v);
      int num_visited = 1;
      while (q.size()) {
        int w = q.front();
        q.pop();
        for (int n : nodes[w].nbs) 
          if (n != ignore and assigned[n] == k and not visited[n]) {
            visited[n] = true;
            ++num_visited;
            q.push(n);
          }
      }
      return num_visited == (int)lot_assigned[k].size() - (ignore != -1);
    }
  return false;
}

void solution::recompute_balance(bool populate_activities) {
  if (populate_activities) {
    for3(j) w[j].assign(ndistr, 0);
    auto& a = assigned;
    for (int i = 0; i < nnodes; ++i)
      if (a[i] != -1) for3(j) w[j][a[i]] += nodes[i].a[j];
  }

  bal = 0;
  for (int k = 0; k < ndistr; ++k)
    for3(j) bal += max(0.0, abs(w[j][k] - mu[j]) - mu_tau[j]) / mu[j];

  if (f_eq(bal, 0.0)) bal = 0.0;
}

void solution::recompute_neighborhoods(bool nbs_ls) {
  num_assigned = 0;
  for (int k = 0; k < ndistr; ++k) {
    lot_nbs[k].clear();
    lot_assigned[k].clear();
  }
  for (int u = 0; u < nnodes; ++u) {
    if (assigned[u] != -1) {
      lot_assigned[assigned[u]].push_back(u);
      ++num_assigned;
    }
    for (int k = 0; k < ndistr; ++k)
      if (test_nb(u, k, nbs_ls)) add_nb(u, k);
  }
}

uint64_t solution::hash() const {
  return boost::hash_range(begin(assigned), end(assigned));
}

void solution::assign_node(int u, int k, bool nbs_ls, bool upd_objs,
                           bool upd_bridges) {
  int p = assigned[u];
  assert(nbs_ls or assigned[u] == -1);
  assert(p != k);
  if (p == -1)
    ++num_assigned; 
  if (p != -1) {
    auto& lap = lot_assigned[p];
    auto f = find(begin(lap), end(lap), u) - begin(lap);
    assert(f >= 0 and f < (int)lap.size());
    assert(not lap.empty());
    swap(lap[f], lap.back());
    lap.pop_back(); 
  }

  if (upd_objs) {
    bal = balance_cost(u, k);
    if (f_eq(bal, 0.0)) bal = 0.0;
  }

  assigned[u] = k;
  lot_assigned[k].push_back(u);

  if (upd_objs) {
    if (p != -1) for3(i) w[i][p] -= nodes[u].a[i];
    for3(i) w[i][k] += nodes[u].a[i];
    diam.update(*this, u, k, p);
  }

  for (int n : nodes[u].nbs) {
    if ((nbs_ls and assigned[n] != k) or assigned[n] == -1) add_nb(n, k);
    if (nbs_ls and p != -1)
      for (int nn : nodes[n].nbs)
        if (test_nb(nn, p, nbs_ls)) add_nb(nn, p);
  }

  if (upd_bridges and use_bridges) {
    if (p != -1) update_bridges(p);
    update_bridges(k);
  }
}

void solution::recompute_all(bool nbs_ls, bool recompute_nbh) {
  num_assigned = 0;
  lot_assigned.assign(ndistr, {});
  for (int i = 0; i < nnodes; ++i)
    if (assigned[i] != -1) {
      ++num_assigned;
      lot_assigned[assigned[i]].push_back(i);
    }
  if (recompute_nbh) recompute_neighborhoods(nbs_ls);
  recompute_balance(true);
  diam.recompute(*this);
  for (int k = 0; k < ndistr; ++k) update_bridges(k);
  assert_correctness(*this, nbs_ls);
}

void assert_correctness(const solution& s, bool ls) {
  (void)s;
  (void)ls;
#ifdef DEBUG_SOLUTION_DS
  for (int i = 0; i < ndistr; ++i) assert(s.lot_is_connected_bfs(i));

  for (int k = 0; k < ndistr; ++k) {
    static set<int> s1, s2;
    s1.clear(), s2.clear();
    for (int n : s.lot_nbs[k]) s1.insert(n);

    for (int i = 0; i < nnodes; ++i)
      if (s.assigned[i] == k)
        for (int n : nodes[i].nbs)
          if (s.assigned[n] == -1 or (ls and s.assigned[n] != k)) s2.insert(n);

    assert(includes(begin(s1), end(s1), begin(s2), end(s2)));
  }

  static vector<vector<int>> lot_assigned;
  lot_assigned.assign(ndistr, {});
  for (int i = 0; i < nnodes; ++i)
    if (s.assigned[i] != -1) lot_assigned[s.assigned[i]].push_back(i);

  for (int k = 0; k < ndistr; ++k) {
    static vector<int> sla;
    sla = s.lot_assigned[k];
    sort(begin(sla), end(sla));
    assert(sla == lot_assigned[k]);
  }
#endif

#ifdef DEBUG_BALANCE
  vector<double> w[3];
  for3(i) w[i].assign(ndistr, 0);
  for (int i = 0; i < nnodes; ++i) {
    if (s.assigned[i] != -1) {
      for3(j) w[j][s.assigned[i]] += nodes[i].a[j];
    }
  }
  for (int i = 0; i < ndistr; ++i) for3(j) assert(f_eq(w[j][i], s.w[j][i]));
#endif

#ifdef DEBUG_DIAMETER
  double md = -1;
  int count = 0;
  for (int k = 0; k < ndistr; ++k) {
    auto& la = s.lot_assigned[k];
    for (int i = 0; i < (int)la.size(); ++i)
      for (int j = i; j < (int)la.size(); ++j)
        consider_dist(md, count, dist(la[i], la[j]), 1);
  }
  assert(f_eq(FAREPS * count + md, s.cmp()));
#endif
}
