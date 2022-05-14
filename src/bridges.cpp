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
#include "bridges.h"
#include "global.h"
#include "solution.h"

namespace {
int root, root_children, counter;
vector<int> low, num, parent;

void find_bridges_rec(solution& s, int u) {
  low[u] = num[u] = counter++;
  for (int v : nodes[u].nbs) {
    if (s.assigned[v] != s.assigned[root]) continue;
    if (num[v] == -1) {
      parent[v] = u;
      if (u == root) ++root_children;
      find_bridges_rec(s, v);
      if (low[v] >= num[u]) s.is_bridge[u] = true;
      low[u] = min(low[u], low[v]);
    } else if (v != parent[u]) {
      low[u] = min(low[u], num[v]);
    }
  }
}
} // unnamed namespace

void find_bridges(solution& s, int k) {
  int V = s.lot_assigned[k].size();
  assert(V >= 1);
  assert((int)s.is_bridge.size() == nnodes);
  num.assign(nnodes, -1);
  low.assign(nnodes, 0);
  parent.assign(nnodes, 0);
  for (int u : s.lot_assigned[k]) s.is_bridge[u] = false;
  root = s.lot_assigned[k][0];
  root_children = counter = 0;
  find_bridges_rec(s, root);
  s.is_bridge[root] = (V == 1 or root_children > 1);
}