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
#include "dynamic_diameter.h"
#include "solution.h"

void dynamic_diameter::init(const solution& s) {
  kmd = -1;
  recompute(s);
}

pair<double, int> dynamic_diameter::move_cost(const solution& s, int u,
                                              int k) const {
  assert(u >= 0 and u < nnodes);
  assert(k >= 0 and k < ndistr);
  int p = s.assigned[u];

  int u_is_far = ::count(begin(far), end(far), u);
  int count = far.size() / 2 - u_is_far;
  bool removing_u_will_reduce_diam = count == 0;
  double FX = removing_u_will_reduce_diam ? -1 : diams[kmd].max_diam;

  auto duk = max_dist_u_in_k(s, u, k);
  consider_dist(FX, count, duk.first, duk.second);
  assert(count != -1 and FX != -1);

  if (p != -1 and removing_u_will_reduce_diam) {
    diam_info diam_p = diams[p];
    assert(diam_p.k == p);
    diam_p.lose_point(s, u);
    consider_dist(FX, count, diam_p.max_diam, diam_p.far.size() / 2);

    for (int i = 0; i < ndistr; ++i)
      if (i != p /* and i != k*/)
        consider_dist(FX, count, diams[i].max_diam, diams[i].far.size() / 2);
  }

#ifdef DEBUG_DIAMETER
  auto r = s;
  r.assign_node(u, k, s.num_assigned == nnodes, false, false);
  r.diam.recompute(r);
  assert(count == (int)r.diam.count());
  assert(f_eq(r.cmp(), FAREPS * count + FX));
#endif
  return {FX, count};
}

pair<double, int> dynamic_diameter::max_dist_u_in_k(const solution& s, int u,
                                                    int k) const {
  double FX = -1;
  int count = 0;

  if (use_calipers and
      diams[k].P.size() > 3) { 
    for (const auto& X : {&diams[k].U, &diams[k].L}) {
      for (const auto& p : *X) consider_dist(FX, count, dist2(u, p.index), 1);
    }
#ifdef DEBUG_DIAMETER
    double FX2 = -1;
    int c2 = 0;
    for (int n : s.lot_assigned[k]) consider_dist(FX2, c2, dist2(u, n), 1);
    assert(f_eq(FX2, FX));
    assert(c2 == count);
#endif
  } else {
    for (int n : s.lot_assigned[k]) consider_dist(FX, count, dist2(u, n), 1);
  }
  assert(count >= 1);
  return {sqrt(FX), count};
}

void dynamic_diameter::recompute(const solution& s) {
  diams.resize(ndistr);
  for (int k = 0; k < ndistr; ++k) diams[k].k = k;

  kmd = 0;
  for (int k = 0; k < ndistr; ++k) {
    if (use_calipers) {
      diams[k].recompute_ABCD_and_P(s, k, -1);
    } else {
      diams[k].P.clear();
      for (int u : s.lot_assigned[k]) diams[k].P.emplace_back(u);
    }
    diams[k].comp_diam();
    assert(s.area(k) <= 1 or diams[k].far.size() >= 1);
  }
  update_kmd_and_far();

  assert(s.num_assigned <= ndistr or far.size() >= 1);
#ifdef DEBUG_DIAMETER
  for (int i = 0; i < (int)far.size(); i += 2)
    for (int j = i + 2; j < (int)far.size(); j += 2)
      assert(not((far[i] == far[j] and far[i + 1] == far[j]) or
                 (far[i] == far[j + 1] and far[i + 1] == far[j])));
#endif
}

void dynamic_diameter::update(const solution& s, int u, int k, int p) {
  assert(k != -1);
  assert(diams[k].k == k and (p == -1 or diams[p].k == p));

  if (p != -1) diams[p].lose_point(s, u);
  diams[k].gain_point(s, u);

#ifdef DEBUG_DIAMETER
  for (int i = 0; i < (int)diams[k].P.size(); ++i)
    for (int j = i + 1; j < (int)diams[k].P.size(); ++j)
      assert(diams[k].P[i] != diams[k].P[j]);
#endif
  update_kmd_and_far();
}

void dynamic_diameter::update_kmd_and_far() {
  kmd = 0, far = diams[0].far;
  for (int k = 1; k < (int)diams.size(); ++k) {
    if (f_eq(diams[k].max_diam, diams[kmd].max_diam)) {
      for (const auto& p : diams[k].far) far.push_back(p);
    } else if (diams[k].max_diam > diams[kmd].max_diam) {
      far = diams[k].far;
      kmd = k;
    }
  }
  assert(far.size() % 2 == 0);
}

void dynamic_diameter::diam_info::recompute_ABCD_and_P(const solution& s, int k,
                                                       int ignore) {
  assert(use_calipers);
  A = B = C = D = point(-1);
  x1 = x2 = y1 = y2 = -1;
  P.clear();

  for (int i : s.lot_assigned[k]) {
    if (i == ignore) continue;
    update_ABCD(i);
  }

  for (int i : s.lot_assigned[k]) {
    if (i == ignore) continue;
    point p(i);
    if (not in_rect(p)) P.push_back(p);
  }
  sort(begin(P), end(P));
}

void dynamic_diameter::diam_info::gain_point(const solution& s,
                                             const point& p) {
  if (use_calipers) {
    if (in_rect(p)) return;
    if (update_ABCD(p)) update_after_rect_change(s);
    add_to_sorted(P, p, false);
  } else {
    P.push_back(p);
  }
  comp_diam();
}

void dynamic_diameter::diam_info::lose_point(const solution& s,
                                             const point& p) {
  if (use_calipers) {
    if (in_rect(p) or P.empty()) return;
    remove_from_sorted(P, p);

    if (p == A or p == B or p == C or p == D) {
      recompute_ABCD_and_P(s, k, p.index);
    }
  } else {
    auto it = find(begin(P), end(P), p);
    assert(it != end(P));
    P.erase(it);
  }
  comp_diam();
}

bool dynamic_diameter::diam_info::update_ABCD(const point& p) {
  assert(use_calipers);
  bool updated = false;
  if (A.index == -1 or p.x - p.y > A.x - A.y) A = p, updated = true;
  if (B.index == -1 or p.x + p.y > B.x + B.y) B = p, updated = true;
  if (C.index == -1 or p.x - p.y < C.x - C.y) C = p, updated = true;
  if (D.index == -1 or p.x + p.y < D.x + D.y) D = p, updated = true;
  if (updated) {
    x1 = max(C.x, D.x);
    x2 = min(A.x, B.x);
    y1 = max(A.y, D.y);
    y2 = min(B.y, C.y);
    if (x1 > x2) swap(x1, x2);
    if (y1 > y2) swap(y1, y2);
    assert(P.size() < 4 or x1 <= x2);
    assert(P.size() < 4 or y1 <= y2);
    assert(not in_rect(A) and not in_rect(B) and not in_rect(C) and
           not in_rect(D));
    degenerate_rect = A == B or A == C or A == D or B == C or B == D or C == D;
    statistic("degenerate rect", degenerate_rect);
  }
  return updated;
}

void dynamic_diameter::diam_info::update_after_rect_change(const solution& s) {
  P.erase(remove_if(begin(P), end(P),
                    [&](const auto& p) { return this->in_rect(p); }),
          end(P));
  for (int i : s.lot_assigned[k]) {
    point p(i);
    if (in_rect(p)) continue;
    add_to_sorted(P, p, false);
  }
#ifdef DEBUG_DIAMETER
  assert(is_sorted(begin(P), end(P)));
#endif
}

void dynamic_diameter::diam_info::comp_diam_ch() {
  if (P.size() < 3) {
    comp_diam_bf();
    return;
  }

  compute_ch_and_calipers(P, U, L, max_diam, far);

#ifdef DEBUG_DIAMETER
  diam_info bf = *this;
  bf.comp_diam_bf();
  assert(f_eq(bf.max_diam, max_diam));
  auto &f1 = bf.far, &f2 = far;
  assert(f1.size() == f2.size());
  for (int i = 0; i < (int)f1.size(); i += 2) {
    bool ok = false;
    for (int j = 0; j < (int)f2.size(); j += 2)
      ok = ok or (f1[i] == f2[j] and f1[i + 1] == f2[j + 1]) or
           (f1[i] == f2[j + 1] and f1[i + 1] == f2[j]);
    assert(ok);
  }
#endif
}

void dynamic_diameter::diam_info::comp_diam_bf() {
  max_diam = -1;
  for (int i = 0; i < (int)P.size(); ++i) {
    const auto& p = P[i];
    for (int j = i + 1; j < (int)P.size(); ++j) {
      const auto& q = P[j];
      double d = dist2(p.index, q.index);
      if (f_eq(d, max_diam)) {
        far.push_back(p.index);
        far.push_back(q.index);
      } else if (d > max_diam) {
        max_diam = d;
        far = {p.index, q.index};
      }
    }
  }
  if (max_diam != -1) max_diam = sqrt(max_diam);
}