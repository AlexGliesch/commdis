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
#include "geometry.h"
#include "global.h"

void compute_ch_and_calipers(const vector<point>& P, vector<point>& U,
                             vector<point>& L, double& max_diameter,
                             vector<int>& far) {
#ifdef DEBUG_GEOMETRY
  assert(is_sorted(begin(P), end(P)));
#endif

  assert(P.size() >= 3);
  U.resize(P.size()), L.resize(P.size());
  int usz = 0, lsz = 0;
  for (const auto& p : P) {
    while (usz > 1 and orientation(U[usz - 2], U[usz - 1], p) <= 0) --usz;
    while (lsz > 1 and orientation(L[lsz - 2], L[lsz - 1], p) >= 0) --lsz;
    U[usz++] = p;
    L[lsz++] = p;
  }
  assert(usz + lsz >= 3);
  U.resize(usz), L.resize(lsz);

  max_diameter = -1;

  int i = 0, j = lsz - 1;
  while (i < usz - 1 or j > 0) {
    double dist = dist2(U[i].index, L[j].index);
    if (f_eq(dist, max_diameter)) {
      far.push_back(U[i].index), far.push_back(L[j].index);
    } else if (dist > max_diameter) {
      max_diameter = dist;
      far = {U[i].index, L[j].index};
    }

    if (i == usz - 1)
      --j;
    else if (j == 0)
      ++i;
    else if ((U[i + 1].y - U[i].y) * (L[j].x - L[j - 1].x) >
             (L[j].y - L[j - 1].y) * (U[i + 1].x - U[i].x))
      ++i;
    else
      --j;
  }

  while (U.size() and (U.front() == L.front() or U.front() == L.back())) {
    swap(U.front(), U.back());
    U.pop_back();
  }
  while (U.size() and (U.back() == L.front() or U.back() == L.back()))
    U.pop_back();

  assert(far.size() >= 2);
  if (max_diameter != -1) max_diameter = sqrt(max_diameter);
}