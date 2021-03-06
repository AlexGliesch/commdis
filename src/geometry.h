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
#include "global.h"

using namespace std;

struct point {
  point(int i = -1) : index(i) {
    if (i != -1) x = nodes[i].x, y = nodes[i].y;
  }
  bool operator==(const point& p) const { return p.index == index; }
  bool operator!=(const point& p) const { return p.index != index; }
  bool operator<(const point& p) const {
    return x < p.x or (x == p.x and y < p.y);
  }
  double x, y;
  int index;
};

inline double orientation(const point& p, const point& q, const point& r) {
  return (q.y - p.y) * (r.x - p.x) - (q.x - p.x) * (r.y - p.y);
}

inline ostream& operator<<(ostream& o, const point& p) {
  o << format("({:.2f}, {:.2f}, {})", p.x, p.y, p.index);
  return o;
}

void compute_ch_and_calipers(const vector<point>& P, vector<point>& U,
                             vector<point>& L, double& max_diameter,
                             vector<int>& far);