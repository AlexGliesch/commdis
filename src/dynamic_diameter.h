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
#include "geometry.h"
#include "global.h"
#define FAREPS (EPS * 2)

struct solution;

struct dynamic_diameter {
  struct diam_info {
    double max_diam = -1; 
    vector<point> P;      
    vector<point> U, L;   
    vector<int> far;      
                          
    int k;                

    point A, B, C, D;
    double x1, x2, y1, y2; 
    bool degenerate_rect;  

    void recompute_ABCD_and_P(const solution& s, int k, int ignore);

    bool update_ABCD(const point& p);

    void gain_point(const solution& s, const point& p);

    void lose_point(const solution& s, const point& p);

    void update_after_rect_change(const solution& s);

    void comp_diam() { use_calipers ? comp_diam_ch() : comp_diam_bf(); }
    void comp_diam_ch();
    void comp_diam_bf();

    bool in_rect(const point& p) const {
      if (P.size() < 4 or degenerate_rect) return false;
      assert(x1 <= x2 and y1 <= y2);
      return p.x > x1 and p.x < x2 and p.y > y1 and p.y < y2;
    }
  };

  void init(const solution& s);

  void recompute(const solution& s);

  pair<double, int> move_cost(const solution& s, int u, int k) const;

  pair<double, int> max_dist_u_in_k(const solution& s, int u, int k) const;

  void update(const solution& s, int u, int k, int p);

  double diameter() const { return lot_diameter(kmd); }

  double lot_diameter(int k) const { return diams[k].max_diam; }

  int count() const { return far.size() / 2; }

  void update_kmd_and_far();

  vector<diam_info> diams; 

  int kmd = -1; 

  vector<int> far; 
};

inline void consider_dist(double& FX, int& count, double dist, int c2) {
  if (f_eq(dist, FX))
    count += c2;
  else if (dist > FX)
    FX = dist, count = c2;
}