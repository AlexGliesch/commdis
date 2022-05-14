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
#include "solution.h"

struct tabu_list;

struct weighted_search {
  weighted_search(const solution& s);

  void choose_next_move(const solution& s, const tabu_list& tabu, int& u,
                        int& k) const;

  static bool better_sol(const solution& a, const solution& b) {
    return f_less(value(a), value(b));
  }

  void update(const solution& s, const tabu_list&, int u, int p, int k);

  void recompute_moves(solution&, const tabu_list&) {} 

private:
  static double alphamix(double bal, double cmp) {
    static const double Alpha = 0.3;
    return Alpha * bal + (1.0 - Alpha) * (cmp / d_max);
  }

  static double value(const solution& s) {
    return alphamix(s.balance(), s.cmp());
  }

  static double move_value(const solution& s, int u, int k) {
    return alphamix(s.balance_cost(u, k), s.cmp_cost(u, k));
  }

  void compute_move_of_lot(const solution& s, const tabu_list& tabu, int k);

  vector<int> cnd;

  queue<int> tabu_q;
};
