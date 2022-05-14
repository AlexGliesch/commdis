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
#include "local_improvement.h"
#include "balance_search.h"
#include "compactness_search.h"
#include "constructive.h"
#include "solution.h"
#include "weighted_search.h"

void tabu_search(solution& s, auto&& nbh) {
  solution inc = s; // The incumbent solution
  int num_moves = 0, last_improved = 0;
  tabu_list tabu;

  while (not timeout()) {
    int u = -1, k = -1;
    nbh.choose_next_move(inc, tabu, u, k);

    if (u != -1) {
      ++num_moves;
      if (not do_ls) tabu.mark_as_tabu(u);
      tabu.advance_iter();

      int p = inc.assigned[u];
      double bal_bef = inc.balance(), cmp_bef = inc.cmp(); 
      (void)bal_bef, (void)cmp_bef;
      inc.assign_node(u, k, true, true, true);
      nbh.update(inc, tabu, u, p, k);
      assert_correctness(inc, true);

      if (do_ls) 
        assert(f_leq(inc.balance(), bal_bef) or f_leq(inc.cmp(), cmp_bef));

      statistic("balanc impr. < EPS", f_eq(inc.balance(), bal_bef));
    } else {
      break;
    }

    if (nbh.better_sol(inc, s)) {
      s = inc;
      last_improved = num_moves;
    }

    if (num_moves - last_improved >= I_max) {
      break; 
    }

    statistic("num pairs", inc.diam.count());
  }
  statistic("search moves", num_moves);
  boost::hash_combine(execution_hash, s.balance() + s.diameter());
}

void local_improvement(solution& s) {
  pr("Start impr.: {}\n", s);

  if (do_weighted) {
    tabu_search(s, weighted_search(s));
    pr("Weighted: {}\n", s);
    return;
  }
  solution inc = s; 
  unordered_set<uint64_t> visited; 

  for (int i = 0; i < A_max and not timeout(); ++i) {
    if (do_cmp_search) {
      tabu_search(inc, compactness_search());
      pr("Optimize cmp.: {}\n", inc);
    }

    if (do_bal_search) {
      tabu_search(inc, balance_search(inc, (1.0 + alpha) * inc.diameter()));
      pr("Optimize bal.: {}\n", inc);
    }

    auto h = inc.hash();
    if (visited.count(h)) {
      break;
    }
    visited.insert(h);
    if (inc < s) s = inc;
  }

  if (s.balance() > 0 and do_bal_search)
    tabu_search(s, balance_search(s, d_max));

  pr("end impr.: {}\n", s);
}
