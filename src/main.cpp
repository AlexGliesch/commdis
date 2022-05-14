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
#include "constructive.h"
#include "global.h"
#include "local_improvement.h"
#include "solution.h"

int repl = 0;
vector<double> times;
int num_feasible;
vector<double> bal, cmp;
vector<double> cbal, ccmp;

void print_bootstrap(double cbal, double ccmp, double bal, double cmp, int num) {
  fmt::print("{} {} {} {:.4f} ", boost::filesystem::path(instance_name).stem().string(),
             nnodes, ndistr, times[num] - (num > 0 ? times[num - 1] : 0));
  fmt::print("{:.4f} {:.4f} {:.4f} {:.4f} ", cbal, ccmp, bal, cmp);
  fmt::print("{} {}", num + 1, random_seed);
  fmt::print("\n");
}

int main(int argc, char** argv) {
  read_cmd_line(argc, argv);
  read_instance();
  start_time = steady_clock::now();
  int best = -1;
  double ttb = 0;
  solution best_sol;

  while (repl++ < max_iter and not timeout()) {
    solution s = constructive();

    s.recompute_all(true, true);
    assert_correctness(s, true);

    double xcbal = s.balance(), xccmp = s.cmp();

    if (not timeout() and not do_only_cons) {
      local_improvement(s);
    }

    assert_correctness(s, not do_only_cons);

    pr("#{}: {}\n", repl, s);

    if (not s.is_incomplete()) {
      num_feasible += f_eq(s.balance(), 0);
      bal.push_back(s.balance()), cmp.push_back(s.cmp());
      cbal.push_back(xcbal), ccmp.push_back(xccmp);
      times.push_back(elapsed_time(start_time));
      if (bootstrap) {
        print_bootstrap(cbal.back(), ccmp.back(), bal.back(), cmp.back(),
                        times.size() - 1);
      }

      if (best == -1 or lexi_less(bal.back(), cmp.back(), bal[best], cmp[best]))
        best = bal.size() - 1, ttb = times[best], best_sol = s;
    }
  }

  if (bootstrap) {
    exit(EXIT_SUCCESS);
  }

  int N = times.size();

  if (irace) {
    if (N > 0) {
      double irace_val;
      if (f_greater(bal[best], 0)) {
        irace_val = 1.0 + bal[best];
      } else {
        irace_val = cmp[best] / d_max;
      }
      fmt::print("{}", irace_val);
    } else {
      fmt::print("{}", double(INT_MAX));
    }
    exit(EXIT_SUCCESS);
  }

  if (not out_file.empty()) {
    ofstream f(out_file);
    if (not f.good()) {
      pr("Error: could not open output file.\n");
      exit(EXIT_FAILURE);
    }
    for (int i = 0; i < (int)best_sol.assigned.size(); ++i) {
      if (i) f << ' ';
      f << best_sol.assigned[i];
    }
    f.close();
  }

  double avg_bal = 0, avg_cmp = 0, avg_cbal = 0, avg_ccmp = 0;
  double sdev_bal = 0, sdev_cmp = 0, sdev_cbal = 0, sdev_ccmp = 0;
  calc_avg_and_sdev(begin(bal), end(bal), avg_bal, sdev_bal);
  calc_avg_and_sdev(begin(cmp), end(cmp), avg_cmp, sdev_cmp);
  calc_avg_and_sdev(begin(cbal), end(cbal), avg_cbal, sdev_cbal);
  calc_avg_and_sdev(begin(ccmp), end(ccmp), avg_ccmp, sdev_ccmp);

  pr("\n- Finished, results:\n");
  pr("Instance: {}\n", boost::filesystem::path(instance_name).stem());
  pr("Replications: {}\n", repl - 1);
  pr("Num. feasible: {}\n", num_feasible);
  pr("T.t.b.: {:.3f}\n", ttb);
  pr("Time: {:.3f}\n", elapsed_time(start_time));
  pr("Balance: {:.3f} (avg.: {:.3f}, std. dev.: {:.3f})\n", bal[best], avg_bal, sdev_bal);
  pr("Compactness: {:.3f} (avg.: {:.3f}, std. dev.: {:.3f})\n", cmp[best], avg_cmp,
     sdev_cmp);
  pr("Balance (cons.): {:.3f} (avg.: {:.3f}, std. dev.: {:.3f})\n", cbal[best], avg_cbal,
     sdev_cbal);
  pr("Compactness (cons.): {:.3f} (avg.: {:.3f}, std. dev.: {:.3f})\n", ccmp[best],
     avg_ccmp, sdev_ccmp);
  pr("Seed: {}\n", random_seed);
  pr("Execution hash: {}\n", execution_hash % size_t(pow(10, 6)));
  pr("\n");
  print_statistics();

  fmt::print("{} {} {} {:.4f} {:.4f} ",
             boost::filesystem::path(instance_name).stem().string(), nnodes, ndistr,
             elapsed_time(start_time), ttb);
  //
  fmt::print("{:.4f} {:.4f} {:.4f} {:.4f} ", cbal[best], ccmp[best], bal[best],
             cmp[best]);
  //
  fmt::print("{} {} {}", repl - 1, num_feasible, random_seed);
  //
  fmt::print("\n");
}
