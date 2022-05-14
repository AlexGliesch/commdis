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
#include "global.h"
#include "dynamic_diameter.h"
#include "solution.h"

steady_clock::time_point start_time;
mt19937 rng;

int time_limit_seconds;
int max_iter;
uint64_t random_seed;
string instance_name;
string out_file;
int tabu_tenure;
double tenure_mult;
int I_max;
int A_max;
double alpha;
string seeds_alg;
string cons_alg;
bool do_cmp_search;
bool do_bal_search;
bool do_pi;
int pi_max;
bool do_only_cons;
bool do_ls;
bool do_weighted;
bool use_bridges;
bool use_calipers;
bool irace;
bool bootstrap;
size_t execution_hash = 0;

double tau[3] = {-1.0, -1.0, -1.0};
double d_max;
array<double, 3> mu;
array<double, 3> mu_tau;
int nnodes, ndistr, nedges;
vector<node> nodes;
vector<int> E;

map<string, pair<double, int>> statistics;
void statistic(const string& s, double value) {
  (void)s, (void)value;
#ifndef NDEBUG
  auto it = statistics.find(s);
  if (it == end(statistics))
    statistics[s] = make_pair(value, 1);
  else {
    it->second.first += value;
    it->second.second += 1;
  }
#endif
}

void print_statistics() {
#ifndef NDEBUG
  pr("- Statistics: \n");
  for (auto& p : statistics) {
    pr("{} - avg: {} val: {} total: {}\n", p.first,
       p.second.first / (double)p.second.second, p.second.first, p.second.second);
  }
#endif
}

void read_instance() {
  if (instance_name.find("/cygdrive/c/") != string::npos)
    instance_name = "C:/" + string(instance_name, 12);

  ifstream f(instance_name);
  if (not f or f.fail()) {
    cout << "Error: invalid instance name: " << instance_name << ". Exiting program.\n";
    exit(EXIT_FAILURE);
  }

  f >> nnodes;
  nodes.resize(nnodes);

  for3(j) mu[j] = 0;
  for (auto& n : nodes) {
    int index;
    f >> index >> n.x >> n.y >> n.a[0] >> n.a[1] >> n.a[2];
    for3(j) mu[j] += n.a[j];
  }
  f >> nedges;
  vector<vector<int>> nbs(nnodes);
  for (int i = 0; i < nedges; ++i) {
    int a, b;
    f >> a >> b;
    nbs[a].push_back(b);
    nbs[b].push_back(a);
  }
  for (int u = 0; u < nnodes; ++u) {
    nodes[u].nbs.i0 = E.size();
    for (auto v : nbs[u])
      E.push_back(v);
    nodes[u].nbs.i1 = E.size();
  }

  double dummy;
  if (ndistr == -1)
    f >> ndistr;
  else
    f >> dummy;
  f >> dummy;
  if (tau[0] == -1)
    f >> tau[0] >> tau[1] >> tau[2];
  else
    f >> dummy >> dummy >> dummy;

  ndistr = min(ndistr, nnodes);
  ndistr = max(ndistr, 2);

  for3(j) mu[j] /= (double)ndistr;

  if (use_calipers) {
    dynamic_diameter::diam_info dd;
    dd.A = dd.B = dd.C = dd.D = point(-1);
    for (int i = 0; i < nnodes; ++i)
      dd.update_ABCD(i);
    for (int i = 0; i < nnodes; ++i) {
      if (not dd.in_rect(i)) dd.P.emplace_back(i);
    }
    sort(begin(dd.P), end(dd.P));
    dd.comp_diam();
    d_max = dd.max_diam;
  } else {
    d_max = -1;
    for (int i = 0; i < nnodes; ++i)
      for (int j = i + 1; j < nnodes; ++j)
        d_max = max(d_max, dist(i, j));
  }

  for (auto& n : nodes) {
    n.value = 0;
    for3(j) n.value += n.a[j] / mu[j];
  }

  for3(j) mu_tau[j] = mu[j] * tau[j];

  tabu_tenure = ndistr * tenure_mult;

  pr("- Instance stats\n");
  pr("Name: {}\n", instance_name);
  pr("Nodes: {}\n", nnodes);
  pr("Edges: {}\n", nedges);
  pr("Lots: {}\n", ndistr);
  pr("tau: ({}, {}, {})\n", tau[0], tau[1], tau[2]);
  pr("mu: ({}, {}, {})\n", mu[0], mu[1], mu[2]);
  pr("d_max: {}\n", d_max);
  pr("tenure: {}\n", tabu_tenure);

#ifndef NDEBUG
  queue<int> q;
  vector<int> visited(nnodes, false);
  visited[0] = true;
  q.push(0);
  while (q.size()) {
    int v = q.front();
    q.pop();
    for (int u : nodes[v].nbs)
      if (not visited[u]) {
        visited[u] = true;
        q.push(u);
      }
  }
  bool connected = (count(begin(visited), end(visited), true) == nnodes);
  pr("Input graph is {}connected\n\n", (connected ? "" : "not "));
  assert(connected);
#endif

#ifndef NDEBUG
  set<pair<double, double>> m;
  for (const auto& n : nodes)
    m.emplace(n.x, n.y);
  assert(m.size() == nodes.size());
#endif
}

void read_cmd_line(int argc, char** argv) {
  namespace po = boost::program_options;
  po::options_description desc("grm17 - commercial districting");

  desc.add_options()("help", "Show help menu.");
  desc.add_options()("in", po::value<string>(&instance_name)->required(),
                     "Input filename.");
  desc.add_options()("time", po::value<int>(&time_limit_seconds)->default_value(2 << 28),
                     "Time limit (seconds). Default: no limit.");
  desc.add_options()("out", po::value<string>(&out_file)->required(), "Output filename.");
  desc.add_options()("seed", po::value<uint64_t>(&random_seed)->default_value(0),
                     "Random seed. If 0, a random random seed will be used.");
  desc.add_options()("iter", po::value<int>(&max_iter)->default_value(INT_MAX),
                     "Maximum number of multistart iterations. "
                     "Default: no limit.");
  desc.add_options()("tau", po::value<double>(&tau[0])->default_value(-1.0),
                     "Balance tolerance. If set, will override the tolerance "
                     "specified by the instance.");
  desc.add_options()("p", po::value<int>(&ndistr)->default_value(-1),
                     "Number of districts. If set, will override the value "
                     "specified by the instance.");
  desc.add_options()("cons-alg", po::value<string>(&cons_alg)->default_value("greedy-d"),
                     "Constructive algorithm. Possible values: {bfs, "
                     "greedy-d, greedy-b}.");
  desc.add_options()("seeds-alg", po::value<string>(&seeds_alg)->default_value("maxdisp"),
                     "Algorithm to find initial centers. Possible values: "
                     "{random, maxdisp}.");
  desc.add_options()("tenure", po::value<double>(&tenure_mult)->default_value(0.5),
                     "Tabu tenure (rel. to k). Default: k.");
  desc.add_options()("alpha", po::value<double>(&alpha)->default_value(0.00),
                     "Maximum relative increase to compactness allowed "
                     "during balance search.");
  desc.add_options()("imax", po::value<int>(&I_max)->default_value(500),
                     "Max. consecutive tabu search iterations without "
                     "improvement.");
  desc.add_options()("amax", po::value<int>(&A_max)->default_value(10),
                     "Maximum iterations of local improvement, i.e., "
                     "alternations between optimizing balance and compactness");
  desc.add_options()("only-cons", "Perform only construction.");
  desc.add_options()("ls", "Perform local search (instead of tabu).");
  desc.add_options()("no-cmp", "Disable compactness search.");
  desc.add_options()("no-bal", "Disable balance search.");
  desc.add_options()("no-pi", "Disable compactness path-improvement.");
  desc.add_options()("pimax", po::value<int>(&pi_max)->default_value(25),
                     "Maximum attemps at path-improvement.");
  desc.add_options()("naive-diameter", "Compute diameters by brute force (O(n^2)).");
  desc.add_options()("naive-connectivity",
                     "Do connectivity tests by brute force (with a BFS).");
  desc.add_options()("weighted",
                     "Search with a weighted obj. fun., instead of alternating "
                     "bal. and cmp.");
  desc.add_options()("bootstrap", "Do bootstrapping (i.e., output results after every "
                                  "iteration).");
  desc.add_options()("irace", "Result for irace (i.e., output is only one value).");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      cout << desc << endl;
      exit(EXIT_SUCCESS);
    } else {
      po::notify(vm);
    }

    do_cmp_search = not vm.count("no-cmp");
    do_bal_search = not vm.count("no-bal");
    do_pi = not vm.count("no-cmp-pi");
    use_bridges = not vm.count("naive-conn");
    use_calipers = not vm.count("naive-diam");
    do_only_cons = vm.count("only-cons");
    do_ls = vm.count("ls");
    do_weighted = vm.count("weighted");
    bootstrap = vm.count("bootstrap");
    irace = vm.count("irace");

    if (bootstrap and irace) {
      pr("Warning: bootstrap and irace both enabled. Defaulting to just "
         "bootstrap.");
      irace = false;
    }

    if (do_cmp_search + do_bal_search == 0) do_only_cons = true;

    tau[1] = tau[2] = tau[0];

    if (random_seed == 0) random_seed = choose_unique_random_seed() % UINT32_MAX;
    execution_hash = random_seed;
    rng.seed(random_seed);

    pr("- Seed\n{}\n\n", random_seed);

    if (seeds_alg != "maxdisp" and seeds_alg != "random") {
      cout << "Error: invalid initial nodes algorithm." << endl;
      exit(EXIT_SUCCESS);
    }
    if (cons_alg != "bfs" and cons_alg != "greedy-b" and cons_alg != "greedy-d") {
      cout << "Error: invalid constructive algorithm." << endl;
      exit(EXIT_SUCCESS);
    }

  } catch (po::error& e) {
    cout << "Error parsing command line: " << e.what() << endl;
    cout << desc << endl;
    exit(EXIT_FAILURE);
  }
}
