#include "pre.h"
#include "delaunay.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>

using namespace std;

// Global constants
// x, y will be within [min_coord, max_coord]
double minXY = 0, maxXY = 1000;

// Global variables
vector<vector<int>> g;     // adjacency graph
vector<XYZ> P;             // vertex locations
vector<double> a1, a2, a3; // vertex activity values
int num_edges;             // current number of edges
mt19937 rng;

// Command line parameters
int n, m, num_districts;
double tau;
string output_filename;
string weights;
using seed_type = uint16_t;
seed_type seed;
bool make_grid;
bool make_tree;
string fusy;
string location;

void read_cmd_line(int argc, char** argv) {
  namespace po = boost::program_options;
  po::options_description desc("Districting instance generator.");
  desc.add_options()("help", "");
  desc.add_options()("n", po::value<int>(&n)->required(), "Number of vertices");
  desc.add_options()(
      "m", po::value<int>(&m)->default_value(numeric_limits<int>::max()),
      "Number of edges");
  desc.add_options()("p", po::value<int>(&num_districts)->required(),
                     "Number of districts.");
  desc.add_options()("out",
                     po::value<string>(&output_filename)->default_value(""),
                     "Output filename; if unset, a default name will be used.");
  desc.add_options()("tau", po::value<double>(&tau)->default_value(0.05),
                     "Default threshold for the three activities.");
  desc.add_options()("seed", po::value<seed_type>(&seed)->default_value(0),
                     "Random seed. if 0, a random value will be used.");
  desc.add_options()("grid", "Generate a grid graph.");
  // TODO change parameter name and desc.
  desc.add_options()("fusy", po::value<string>(&fusy)->default_value(""),
                     "Read a \"fusy\" graph from input.");
  desc.add_options()("tree", "Generate a tree graph.");
  desc.add_options()("min-xy", po::value<double>(&minXY),
                     "Minimum (x,y) coordinate.");
  desc.add_options()("max-xy", po::value<double>(&maxXY),
                     "Maximum (x,y) coordinate.");
  desc.add_options()(
      "location", po::value<string>(&location)->default_value("uniform"),
      "Node location strategy, in: [uniform, jitter-grid, poisson].");
  desc.add_options()("weights",
                     po::value<string>(&weights)->default_value("DT"),
                     "How unit activities are determined: \n"
                     "DS: activities are selected in the ranges [4,20], "
                     "[15,400] and [15,100], as in [RMF09].\n"
                     "DT: activities are determined by the sum of 68 random "
                     "variables, as in [RMF09].\n"
                     "RR: same as DS, but the range limits [v_min, v_max] for"
                     " each activity are determined randomly: v_min \\in "
                     "[1,100] and v_max \\in [v_min, 5000].");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      cout << desc << endl;
      exit(EXIT_SUCCESS);
    } else {
      po::notify(vm);
    }

    make_grid = vm.count("grid");
    make_tree = vm.count("tree");

    if (m < (n - 1)) {
      cout << "Number of edges (m) must be at least n-1." << endl;
      exit(EXIT_FAILURE);
    }

    if (seed == 0) seed = choose_unique_random_seed();
    rng.seed(seed);
    cout << "- Seed " << seed << endl;

    if (weights != "RR" and weights != "DS" and weights != "DT") {
      cout << "Invalid 'weights' parameter. See --help." << endl;
      exit(EXIT_FAILURE);
    }

    if (location != "uniform" and location != "jitter-grid" and
        location != "poisson") {
      cout << "Invalid 'location' parameter. See --help." << endl;
      exit(EXIT_FAILURE);
    }
  } catch (po::error& e) {
    cout << "Error parsing command line: " << e.what() << endl;
    cout << desc << endl;
    exit(EXIT_FAILURE);
  }
}

void assert_connected() {
  queue<int> q;
  vector<bool> visited(n, false);
  q.push(0);
  visited[0] = true;
  while (q.size()) {
    int v = q.front();
    q.pop();
    for (int n : g[v])
      if (not visited[n]) {
        visited[n] = true;
        q.push(n);
      }
  }
  assert(all_of(begin(visited), end(visited), [](bool b) { return b; }));
}

bool same_point(const XYZ& a, const XYZ& b) {
  return abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6;
}

int find_index(const XYZ& a) {
  int i = -1;
  for (; i < n; ++i)
    if (same_point(a, P[i])) break;
  return i;
}

void add_edge(int p1, int p2) {
  if (find(begin(g[p1]), end(g[p1]), p2) == end(g[p1])) g[p1].push_back(p2);
  if (find(begin(g[p2]), end(g[p2]), p1) == end(g[p2])) g[p2].push_back(p1);
}

set<pair<int, int>> edge_set() {
  set<pair<int, int>> edges;
  for (int i = 0; i < n; ++i)
    for (int n : g[i]) edges.emplace(min(i, n), max(i, n));
  return (edges);
}

vector<pair<int, int>> compute_edges_without_ST() {
  pr("Edge set: \n");
  auto edges = edge_set();
  num_edges = edges.size();
  pr("num_edges: {}\n", num_edges);

  pr("BFS\n");
  queue<int> q;
  vector<bool> visited(n, false);
  q.push(0);
  visited[0] = true;
  while (q.size()) {
    int v = q.front();
    q.pop();
    for (int n : g[v])
      if (not visited[n]) {
        visited[n] = true;
        q.push(n);
        edges.erase(make_pair(min(v, n), max(v, n)));
      }
  }
  assert((int)edges.size() == num_edges - (n - 1));

  pr("BFS done\n");

  vector<pair<int, int>> ewst;
  for (auto& e : edges) ewst.push_back(e);
  return ewst;
}

void DS_value(double& r1, double& r2, double& r3) {
  r1 = rand_double(4, 20);
  r2 = rand_double(15, 400);
  r3 = rand_double(15, 100);
}

int random_variable(const map<int, double>& m) {
  double r = rand_double(0, 1), c = 0;
  for (auto& p : m) {
    c += p.second;
    if (r < c) return p.first;
  }
  return 0;
}

void DT_value(double& r1, double& r2, double& r3) {
  r1 = r2 = r3 = 0;
  map<int, double> m1 = {{0, 0.15}, {3, 0.15}, {1, 0.35}, {1, 0.35}};
  map<int, double> m2 = {{1, 0.01}, {12, 0.01}, {2, 0.03}, {11, 0.03},
                         {3, 0.06}, {10, 0.06}, {4, 0.1},  {9, 0.1},
                         {5, 0.21}, {8, 0.12},  {6, 0.18}, {7, 0.18}};

  for (int i = 0; i < 68; ++i) {
    int d1 = random_variable(m1);
    if (d1 > 0) {
      int d2 = random_variable(m2), d3 = random_variable(m2);
      r1 += d1, r2 += d2, r3 += d3;
    }
  }
}

void RR_value(double& r1, double& r2, double& r3) {
  bool chosen_ranges = false;
  double min1, max1, min2, max2, min3, max3;
  if (not chosen_ranges) {
    double minmin = 5, minmax = 100;
    double maxmax = 5000;
    chosen_ranges = true;
    min1 = rand_double(minmin, minmax);
    max1 = rand_double(min1, maxmax);
    min2 = rand_double(minmin, minmax);
    max2 = rand_double(min2, maxmax);
    min3 = rand_double(minmin, minmax);
    max3 = rand_double(min3, maxmax);
  }
  r1 = rand_double(min1, max1);
  r2 = rand_double(min2, max2);
  r3 = rand_double(min3, max3);
}

void write_output_file() {
  auto edges = edge_set();

  if (output_filename.empty()) {
    stringstream ss;
    if (make_grid)
      ss << "grid-";
    else if (fusy.size())
      ss << "fusy-";
    else if (make_tree)
      ss << "tree-";
    else
      ss << "del-";
    ss << "n" << n << "-k" << num_districts << "-s" << seed << ".in";
    output_filename = ss.str();
  }

  ofstream f(output_filename);
  pr("Write nodes\n");
  f << n << endl;
  for (int i = 0; i < n; ++i) {
    f << i << " " << P[i].x << " " << P[i].y << " " << a1[i] << " " << a2[i]
      << " " << a3[i] << endl;
  }
  pr("Write edges\n");
  f << edges.size() << endl;
  for (auto& e : edges) f << e.first << " " << e.second << endl;
  f << num_districts << " " << 50 << " " << tau << " " << tau << " " << tau << endl;
}

void generate_points() {
  pr("- Generating random points\n");
  P.resize(n + 3);
  if (location == "uniform") {
    for (auto& p : P) {
      p.x = rand_double(minXY, maxXY);
      p.y = rand_double(minXY, maxXY);
    }
  } else if (location == "poisson") {
    double dist = 3 * (maxXY - minXY) / double(P.size());
    for (int i = 0; i < (int)P.size(); ++i) {
      double minD;
      do {
        P[i].x = rand_double(minXY, maxXY);
        P[i].y = rand_double(minXY, maxXY);
        minD = numeric_limits<double>::max();
        for (int j = 0; j < i; ++j)
          minD = min(minD, pow(P[i].x - P[j].x, 2) + pow(P[i].y - P[j].y, 2));
      } while (sqrt(minD) <= dist);
    }
  } else if (location == "jitter-grid") {
    double side = (maxXY - minXY) / sqrt(P.size());
    for (int i = 0; i < (int)P.size(); ++i) {
      int yi = i / (int)sqrt(P.size());
      int xi = i % (int)sqrt(P.size());
      P[i].x = rand_double(xi * side, (xi + 1) * side);
      P[i].y = rand_double(yi * side, (yi + 1) * side);
    }
  }
}

void graph_delaunay() {
  // Generate random points
  generate_points();

  // Generate triangulation
  pr("- Triangulating\n");
  vector<ITRIANGLE> tri(n * 3);
  int ntri = 0;
  qsort(&P[0], n, sizeof(XYZ), XYZCompare);
  cout << "Triangulate" << endl << flush;
  Triangulate(n, &P[0], &tri[0], ntri);
  cout << "ntri: " << ntri << endl;

  // Populate the graph g with the edges from the triangulation
  pr("- Populating graph\n");
  g.assign(n, {});
  for (int i = 0; i < ntri; ++i) {
    auto& t = tri[i];
    assert(t.p1 != t.p2 and t.p2 != t.p3 and t.p1 != t.p3);
    assert(t.p1 >= 0 and t.p1 < n);
    assert(t.p2 >= 0 and t.p2 < n);
    assert(t.p3 >= 0 and t.p3 < n);
    add_edge(t.p1, t.p2);
    add_edge(t.p2, t.p3);
    add_edge(t.p1, t.p3);
  }
  P.resize(n);
  pr("- Triangulation done\n");
}

void graph_grid(int neighborhood_size = 4) {
  int l = sqrt(n);
  n = l * l;
  g.assign(l * l, {});
  P.resize(l * l);

  for (int i = 0; i < l; ++i)
    for (int j = 0; j < l; ++j) {
      P[i * l + j].x = i;
      P[i * l + j].y = j;

      static int di[] = {0, 0, 1, -1, 1, 1, -1, -1},
                 dj[] = {-1, 1, 0, 0, 1, -1, 1, -1};
      for (int k = 0; k < neighborhood_size; ++k) {
        int ii = i + di[k], jj = j + dj[k];
        if (ii >= 0 and ii < l and jj >= 0 and jj < l)
          add_edge(i * l + j, ii * l + jj);
      }
    }
}

void boost_get_embedding(int n, const set<pair<int, int>>& edge_list) {
  using namespace boost;

  struct coord_t {
    size_t x, y;
  };

  using graph =
      adjacency_list<vecS, vecS, undirectedS, property<vertex_index_t, int>,
                     property<edge_index_t, int>>;

  graph g(n);
  for (const auto& p : edge_list) boost::add_edge(p.first, p.second, g);

  std::vector<int> component(num_vertices(g));
  int cc = connected_components(g, &component[0]);
  if (cc != 1) {
    pr("Error: graph is not connected.\n");
    exit(EXIT_FAILURE);
  } else
    pr("Graph is connected.\n");

  // Initialize the interior edge index
  property_map<graph, edge_index_t>::type e_index = get(edge_index, g);
  graph_traits<graph>::edges_size_type edge_count = 0;
  graph_traits<graph>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    put(e_index, *ei, edge_count++);

  // Test for planarity; compute the planar embedding as a side-effect
  using vec_t = vector<graph_traits<graph>::edge_descriptor>;
  vector<vec_t> embedding(num_vertices(g));

  bool is_planar = boyer_myrvold_planarity_test(
      boyer_myrvold_params::graph = g,
      boyer_myrvold_params::embedding = &embedding[0]);

  if (not is_planar) {
    pr("Error: graph is not planar.\n");
    exit(EXIT_FAILURE);
  } else
    pr("Graph is planar.\n");

  pr("make_biconnected_planar, num_edges: {}\n", boost::num_edges(g));
  make_biconnected_planar(g, &embedding[0]);

  // Re-initialize the edge index, since we just added a few edges
  edge_count = 0;
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    put(e_index, *ei, edge_count++);

  // Test for planarity again; compute the planar embedding as a side-effect
  is_planar = boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                           boyer_myrvold_params::embedding =
                                               &embedding[0]);

  pr("make_maximal_planar, num_edges: {}\n", boost::num_edges(g));
  make_maximal_planar(g, &embedding[0]);
  pr("max num_edges: {}\n", boost::num_edges(g));

  // Re-initialize the edge index, since we just added a few edges
  edge_count = 0;
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    put(e_index, *ei, edge_count++);

  is_planar = boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                           boyer_myrvold_params::embedding =
                                               &embedding[0]);

  // Find a canonical ordering
  vector<graph_traits<graph>::vertex_descriptor> ordering;
  planar_canonical_ordering(g, &embedding[0], back_inserter(ordering));

  // Set up a property map to hold the mapping from vertices to coord_t's
  using straight_line_drawing_storage_t = vector<coord_t>;
  using straight_line_drawing_t =
      boost::iterator_property_map<straight_line_drawing_storage_t::iterator,
                                   property_map<graph, vertex_index_t>::type>;

  straight_line_drawing_storage_t straight_line_drawing_storage(
      num_vertices(g));
  straight_line_drawing_t straight_line_drawing(
      straight_line_drawing_storage.begin(), get(vertex_index, g));

  // Compute the straight line drawing
  chrobak_payne_straight_line_drawing(g, embedding, ordering.begin(),
                                      ordering.end(), straight_line_drawing);

  graph_traits<graph>::vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
    coord_t coord(get(straight_line_drawing, *vi));
    P[*vi].x = coord.x;
    P[*vi].y = coord.y;
  }

  // Verify that the drawing is actually a plane drawing
  if (not is_straight_line_drawing(g, straight_line_drawing)) {
    pr("Is NOT a plane drawing.\n");
    exit(EXIT_FAILURE);
  } else
    pr("Is a plane drawing.\n");
}

void random_walk(int sz) {
  // Assume g is connected
  pr("random_walk({}), g: {}\n", sz, g.size());
  assert(sz <= (int)g.size());
  vector<bool> visited(g.size(), false);
  int start = rand_int(0, g.size() - 1);
  int v = start;
  int num_visited = 0, num_moves = 0;
  while (num_visited < sz) {
    ++num_moves;
    if (not visited[v]) visited[v] = true, ++num_visited;
    //     if (rand_double(0.0, 1.0) < 0.15)
    //       v = start;
    //     else
    v = g[v][rand_int(0, g[v].size() - 1)];
  }
  assert(count(begin(visited), end(visited), true) == sz);
  pr("random_walk, moves: {}, visited: {}\n", num_moves, num_visited);

  unordered_map<int, int> index_map;
  int x = 0;
  for (int i = 0; i < (int)g.size(); ++i)
    if (visited[i]) index_map[i] = x++;

  assert(x == sz);
  vector<vector<int>> g2(sz);
  for (int i = 0; i < (int)g.size(); ++i)
    if (visited[i])
      for (int j : g[i])
        if (visited[j]) g2[index_map[i]].push_back(index_map[j]);
  g = move(g2);
}

void graph_fusy(const string& filename) {
  ifstream f(filename);
  string dummy;
  f >> dummy >> dummy;

  int i, j;
  unordered_map<int, int> vertex_indices;
  int v_count = 0;
  g.clear();
  while (f >> i >> j) {
    if (i == 0 and j == 0) break;
    for (int k : {i, j})
      if (vertex_indices.count(k) == 0) {
        vertex_indices[k] = v_count++;
        g.emplace_back();
      }
    int vi = vertex_indices[i], vj = vertex_indices[j];
    assert((int)g.size() > max(vi, vj));
    add_edge(vi, vj);
  }
  pr("Read {} vertices.\n", g.size());

  if ((int)g.size() <= n) {
    pr("Warning: specified n is larger than graph size, bounding n={}\n",
       g.size());
    n = g.size();
  } else {
    random_walk(n);
    assert(n == (int)g.size());
  }

  P.resize(n);
  set<pair<int, int>> edge_list;
  for (int i = 0; i < n; ++i)
    for (int j : g[i]) edge_list.emplace(min(i, j), max(i, j));

  boost_get_embedding(n, edge_list);

  double maxy = -1, maxx = -1, miny = DBL_MAX, minx = DBL_MAX;
  for (auto& p : P) {
    maxx = max(p.x, maxx);
    minx = min(p.x, minx);
    maxy = max(p.y, maxy);
    miny = min(p.y, miny);
  }

  double num = 16;
  double stepx = (maxx + minx) / num;
  double stepy = (maxy + miny) / num;
  vector<vector<int>> q(num, vector<int>(num, 0));
  for (auto& p : P) {
    for (int i = 0; i < num; ++i)
      for (int j = 0; j < num; ++j) {
        double x = minx + i * stepx;
        double y = miny + j * stepy;
        q[i][j] +=
            (p.x >= x and p.x < x + stepx and p.y >= y and p.y < y + stepy);
      }
  }

  pr("q: \n");
  for (int i = 0; i < num; ++i) {
    for (int j = 0; j < num; ++j)
      cout << setw(4) << setfill(' ') << q[j][i] << " ";
    cout << endl;
  }
}

int main(int argc, char** argv) {
  // Read command line arguments
  read_cmd_line(argc, argv);

  pr("- Parameters\n");
  pr("n: {}\n", n);
  pr("m: {}\n", m);
  pr("districts: {}\n", num_districts);
  pr("weights: {}\n", weights);
  pr("output file: {}\n", output_filename);

  if (fusy.size()) {
    graph_fusy(fusy);
  } else if (make_grid) {
    graph_grid();
  } else {
    if (make_tree) m = n - 1;
    graph_delaunay();
  }

  // Find some spanning tree of the triangulation (this is done by a BFS)
  pr("- Finding spanning tree\n");
  auto edges = compute_edges_without_ST();

  // Remove edges not in the spanning tree until the desired number of edges 'm'
  // is met
  pr("- Removing edges\n");
  while (num_edges > m and edges.size() > 0) {
    --num_edges;
    int e = rand_int(0, edges.size() - 1);
    int p1 = edges[e].first, p2 = edges[e].second;
    g[p1].erase(find(begin(g[p1]), end(g[p1]), p2));
    g[p2].erase(find(begin(g[p2]), end(g[p2]), p1));
    edges.erase(begin(edges) + e);
  }

  // Assert that the resulting graph is connected
  assert_connected();

  // Determine activity values
  pr("- Determining activity values\n");
  a1.resize(n), a2.resize(n), a3.resize(n);
  for (int i = 0; i < n; ++i) {
    if (weights == "DS")
      DS_value(a1[i], a2[i], a3[i]);
    else if (weights == "DT")
      DT_value(a1[i], a2[i], a3[i]);
    else if (weights == "RR")
      RR_value(a1[i], a2[i], a3[i]);
  }

  // Write output file
  pr("- Writing output\n");
  write_output_file();
}
