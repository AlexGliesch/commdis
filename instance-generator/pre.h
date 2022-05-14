#pragma once
#include <fmt/printf.h>
#include <algorithm>
#include <array>
#include <boost/filesystem.hpp>
#include <boost/functional/hash.hpp>
#include <boost/heap/binomial_heap.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/pairing_heap.hpp>
#include <boost/program_options.hpp>
#include <cassert>
#include <cfloat>
#include <chrono>
#include <ciso646>
#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <queue>
#include <random>
#include <set>
#include <stack>
#include <string>
#include <tuple>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#define mp make_pair
#define tstr to_string

// near and far are defined by minwindef.h (for some reason) and breaks my code
#undef near
#undef far

using namespace std;

/// Float handling
#define EPS 1e-8
inline bool f_eq(double a, double b) { return abs(a - b) < EPS; }
inline bool f_leq(double a, double b) { return a < b or f_eq(a, b); }
inline bool f_less(double a, double b) { return a < b and not f_eq(a, b); }
inline bool f_geq(double a, double b) { return f_leq(b, a); }
inline bool f_greater(double a, double b) { return f_less(b, a); }

/// Time-tracking
using namespace std::chrono;
extern steady_clock::time_point start_time;
// Elapsed time since time_point
inline double elapsed_time(const steady_clock::time_point& start) {
  return duration_cast<duration<double, milli>>(steady_clock::now() - start)
             .count() /
         1000.0;
}

/// Random number generation
extern mt19937 rng;
// Random inter in range [from, to] inclusive
inline int rand_int(int from, int to) {
  static uniform_int_distribution<int> d;
  return d(rng, decltype(d)::param_type{from, to});
}
// See above, for double, also inclusive
inline double rand_double(double from, double to) {
  static std::uniform_real_distribution<double> d;
  return d(rng, decltype(d)::param_type{from, to});
}
// Choose a unique random seed based on system clock(), time() and getpid()
inline uint64_t choose_unique_random_seed() {
  uint64_t a = (uint64_t)clock(), b = (uint64_t)time(nullptr),
           c = (uint64_t)getpid();
  a = (a - b - c) ^ (c >> 13);
  b = (b - c - a) ^ (a << 8);
  c = (c - a - b) ^ (b >> 13);
  a = (a - b - c) ^ (c >> 12);
  b = (b - c - a) ^ (a << 16);
  c = (c - a - b) ^ (b >> 5);
  a = (a - b - c) ^ (c >> 3);
  b = (b - c - a) ^ (a << 10);
  c = (c - a - b) ^ (b >> 15);
  return c;
}

// Select k elements from v and put them in 'result', using reservoir sampling
template <typename T>
inline void choice(const vector<T>& v, vector<T>& result, int k) {
  int n = v.size();
  result.resize(k);
  for (int i = 0; i < k; ++i) result[i] = v[i];
  for (int i = k + 1; i < n; ++i) {
    int j = rand_int(0, i);
    if (j < k) result[j] = v[i];
  }
}

// Class to handle reservoir sampling more easily
// Could we be using ints here instead of double?
struct reservoir_sampling {
  bool consider() { return (1.0 / ++num) >= rand_double(0.0, 1.0); }
  void reset() { num = 0.0; }
  double num = 0.0;
};

/// Message logging, debugging
inline void empty_func(const auto&, ...) {}
#ifdef NDEBUG
#define pr empty_func
#else
#define pr fmt::print
#endif
using fmt::format;

// Print a range "v" to a string
inline string fmt_range(const auto& v, const string& sep = " ") {
  string s;
  for (auto it = begin(v); it != end(v); ++it) {
    s += format("{}", *it);
    if (next(it) != end(v)) s += sep;
  }
  return (s);
}

struct debug {
  template <typename... Args> debug(const char* msg, const Args&... args) {
    s = fmt::format(msg, args...);
    pr("Call {}\n", s);
  }
  ~debug() { pr("Returning {}\n", s); }
  string s;
};

/// Utility STL
template <typename T> int num_duplicates(const vector<T>& v) {
  auto r2 = v;
  sort(begin(r2), end(r2));
  auto it = unique(begin(r2), end(r2));
  r2.resize(it - begin(r2));
  return v.size() - r2.size();
}

template <typename T> bool is_unique(const vector<T>& v) {
  return num_duplicates(v) == 0;
}

// Remove value t from a sorted vector. If t appears more than once in v, only
// one instance is removed.
template <typename T> void remove_from_sorted(vector<T>& v, const T& t) {
  auto it = lower_bound(begin(v), end(v), t);
  if (*it == t) v.erase(it);
}

// Insert value t in a sorted vector v. If allow_duplicates=false, a value will
// not be inserted if it already appears in the range.
template <typename T>
void add_to_sorted(vector<T>& v, const T& t, bool allow_duplicates) {
  auto ub = upper_bound(begin(v), end(v), t);
  if (allow_duplicates or ub == begin(v) or *(ub - 1) != t) v.insert(ub, t);
}

template <typename P1, typename P2>
ostream& operator<<(ostream& o, const pair<P1, P2>& p) {
  o << "(" << p.first << ", " << p.second << ")";
  return o;
}

void calc_avg_and_sdev(auto first, auto last, double& avg, double& sdev) {
  size_t sz = last - first;
  if (sz == 0) {
    avg = sdev = -1;
    return;
  }
  avg = accumulate(first, last, 0.0) / double(sz);
  sdev = 0;
  for (auto x = first; x != last; ++x) sdev += pow(*x - avg, 2);
  sdev = sqrt(sdev / double(sz));
}
