#pragma once
#include "distpt.h"
#include "regdata.h"
#include <cmath>
#include <limits>
#include <vector>

class MinRes {
  friend class MinGrid;
  friend class MinResTest;

public:
  MinRes(size_t n);
  void resize(size_t n);
  void minimize(const PtVector &dps, const double tau);
  void combine(const MinRes &res1, const MinRes &res2);

private:
  size_t n;                        // "Visible" length allocated for the vectors
  size_t capacity;                 // Length of the vectors
  std::vector<double> epsilon_cum; // Cumulative "budget" spent
  std::vector<double> Sigma_cum;   // Cumulative increase in sum of squares
  std::vector<double> delta_cum;   // Difference at the i-th iteration
  std::vector<double> kappa_cum;   // Coefficient of SS growth in delta
};

class MinGrid {
  friend class MinGridTest; // NOTE: For unit tests only
public:
  MinGrid(size_t m);
  MinGrid(std::initializer_list<double> init_list);
  void grid(const MinRes &res, const double gamma);
  size_t size() const;

  double &operator[](std::size_t index);
  const double &operator[](std::size_t index) const;

private:
  const size_t m;        // Length of the grid
  std::vector<double> g; // Underlying grid
};

class MinTree {
  friend class MinTreeTest; // NOTE: For unit tests only
public:
  MinTree(const RegData &data, const size_t m, const double gamma);

  const std::vector<std::vector<double>> &get_beta_min() const;
  const std::vector<std::vector<double>> &get_beta_max() const;

  void remean(const std::vector<double> &taus);
  void minimize();
  void regress();

private:
  const RegData &data;                       // Regression data
  const size_t g;                            // Number of groups
  const size_t m;                            // Length of the grid
  const double gamma;                        // Step size
  size_t i_tau;                              // First unchanged tau
  std::vector<double> taus;                  // Mean changes in risk
  std::vector<MinRes> ress;                  // Objects to minimize
  std::vector<MinRes> combs;                 // Combining results
  MinGrid grid;                              // Grid of results
  std::vector<std::vector<double>> beta_min; // Minimizing coefficients
  std::vector<std::vector<double>> beta_max; // Maximizing coefficients
};
