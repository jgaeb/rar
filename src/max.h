#pragma once
#include "regdata.h"
#include <cmath>
#include <limits>
#include <vector>

class MaxRes {
  friend class MaxGrid;
  friend class MaxResTest; // NOTE: For unit tests only
public:
  MaxRes(size_t n);
  void resize(size_t n);
  void maximize(const std::vector<double> &e, const std::vector<double> &lwr,
                const std::vector<double> &upr, const double tau);

private:
  size_t n;                        // "Visible" length allocated for the vectors
  size_t capacity;                 // Actual length allocated for the vectors
  std::vector<double> epsilon_cum; // Cumulative "budget" spent
  std::vector<double> Sigma_cum;   // Cumulative increase in sum of squares
  std::vector<double> delta_cum;   // Difference at i-th iteration
};

class MaxGrid {
  friend class MaxGridTest; // NOTE: For unit tests only
public:
  MaxGrid(size_t m);
  MaxGrid(std::initializer_list<double> init_list);
  void grid(const MaxRes &res, const double gamma);
  void combine(const MaxGrid &g1, const MaxGrid &g2);
  size_t size() const;

  double &operator[](std::size_t index);
  const double &operator[](std::size_t index) const;

private:
  const size_t m;        // Length of the grid
  std::vector<double> g; // Underlying grid
};

class MaxTree {
  friend class MaxTreeTest; // NOTE: For unit tests only
public:
  MaxTree(const RegData &data, const size_t m, const double gamma);

  const std::vector<std::vector<double>> &get_beta_min() const;
  const std::vector<std::vector<double>> &get_beta_max() const;

  void remean(const std::vector<double> &taus);
  void maximize();
  void regress();

private:
  const RegData &data;                       // Regression data
  const size_t g;                            // Number of groups
  const size_t m;                            // Length of the grid
  const double gamma;                        // Step size
  size_t i_tau;                              // First unchanged tau
  std::vector<double> taus;                  // Mean changes in risk
  std::vector<MaxGrid> res_grids;            // Max results
  std::vector<MaxGrid> comb_grids;           // Combining results
  std::vector<MaxRes> ress;                  // Objects to maximize
  std::vector<std::vector<double>> beta_min; // Minimizing coefficients
  std::vector<std::vector<double>> beta_max; // Maximizing coefficients
};
