#pragma once
#include "distpt.h"
#include <cmath>
#include <numeric>
#include <vector>

struct RegData {
  const size_t g;                               // Number of groups
  const std::vector<std::vector<double>> &es;   // Risk vectors
  const std::vector<std::vector<double>> &lwrs; // Lower bounds
  const std::vector<std::vector<double>> &uprs; // Upper bounds
  const std::vector<PtVector> dps;              // DistPt vectors of above data
  const std::vector<size_t> n;                  // Size of each group
  const std::vector<double> sigma;              // Search rate in each group
  const std::vector<double> rho;                // Average risk in each group
  const double numerator; // $$\sum_{j=1}^m \sigma_j t_j - (1 - \sigma_j) r_j$$
  const double ss;        // Sum of squared risk

  RegData(const std::vector<std::vector<double>> &es,
          const std::vector<std::vector<double>> &lwrs,
          const std::vector<std::vector<double>> &uprs);

  template <typename T>
  void regress(const T &grid, const std::vector<double> &taus,
               std::vector<std::vector<double>> &beta_min,
               std::vector<std::vector<double>> &beta_max) const;

private:
  std::vector<PtVector>
  compute_dps(const std::vector<std::vector<double>> &es,
              const std::vector<std::vector<double>> &lwrs,
              const std::vector<std::vector<double>> &uprs) const;
  std::vector<size_t>
  compute_n(const std::vector<std::vector<double>> &es) const;
  std::vector<double>
  compute_sigma(const std::vector<std::vector<double>> &es) const;
  std::vector<double>
  compute_rho(const std::vector<std::vector<double>> &es) const;
  double compute_numerator(const std::vector<std::vector<double>> &es) const;
  double compute_ss(const std::vector<std::vector<double>> &es) const;
};

/*
 * @brief Calculate the regression coefficients and compare to existing records.
 * @param grid Grid of changes in sum of squares.
 * @param taus Vector of group-level average risks.
 * @param beta_min Minimum coefficients.
 * @param beta_max Maximum coefficients.
 * @return None. (The regression coefficients are updated in place.)
 */
template <typename T>
void RegData::regress(const T &grid, const std::vector<double> &taus,
                      std::vector<std::vector<double>> &beta_min,
                      std::vector<std::vector<double>> &beta_max) const {
  // Calculate the updated numerator, denominator, and rho (as `coef`).
  std::vector<double> coef = rho;
  for (size_t j = 0; j < g; ++j)
    coef[j] += (1 - sigma[j]) * taus[j];
  double denominator = ss;
  for (size_t j = 0; j < g; ++j)
    denominator -= n[j] * std::pow(coef[j], 2);
  double numerator_ = numerator;
  for (size_t j = 0; j < g; ++j)
    numerator_ += n[j] * (1 - sigma[j]) * sigma[j] * taus[j];

  // Iterate through the grid and update the betas, continuing if the budget is
  // infeasible.
  double denominator_;
  double scale;
  double beta_0;
  double beta;
  for (size_t i = 0; i < grid.size(); ++i) {
    if (!std::isfinite(grid[i])) {
      continue;
    }
    denominator_ = denominator + grid[i];
    scale = numerator_ / denominator_;
    beta_0 = sigma[0] + scale * coef[0];
    for (size_t j = 1; j < g; ++j) {
      beta = sigma[j] + scale * coef[j] - beta_0;
      if (i == 0) {
        beta_min[j - 1][0] = std::min(beta, beta_min[j - 1][0]);
        beta_max[j - 1][0] = std::max(beta, beta_max[j - 1][0]);
      } else {
        beta_min[j - 1][i] =
            std::min({beta, beta_min[j - 1][i], beta_min[j - 1][i - 1]});
        beta_max[j - 1][i] =
            std::max({beta, beta_max[j - 1][i], beta_max[j - 1][i - 1]});
      }
    }
  }
}
