#include "sens.h"
#include <testthat.h>

context("Sensitivity Analysis") {
  test_that("Sensitivity Analysis Runs") {
    std::vector<double> e = {0.0, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.4, 0.0, 0.1,
                             0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.1, 0.2,
                             0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    std::vector<double> lwr = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> upr = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    std::vector<double> lwr_iter = {-0.15, -0.25};
    std::vector<double> upr_iter = {0.85, 0.75};
    std::vector<int> grp = {1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    std::vector<bool> a = {false, false, false, false, false, false, false,
                           false, true,  true,  true,  true,  true,  true,
                           true,  true,  true,  true,  true,  true,  true,
                           true,  true,  true,  true,  true,  true,  true};

    double beta = (10.0 / 14.0 - 10.0 / 14.0)
                  // sigma_0   * rho_0 - r_0 + sigma_1   * rho_1 - r_1
                  + (10.0 / 14.0 * 5.1 - 4.5 + 10.0 / 14.0 * 6.5 - 5.5)
                        // \sum R^2 - rho_0^2 / n_0 - rho_1^2 / n_1
                        / (7.14 - 5.1 * 5.1 / 14.0 - 6.5 * 6.5 / 14.0)
                        // rho_1 / n_1 - rho_0 / n_0
                        * (6.5 / 14.0 - 5.1 / 14.0);

    std::vector<std::vector<double>> beta_min;
    std::vector<std::vector<double>> beta_max;

    std::tie<std::vector<std::vector<double>>,
             std::vector<std::vector<double>>>(beta_min, beta_max) =
        sens(e, lwr, upr, lwr_iter, upr_iter, grp, a, 0.1, 0.1, 10, 10, 2);

    expect_true(std::abs(beta_min[0][0] - beta) < 1e-10);
    expect_true(std::abs(beta_max[0][0] - beta) < 1e-10);

    for (size_t i = 0; i < 9; ++i) {
      expect_true(beta_min[0][i + 1] <= beta_min[0][i]);
      expect_true(beta_max[0][i + 1] >= beta_max[0][i]);
    }
  }

  test_that("Results don't depend on the number of threads") {
    std::vector<double> e = {0.0, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.4, 0.0, 0.1,
                             0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.1, 0.2,
                             0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    std::vector<double> lwr = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> upr = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    std::vector<double> lwr_iter = {-0.15, -0.25};
    std::vector<double> upr_iter = {0.85, 0.75};
    std::vector<int> grp = {1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    std::vector<bool> a = {false, false, false, false, false, false, false,
                           false, true,  true,  true,  true,  true,  true,
                           true,  true,  true,  true,  true,  true,  true,
                           true,  true,  true,  true,  true,  true,  true};

    double beta = (10.0 / 14.0 - 10.0 / 14.0)
                  // sigma_0   * rho_0 - r_0 + sigma_1   * rho_1 - r_1
                  + (10.0 / 14.0 * 5.1 - 4.5 + 10.0 / 14.0 * 6.5 - 5.5)
                        // \sum R^2 - rho_0^2 / n_0 - rho_1^2 / n_1
                        / (7.14 - 5.1 * 5.1 / 14.0 - 6.5 * 6.5 / 14.0)
                        // rho_1 / n_1 - rho_0 / n_0
                        * (6.5 / 14.0 - 5.1 / 14.0);

    std::vector<std::vector<double>> beta_min0, beta_max0, beta_min1, beta_max1;

    std::tie<std::vector<std::vector<double>>,
             std::vector<std::vector<double>>>(beta_min0, beta_max0) =
        sens(e, lwr, upr, lwr_iter, upr_iter, grp, a, 0.1, 0.1, 10, 10, 1);
    std::tie<std::vector<std::vector<double>>,
             std::vector<std::vector<double>>>(beta_min1, beta_max1) =
        sens(e, lwr, upr, lwr_iter, upr_iter, grp, a, 0.1, 0.1, 10, 10, 4);

    for (size_t i = 0; i < beta_min0[0].size(); ++i) {
      expect_true(std::abs(beta_min0[0][i] - beta_min1[0][i]) < 1e-10);
      expect_true(std::abs(beta_max0[0][i] - beta_max1[0][i]) < 1e-10);
    }
  }
}
