#include "max.h"
#include <testthat.h>

/*
 * @brief Test class for MaxRes to expose private members.
 */
class MaxResTest {
public:
  MaxResTest(const MaxRes &res) : res(res){};
  const size_t &get_n() const { return res.n; };
  const size_t &get_capacity() const { return res.capacity; };
  const std::vector<double> &get_epsilon_cum() const {
    return res.epsilon_cum;
  };
  const std::vector<double> &get_Sigma_cum() const { return res.Sigma_cum; };
  const std::vector<double> &get_delta_cum() const { return res.delta_cum; };

private:
  const MaxRes &res;
};

context("Maximization Results (MaxRes)") {
  auto e0 = std::vector<double>{0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                                0.6, 0.7, 0.8, 0.9, 1.0};
  auto lwr0 = std::vector<double>(11, 0.0);
  auto upr0 = std::vector<double>(11, 1.0);
  MaxRes res0_(11);
  MaxResTest res0(res0_);
  res0_.maximize(e0, lwr0, upr0, 0.0);

  test_that("Base Case: Epsilon_cum is increasing") {
    for (int i = 0; i < 10; i++) {
      expect_true(res0.get_epsilon_cum()[i] <= res0.get_epsilon_cum()[i + 1]);
    }
  }

  test_that("Base Case: Sigma_cum is increasing") {
    for (int i = 0; i < 10; i++) {
      expect_true(res0.get_Sigma_cum()[i] <= res0.get_Sigma_cum()[i + 1]);
    }
  }

  test_that("Base Case: Vectors have the correct length") {
    expect_true(res0.get_n() == 11);
    expect_true(res0.get_epsilon_cum().size() == res0.get_capacity());
    expect_true(res0.get_Sigma_cum().size() == res0.get_capacity());
    expect_true(res0.get_delta_cum().size() == res0.get_capacity());
  }

  test_that("Base Case: Cumulative sums are correct") {
    expect_true(std::abs(res0.get_epsilon_cum()[10] - 2) < 1e-10);
    expect_true(std::abs(res0.get_Sigma_cum()[10] - 1.4) < 1e-10);
  }

  // Test base case where mean is adjusted upward
  res0_.maximize(e0, lwr0, upr0, 0.5);

  test_that("Mean increased: Vectors have the correct length") {
    expect_true(res0.get_n() == 1);
    expect_true(res0.get_epsilon_cum().size() == res0.get_capacity());
    expect_true(res0.get_Sigma_cum().size() == res0.get_capacity());
    expect_true(res0.get_delta_cum().size() == res0.get_capacity());
  }

  test_that("Mean increased: Cumulative sums are correct") {
    expect_true(std::abs(res0.get_epsilon_cum()[0] - 5.5) < 1e-10);
    expect_true(std::abs(res0.get_Sigma_cum()[0] - 7.15) < 1e-10);
  }

  // Test base case where mean is adjusted upward partially
  res0_.maximize(e0, lwr0, upr0, 0.25);
  test_that("Mean partially increased: Vectors have the correct length") {
    expect_true(res0.get_n() == 4);
    expect_true(res0.get_epsilon_cum().size() == res0.get_capacity());
    expect_true(res0.get_Sigma_cum().size() == res0.get_capacity());
    expect_true(res0.get_delta_cum().size() == res0.get_capacity());
  }

  test_that("Mean partially increased: Cumulative sums are correct") {
    expect_true(std::abs(res0.get_epsilon_cum()[0] - 2.75) < 1e-10);
    expect_true(std::abs(res0.get_Sigma_cum()[0] - 4.1025) < 1e-10);
    expect_true(std::abs(res0.get_epsilon_cum()[3] - 2.95) < 1e-10);
    expect_true(std::abs(res0.get_Sigma_cum()[3] - 4.2125) < 1e-10);
  }

  // Test base case where mean is adjusted downward
  res0_.maximize(e0, lwr0, upr0, -0.5);

  test_that("Mean decreased: Vectors have the correct length") {
    expect_true(res0.get_n() == 1);
    expect_true(res0.get_epsilon_cum().size() == res0.get_capacity());
    expect_true(res0.get_Sigma_cum().size() == res0.get_capacity());
    expect_true(res0.get_delta_cum().size() == res0.get_capacity());
  }

  test_that("Mean decreased: Cumulative sums are correct") {
    expect_true(std::abs(res0.get_epsilon_cum()[0] - 5.5) < 1e-10);
    expect_true(std::abs(res0.get_Sigma_cum()[0] + 3.85) < 1e-10);
  }

  // Test base case where mean is adjusted downward partially
  res0_.maximize(e0, lwr0, upr0, -0.25);

  test_that("Mean partially decreased: Vectors have the correct length") {
    expect_true(res0.get_n() == 4);
    expect_true(res0.get_epsilon_cum().size() == res0.get_capacity());
    expect_true(res0.get_Sigma_cum().size() == res0.get_capacity());
    expect_true(res0.get_delta_cum().size() == res0.get_capacity());
  }

  test_that("Mean partially decreased: Cumulative sums are correct") {
    expect_true(std::abs(res0.get_epsilon_cum()[0] - 2.75) < 1e-10);
    expect_true(std::abs(res0.get_Sigma_cum()[0] + 1.3975) < 1e-10);
    expect_true(std::abs(res0.get_epsilon_cum()[3] - 2.95) < 1e-10);
    expect_true(std::abs(res0.get_Sigma_cum()[3] + 1.2875) < 1e-10);
  }

  // Test alternative with known solution
  auto e1 = std::vector<double>{0.0, 0.25, 0.5, 0.75, 1.0};
  auto lwr1 = std::vector<double>(5, 0.0);
  auto upr1 = std::vector<double>(5, 1.0);
  MaxRes res1_(5);
  res1_.maximize(e1, lwr1, upr1, 0.0);
  MaxResTest res1(res1_);

  test_that("Alternate Case: Vectors are correct") {
    auto epsilon_cum = std::vector<double>{0.0, 0.0, 0.0, 1.0 / 2.0, 1.0 / 2.0};
    auto Sigma_cum = std::vector<double>{0.0, 0.0, 0.0, 3.0 / 8.0, 3.0 / 8.0};
    auto delta_cum =
        std::vector<double>{1.0, 3.0 / 4.0, 1.0 / 2.0, 1.0 / 2.0, 0.0};
    expect_true(res1.get_n() == 5);
    for (int i = 0; i < 5; i++) {
      expect_true(std::abs(res1.get_epsilon_cum()[i] - epsilon_cum[i]) < 1e-10);
      expect_true(std::abs(res1.get_Sigma_cum()[i] - Sigma_cum[i]) < 1e-10);
      expect_true(std::abs(res1.get_delta_cum()[i] - delta_cum[i]) < 1e-10);
    }
  }

  // Test with non-trivial upper and lower bounds
  auto e2 = std::vector<double>{1.0 / 14.0, 4.0 / 14.0, 10.0 / 14.0};
  auto lwr2 = std::vector<double>{0.0 / 14.0, 2.0 / 14.0, 6.0 / 14.0};
  auto upr2 = std::vector<double>{2.0 / 14.0, 6.0 / 14.0, 14.0 / 14.0};
  MaxRes res2_(3);
  res2_.maximize(e2, lwr2, upr2, 1.0 / 42.0);
  MaxResTest res2(res2_);

  test_that("Bounds Case: Vectors are correct") {
    auto epsilon_cum = std::vector<double>{1.0 / 14.0, 3.0 / 14.0, 7.0 / 14.0};
    auto Sigma_cum =
        std::vector<double>{3.0 / 28.0, 43.0 / 196.0, 83.0 / 196.0};
    auto delta_cum = std::vector<double>{10.0 / 14.0, 8.0 / 14.0, 0.0};

    expect_true(res2.get_n() == 3);
    for (int i = 0; i < 3; i++) {
      expect_true(std::abs(res2.get_epsilon_cum()[i] - epsilon_cum[i]) < 1e-10);
      expect_true(std::abs(res2.get_Sigma_cum()[i] - Sigma_cum[i]) < 1e-10);
      expect_true(std::abs(res2.get_delta_cum()[i] - delta_cum[i]) < 1e-10);
    }
  }
}

class MaxGridTest {
public:
  MaxGridTest(const MaxGrid &g) : g(g){};
  const size_t &get_m() const { return g.m; };
  const std::vector<double> &get_g() const { return g.g; };

private:
  const MaxGrid &g;
};

context("Maximization Grids (MaxGrid)") {
  auto e = std::vector<double>{0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                               0.6, 0.7, 0.8, 0.9, 1.0};
  auto lwr = std::vector<double>(11, 0.0);
  auto upr = std::vector<double>(11, 1.0);
  MaxRes res(11);
  res.maximize(e, lwr, upr, 0.0);

  test_that("Grid construction works correctly") {
    MaxGrid g_(11);
    g_.grid(res, 0.2);
    MaxGridTest g(g_);
    std::vector<double> expected_g = {0,    0.18, 0.32, 0.50, 0.60, 0.74,
                                      0.92, 0.98, 1.08, 1.22, 1.40};
    for (int i = 0; i < 11; ++i) {
      expect_true(std::abs(g.get_g()[i] - expected_g[i]) < 1e-10);
    }
  }

  test_that("Grid combination works correctly") {
    MaxGrid g1 = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    MaxGrid g2 = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    MaxGrid g3 = {-INFINITY, -INFINITY, 0.04, 0.09, 0.16, 0.25,
                  0.36,      0.49,      0.64, 0.81, 1.0};
    MaxGrid g(11);
    g.combine(g1, g2);
    for (int i = 0; i < 11; ++i) {
      expect_true(abs(g1[i] - g[i]) < 1e-10);
    }
    g.combine(g1, g3);
    for (int i = 0; i < 11; ++i) {
      expect_true(g1[i] >= g[i]);
    }

    MaxGrid g4 = {0.0, 0.2, 0.6, 0.9, 1.0};
    MaxGrid g5 = {-INFINITY, -INFINITY, 0.4, 1.0, 1.3};
    MaxGrid g_prime(5);

    g_prime.combine(g4, g5);
    std::vector<double> g_prime_expected = {-INFINITY, -INFINITY, 0.4, 1.0,
                                            1.3};
    for (int i = 0; i < 2; ++i) { // NOTE: -INFINITY + INFINITY is NaN
      expect_true(g_prime_expected[i] == g_prime[i]);
    }
    for (int i = 2; i < 5; ++i) {
      expect_true(std::abs(g_prime_expected[i] - g_prime[i]) < 1e-10);
    }
  }

  test_that("Grids formed from MaxRes objects of length one work correctly") {
    MaxGrid g(21);
    res.maximize(e, lwr, upr, -0.5);
    g.grid(res, 0.5);
    std::vector<double> g_expected = {
        -INFINITY, -INFINITY, -INFINITY, -INFINITY, -INFINITY, -INFINITY,
        -INFINITY, -INFINITY, -INFINITY, -INFINITY, -INFINITY, -INFINITY,
        -3.85,     -3.85,     -3.85,     -3.85,     -3.85,     -3.85,
        -3.85,     -3.85,     -3.85};
    for (int i = 0; i < 12; ++i) { // NOTE: -INFINITY + INFINITY is NaN
      expect_true(g_expected[i] == g[i]);
    }
    for (int i = 12; i < 21; ++i) {
      expect_true(std::abs(g_expected[i] - g[i]) < 1e-10);
    }
    res.maximize(e, lwr, upr, 0.5);
    g.grid(res, 0.5);
    g_expected = {-INFINITY, -INFINITY, -INFINITY, -INFINITY, -INFINITY,
                  -INFINITY, -INFINITY, -INFINITY, -INFINITY, -INFINITY,
                  -INFINITY, -INFINITY, 7.15,      7.15,      7.15,
                  7.15,      7.15,      7.15,      7.15,      7.15,
                  7.15};
    for (int i = 0; i < 12; ++i) { // NOTE: -INFINITY + INFINITY is NaN
      expect_true(g_expected[i] == g[i]);
    }
    for (int i = 12; i < 21; ++i) {
      expect_true(std::abs(g_expected[i] - g[i]) < 1e-10);
    }
  }

  test_that("Grids formed with zero budget work correctly") {
    MaxGrid g(1);
    res.maximize(e, lwr, upr, 0.0);
    g.grid(res, 0.0);
    expect_true(g[0] == 0.0);
    res.maximize(e, lwr, upr, 0.5);
    g.grid(res, 0.0);
    expect_true(g[0] == -INFINITY);
    res.maximize(e, lwr, upr, -0.5);
    g.grid(res, 0.0);
    expect_true(g[0] == -INFINITY);
  }

  test_that("Grids formed when the mean is changed begin with inf.") {
    MaxGrid g(21);
    res.maximize(e, lwr, upr, 1e-6);
    g.grid(res, 0.5);
    expect_true(g[0] == -INFINITY);
    res.maximize(e, lwr, upr, -1e-6);
    g.grid(res, 0.5);
    expect_true(g[0] == -INFINITY);
  }
}

class MaxTreeTest {
public:
  MaxTreeTest(const MaxTree &t) : t(t){};
  const size_t &get_m() const { return t.m; };
  const double &get_gamma() const { return t.gamma; };
  const size_t &get_i_tau() const { return t.i_tau; };
  const std::vector<double> &get_taus() const { return t.taus; };
  const std::vector<MaxGrid> &get_res_grids() const { return t.res_grids; };
  const std::vector<MaxGrid> &get_comb_grids() const { return t.comb_grids; };
  const std::vector<MaxRes> &get_ress() const { return t.ress; };
  const std::vector<std::vector<double>> &get_beta_min() const {
    return t.beta_min;
  };
  const std::vector<std::vector<double>> &get_beta_max() const {
    return t.beta_max;
  };

private:
  const MaxTree &t;
};

context("Maximization trees (MaxTree)") {
  test_that("Dimensions of underlying strata and grids are correct") {
    // First stratum has three people, second has two, third has one, fourth has
    // seven.
    auto es = std::vector<std::vector<double>>{
        {0.0, 0.1, 0.2},
        {0.3, 0.4},
        {0.5},
        {0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9}};
    auto lwrs =
        std::vector<std::vector<double>>{{0.0, 0.0, 0.0},
                                         {0.0, 0.0},
                                         {0.0},
                                         {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    auto uprs =
        std::vector<std::vector<double>>{{1.0, 1.0, 1.0},
                                         {1.0, 1.0},
                                         {1.0},
                                         {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};
    RegData data(es, lwrs, uprs);
    MaxTree t_(data, 11, 0.1);
    MaxTreeTest t(t_);

    expect_true(t.get_m() == 11);
    expect_true(t.get_gamma() == 0.1);
    expect_true(t.get_i_tau() == 1);
    expect_true(t.get_taus().size() == 2);
    expect_true(t.get_res_grids().size() == 4);
    expect_true(t.get_comb_grids().size() == 3);
    expect_true(t.get_ress().size() == 4);
    expect_true(t.get_beta_min().size() == 1);
    expect_true(t.get_beta_min()[0].size() == 11);
    expect_true(t.get_beta_max().size() == 1);
    expect_true(t.get_beta_max()[0].size() == 11);
  }

  test_that("Optimization is correct") {
    std::vector<double> e{0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                          0.6, 0.7, 0.8, 0.9, 1.0};
    const std::vector<decltype(e)> es{e, e, e, e};
    decltype(es) lwrs(4, decltype(e)(11, 0.0));
    decltype(es) uprs(4, decltype(e)(11, 1.0));
    RegData data(es, lwrs, uprs);
    MaxTree t_(data, 11, 0.2);
    t_.maximize();
    MaxTreeTest t(t_);

    expect_true(t.get_taus().size() == 2);
    std::vector<double> expected_taus = {0.0, 0.0};
    for (int i = 0; i < 2; ++i)
      expect_true(std::abs(t.get_taus()[i] - expected_taus[i]) < 1e-10);
    std::vector<double> expected_res_grid = {0.0,  0.18, 0.32, 0.50, 0.60, 0.74,
                                             0.92, 0.98, 1.08, 1.22, 1.40};
    for (int i = 0; i < 11; ++i)
      for (int j = 0; j < 3; ++j)
        expect_true(std::abs(t.get_res_grids()[j][i] - expected_res_grid[i]) <
                    1e-10);
    std::vector<double> expected_comb_grid2 = {
        0.0, 0.18, 0.36, 0.50, 0.68, 0.82, 1.0, 1.1, 1.24, 1.42, 1.52};
    for (int i = 0; i < 11; ++i)
      expect_true(std::abs(t.get_comb_grids()[2][i] - expected_comb_grid2[i]) <
                  1e-10);
    std::vector<double> expected_comb_grid1 = {
        0.0, 0.18, 0.36, 0.54, 0.68, 0.86, 1.0, 1.18, 1.32, 1.5, 1.6};
    for (int i = 0; i < 11; ++i)
      expect_true(std::abs(t.get_comb_grids()[1][i] - expected_comb_grid1[i]) <
                  1e-10);
    std::vector<double> expected_comb_grid0 = {
        0.0, 0.18, 0.36, 0.54, 0.72, 0.86, 1.04, 1.18, 1.36, 1.5, 1.68};
    for (int i = 0; i < 11; ++i)
      expect_true(std::abs(t.get_comb_grids()[0][i] - expected_comb_grid0[i]) <
                  1e-10);
  }

  test_that("Regression coefficients are computed correctly.") {
    auto es = std::vector<std::vector<double>>{
        {0.0, 0.1, 0.2, 0.3},
        {0.1, 0.2, 0.3, 0.4},
        {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},
        {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}};
    auto lwrs = std::vector<std::vector<double>>{
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    auto uprs = std::vector<std::vector<double>>{
        {1.0, 1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};
    RegData data(es, lwrs, uprs);
    MaxTree t_(data, 10, 0.6);
    t_.maximize();
    t_.regress();
    MaxTreeTest t(t_);
    //             sigma_1     - sigma_0
    double beta0 = (10.0 / 14.0 - 10.0 / 14.0)
                   // sigma_0   * rho_0 - r_0 + sigma_1   * rho_1 - r_1
                   + (10.0 / 14.0 * 5.1 - 4.5 + 10.0 / 14.0 * 6.5 - 5.5)
                         // \sum R^2 - rho_0^2 / n_0 - rho_1^2 / n_1
                         / (7.14 - 5.1 * 5.1 / 14.0 - 6.5 * 6.5 / 14.0)
                         // rho_1 / n_1 - rho_0 / n_0
                         * (6.5 / 14.0 - 5.1 / 14.0);
    expect_true(std::abs(t.get_beta_min()[0][0] - beta0) < 1e-10);
    expect_true(std::abs(t.get_beta_max()[0][0] - beta0) < 1e-10);
    for (size_t i = 0; i < 9; ++i) {
      expect_true(t.get_beta_min()[0][i + 1] <= t.get_beta_min()[0][i]);
      expect_true(t.get_beta_max()[0][i + 1] >= t.get_beta_max()[0][i]);
    }

    // Increase the mean of the first stratum by 0.15
    t_.remean(std::vector<double>{0.15, 0.0});
    t_.maximize();
    t_.regress();

    // Recalculate beta, noting that the first stratum has 0.15 added to it,
    // which increases rho0 by 0.15 * 4 = 0.6, and increases the sum of squares
    // by 0.72
    //              sigma_1     - sigma_0
    double beta1 = (10.0 / 14.0 - 10.0 / 14.0)
                   // sigma_0   * rho_0 - r_0 + sigma_1   * rho_1 - r_1
                   + (10.0 / 14.0 * 5.7 - 4.5 + 10.0 / 14.0 * 6.5 - 5.5)
                         // \sum R^2 - rho_0^2 / n_0 - rho_1^2 / n_1
                         / (7.86 - 5.7 * 5.7 / 14.0 - 6.5 * 6.5 / 14.0)
                         // rho_1 / n_1 - rho_0 / n_0
                         * (6.5 / 14.0 - 5.7 / 14.0);
    expect_true(std::abs(t.get_beta_min()[0][0] - beta0) < 1e-10);
    expect_true(std::abs(t.get_beta_max()[0][0] - beta0) < 1e-10);
    expect_true(std::abs(t.get_beta_max()[0][1] - beta1) < 1e-10);
    for (size_t i = 0; i < 9; ++i) {
      expect_true(t.get_beta_min()[0][i + 1] <= t.get_beta_min()[0][i]);
      expect_true(t.get_beta_max()[0][i + 1] >= t.get_beta_max()[0][i]);
    }
  }

  test_that("Remeaning works correctly.") {
    auto es = std::vector<std::vector<double>>{
        {0.0, 0.1, 0.2, 0.3},
        {0.1, 0.2, 0.3, 0.4},
        {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},
        {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}};
    auto lwrs = std::vector<std::vector<double>>{
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    auto uprs = std::vector<std::vector<double>>{
        {1.0, 1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};
    RegData data(es, lwrs, uprs);
    MaxTree t_(data, 10, 0.2);
    MaxTreeTest t(t_);
    t_.remean(std::vector<double>{0.0, 0.1});
    expect_true(t.get_taus()[0] == 0.0);
    expect_true(t.get_taus()[1] == 0.1);
    expect_true(t.get_i_tau() == 1);
    t_.remean(std::vector<double>{0.1, 0.1});
    expect_true(t.get_taus()[0] == 0.1);
    expect_true(t.get_taus()[1] == 0.1);
    expect_true(t.get_i_tau() == 0);
  }
}
