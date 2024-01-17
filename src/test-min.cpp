#include "min.h"
#include <iostream>
#include <testthat.h>

/*
 * @brief Test class for MinRes to expose private members
 */
class MinResTest {
public:
  MinResTest(const MinRes &res) : res(res) {}
  const size_t &get_n() { return res.n; }
  const size_t &get_capacity() { return res.capacity; }
  const std::vector<double> &get_epsilon_cum() { return res.epsilon_cum; }
  const std::vector<double> &get_Sigma_cum() { return res.Sigma_cum; }
  const std::vector<double> &get_delta_cum() { return res.delta_cum; }
  const std::vector<double> &get_kappa_cum() { return res.kappa_cum; }

private:
  const MinRes &res;
};

context("Minimization Results (MinRes)") {
  auto e0 = std::vector<double>{0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                                0.6, 0.7, 0.8, 0.9, 1.0};
  auto lwr0 = std::vector<double>(11, 0.0);
  auto upr0 = std::vector<double>(11, 1.0);
  auto dps0 = PtVector(e0, lwr0, upr0);
  MinRes res0_(33);
  res0_.minimize(dps0, 0.0);
  MinResTest res0(res0_);

  test_that("Base Case: Epsilon_cum is increasing") {
    for (int i = 0; i < 32; i++) {
      expect_true(res0.get_epsilon_cum()[i] <= res0.get_epsilon_cum()[i + 1]);
    }
  }

  test_that("Base Case: Sigma_cum is decreasing") {
    for (int i = 0; i < 32; i++) {
      expect_true(res0.get_Sigma_cum()[i] >= res0.get_Sigma_cum()[i + 1]);
    }
  }

  test_that("Base Caes: Vectors have the correct length") {
    expect_true(res0.get_n() == 33);
    expect_true(res0.get_epsilon_cum().size() == res0.get_capacity());
    expect_true(res0.get_Sigma_cum().size() == res0.get_capacity());
    expect_true(res0.get_delta_cum().size() == res0.get_capacity());
    expect_true(res0.get_kappa_cum().size() == res0.get_capacity());
  }

  test_that("Base Case: Cumulative sums are correct") {
    expect_true(std::abs(res0.get_epsilon_cum()[32] - 3) < 1e-10);
    expect_true(std::abs(res0.get_Sigma_cum()[32] + 1.1) < 1e-10);
  }

  // Test base case where mean is adjusted upward
  res0_.minimize(dps0, 0.5);

  test_that("Mean increased: Vectors have the correct length") {
    expect_true(res0.get_epsilon_cum().size() == res0.get_capacity());
    expect_true(res0.get_Sigma_cum().size() == res0.get_capacity());
    expect_true(res0.get_delta_cum().size() == res0.get_capacity());
    expect_true(res0.get_kappa_cum().size() == res0.get_capacity());
  }

  test_that("Mean increased: Cumulative sums are correct") {
    expect_true(std::abs(res0.get_epsilon_cum()[res0.get_n() - 1] - 5.5) <
                1e-10);
    expect_true(std::abs(res0.get_Sigma_cum()[res0.get_n() - 1] - 7.15) <
                1e-10);
  }

  // Test base case where mean is adjusted downward
  res0_.minimize(dps0, -0.5);

  test_that("Mean decreased: Vectors have the correct length") {
    expect_true(res0.get_epsilon_cum().size() == res0.get_capacity());
    expect_true(res0.get_Sigma_cum().size() == res0.get_capacity());
    expect_true(res0.get_delta_cum().size() == res0.get_capacity());
    expect_true(res0.get_kappa_cum().size() == res0.get_capacity());
  }

  test_that("Mean decreased: Cumulative sums are correct") {
    expect_true(std::abs(res0.get_epsilon_cum()[res0.get_n() - 1] - 5.5) <
                1e-10);
    expect_true(std::abs(res0.get_Sigma_cum()[res0.get_n() - 1] + 3.85) <
                1e-10);
  }

  // Test alternative with known solution
  auto e1 = std::vector<double>{0.0, 0.25, 0.5, 0.75, 1.0};
  auto lwr1 = std::vector<double>(5, 0.0);
  auto upr1 = std::vector<double>(5, 1.0);
  auto dps1 = PtVector(e1, lwr1, upr1);
  MinRes res1_(15);
  res1_.minimize(dps1, 0.0);
  MinResTest res1(res1_);

  test_that("Alternate Case: Vectors are correct") {
    auto epsilon_cum = std::vector<double>{0.0, 0.5, 1.5};
    auto Sigma_cum = std::vector<double>{0.0, -3.0 / 8.0, -5.0 / 8.0};
    auto delta_cum = std::vector<double>{1.0, 0.5, 0.0};
    auto kappa_cum = std::vector<double>{0.5, 1.0, 0.0};
    expect_true(res1.get_n() == 15);
    for (int i = 0; i < 3; i++) {
      expect_true(std::abs(res1.get_epsilon_cum()[2 * i + 10] -
                           epsilon_cum[i]) < 1e-10);
      expect_true(std::abs(res1.get_Sigma_cum()[2 * i + 10] - Sigma_cum[i]) <
                  1e-10);
      expect_true(std::abs(res1.get_delta_cum()[2 * i + 10] - delta_cum[i]) <
                  1e-10);
      expect_true(std::abs(res1.get_kappa_cum()[2 * i + 10] - kappa_cum[i]) <
                  1e-10);
    }
  }

  // Test with non-trivial upper and lower bounds
  auto e2 = std::vector<double>{1.0 / 14.0, 4.0 / 14.0, 10.0 / 14.0};
  auto lwr2 = std::vector<double>{0.0 / 14.0, 2.0 / 14.0, 6.0 / 14.0};
  auto upr2 = std::vector<double>{2.0 / 14.0, 6.0 / 14.0, 14.0 / 14.0};
  auto dps2 = PtVector(e2, lwr2, upr2);
  MinRes res2_(9);
  res2_.minimize(dps2, 1.0 / 42.0);
  MinResTest res2(res2_);

  test_that("Bounds Case: Vectors are correct") {
    auto epsilon_cum = std::vector<double>{1.0 / 14.0, 5.0 / 14.0};
    auto Sigma_cum = std::vector<double>{3.0 / 196.0, -13.0 / 196.0};
    auto delta_cum = std::vector<double>{3.0 / 7.0, 0.0};
    auto kappa_cum = std::vector<double>{0.5, 0.0};

    std::vector<size_t> offset = {4, 1};
    for (int i = 0; i < 2; i++) {
      size_t j = res2.get_n() - offset[i];
      expect_true(std::abs(res2.get_epsilon_cum()[j] - epsilon_cum[i]) < 1e-10);
      expect_true(std::abs(res2.get_Sigma_cum()[j] - Sigma_cum[i]) < 1e-10);
      expect_true(std::abs(res2.get_delta_cum()[j] - delta_cum[i]) < 1e-10);
      expect_true(std::abs(res2.get_kappa_cum()[j] - kappa_cum[i]) < 1e-10);
    }
  }

  // Test with mean shifted up only partially
  auto e3 = std::vector<double>{0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                                0.6, 0.7, 0.8, 0.9, 1.0};
  auto lwr3 = std::vector<double>(11, 0.0);
  auto upr3 = std::vector<double>(11, 1.0);
  auto dps3 = PtVector(e3, lwr3, upr3);
  auto res3_ = MinRes(33);
  auto res3 = MinResTest(res3_);
  res3_.minimize(dps3, 0.25);

  test_that("Mean shifted up only partially: Vectors have the correct length") {
    expect_true(res3.get_epsilon_cum().size() == res3.get_capacity());
    expect_true(res3.get_Sigma_cum().size() == res3.get_capacity());
    expect_true(res3.get_delta_cum().size() == res3.get_capacity());
    expect_true(res3.get_kappa_cum().size() == res3.get_capacity());
    expect_true(res3.get_n() == 16);
  }

  test_that("Mean shifted up only partially: Cumulative sums are correct") {
    expect_true(std::abs(res3.get_epsilon_cum()[0] - 2.75) < 1e-10);
    expect_true(std::abs(res3.get_Sigma_cum()[0] - 6861.0 / 2800.0) < 1e-10);
    expect_true(std::abs(res3.get_epsilon_cum()[15] - 3.65) < 1e-10);
    expect_true(std::abs(res3.get_Sigma_cum()[15] - 935.0 / 400.0) < 1e-10);
  }

  // Test with mean shifted down only partially
  res3_.minimize(dps3, -0.25);

  test_that(
      "Mean shifted down only partially: Vectors have the correct length") {
    expect_true(res3.get_epsilon_cum().size() == res3.get_capacity());
    expect_true(res3.get_Sigma_cum().size() == res3.get_capacity());
    expect_true(res3.get_delta_cum().size() == res3.get_capacity());
    expect_true(res3.get_kappa_cum().size() == res3.get_capacity());
    expect_true(res3.get_n() == 16);
  }

  test_that("Mean shifted down only partially: Cumulative sums are correct") {
    expect_true(std::abs(res3.get_epsilon_cum()[0] - 2.75) < 1e-10);
    expect_true(std::abs(res3.get_Sigma_cum()[0] - (-8539.0 / 2800.0)) < 1e-10);
    expect_true(std::abs(res3.get_epsilon_cum()[15] - 3.65) < 1e-10);
    expect_true(std::abs(res3.get_Sigma_cum()[15] - (-3.1625)) < 1e-10);
  }

  // Test MinRes combination
  auto e4 = std::vector<double>{0.0, 0.25, 0.5, 0.75, 1.0};
  auto lwr4 = std::vector<double>(5, 0.0);
  auto upr4 = std::vector<double>(5, 1.0);
  auto dps4 = PtVector(e4, lwr4, upr4);
  auto res4_ = MinRes(15);
  res4_.minimize(dps4, 0.0);

  auto e5 = std::vector<double>{0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0};
  auto lwr5 = std::vector<double>(4, 0.0);
  auto upr5 = std::vector<double>(4, 1.0);
  auto dps5 = PtVector(e5, lwr5, upr5);
  auto res5_ = MinRes(12);
  res5_.minimize(dps5, 0.0);

  auto res6_ = MinRes(27);
  auto res6 = MinResTest(res6_);
  res6_.combine(res4_, res5_);

  test_that("Combined MinRes has the correct length and values.") {
    auto epsilon_cum = std::vector<double>{0.0, 1.0, 1.5, 17.0 / 6.0};
    auto Sigma_cum =
        std::vector<double>{0.0, -0.75, -23.0 / 24.0, -85.0 / 72.0};
    auto delta_cum = std::vector<double>{1.0, 0.5, 1.0 / 3.0, 0.0};
    auto kappa_cum = std::vector<double>{1.0, 1.5, 2.0, 0.0};
    expect_true(res6.get_n() == 4);
    for (int i = 0; i < 4; i++) {
      expect_true(std::abs(res6.get_epsilon_cum()[i] - epsilon_cum[i]) < 1e-10);
      expect_true(std::abs(res6.get_Sigma_cum()[i] - Sigma_cum[i]) < 1e-10);
      expect_true(std::abs(res6.get_delta_cum()[i] - delta_cum[i]) < 1e-10);
      expect_true(std::abs(res6.get_kappa_cum()[i] - kappa_cum[i]) < 1e-10);
    }
  }
}

class MinGridTest {
public:
  MinGridTest(const MinGrid &g) : g(g){};
  const size_t &get_m() const { return g.m; };
  const std::vector<double> &get_g() const { return g.g; };

private:
  const MinGrid &g;
};

context("Minimization Grids (MinGrid)") {
  auto e = std::vector<double>{0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                               0.6, 0.7, 0.8, 0.9, 1.0};
  auto lwr = std::vector<double>(11, 0.0);
  auto upr = std::vector<double>(11, 1.0);
  auto dps = PtVector(e, lwr, upr);
  MinRes res(33);
  res.minimize(dps, 0.0);

  test_that("Grid construction works correctly") {
    MinGrid g_(26);
    g_.grid(res, 0.12);
    MinGridTest g(g_);
    // Calculate changes every 0.001 change in the threshold
    std::vector<double> Sigma_1, Sigma_2, Sigma_3, Sigma_4, Sigma_5;
    for (int k = 0; k <= 100; ++k) {
      double i = 0.001 * k;
      double j = 1.0 - i;
      Sigma_1.push_back(std::pow(i, 2) + std::pow(j, 2) - 1);
    }
    for (int k = 0; k <= 100; ++k) {
      double i = 0.001 * k + 0.1;
      double j = 1.0 - i;
      Sigma_2.push_back(2 * (std::pow(i, 2) + std::pow(j, 2)) -
                        (1 + std::pow(0.9, 2) + std::pow(0.1, 2)));
    }
    for (int k = 0; k <= 100; ++k) {
      double i = 0.001 * k + 0.2;
      double j = 1.0 - i;
      Sigma_3.push_back(3 * (std::pow(i, 2) + std::pow(j, 2)) -
                        (1 + std::pow(0.9, 2) + std::pow(0.8, 2) +
                         std::pow(0.1, 2) + std::pow(0.2, 2)));
    }
    for (int k = 0; k <= 100; ++k) {
      double i = 0.001 * k + 0.3;
      double j = 1.0 - i;
      Sigma_4.push_back(4 * (std::pow(i, 2) + std::pow(j, 2)) -
                        (1 + std::pow(0.9, 2) + std::pow(0.8, 2) +
                         std::pow(0.7, 2) + std::pow(0.1, 2) +
                         std::pow(0.2, 2) + std::pow(0.3, 2)));
    }
    for (int k = 0; k <= 100; ++k) {
      double i = 0.001 * k + 0.4;
      double j = 1.0 - i;
      Sigma_5.push_back(5 * (std::pow(i, 2) + std::pow(j, 2)) -
                        (1 + std::pow(0.9, 2) + std::pow(0.8, 2) +
                         std::pow(0.7, 2) + std::pow(0.6, 2) +
                         std::pow(0.1, 2) + std::pow(0.2, 2) +
                         std::pow(0.3, 2) + std::pow(0.4, 2)));
    }

    // Each step consumes 0.06 of the budget. That means every 60 indices in the
    // first "round," every 30 in the second, every 20 in the third, every 15 in
    // the fourth, and every 12 in the fifth.
    std::vector<double> expected_g;
    for (size_t i = 0; i < Sigma_1.size(); i += 60) {
      expected_g.push_back(Sigma_1[i]);
    }
    for (size_t i = 10; i < Sigma_2.size(); i += 30) {
      expected_g.push_back(Sigma_2[i]);
    }
    for (size_t i = 20; i < Sigma_3.size(); i += 20) {
      expected_g.push_back(Sigma_3[i]);
    }
    for (size_t i = 15; i < Sigma_4.size(); i += 15) {
      expected_g.push_back(Sigma_4[i]);
    }
    for (size_t i = 4; i < Sigma_5.size(); i += 12) {
      expected_g.push_back(Sigma_5[i]);
    }
    for (int i = 0; i < 26; ++i) {
      expect_true(std::abs(g.get_g()[i] - expected_g[i]) < 1e-10);
    }
  }

  test_that("Grids formed from MinRes objects of length one work correctly") {
    MinGrid g(21);
    std::vector<double> e1 = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    std::vector<double> l1 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> u1 = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    auto dps1 = PtVector(e1, l1, u1);
    res.minimize(dps1, -0.5);
    g.grid(res, 0.5);
    std::vector<double> g_expected = {
        INFINITY, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY,
        INFINITY, INFINITY, INFINITY, -2.85,    -2.85,    -2.85,    -2.85,
        -2.85,    -2.85,    -2.85,    -2.85,    -2.85,    -2.85,    -2.85};
    for (int i = 0; i < 10; ++i) {
      expect_true(g_expected[i] == g[i]);
    }
    for (int i = 10; i < 21; ++i) {
      expect_true(std::abs(g_expected[i] - g[i]) < 1e-10);
    }
    res.minimize(dps1, 0.5);
    g.grid(res, 0.5);
    g_expected = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY,
                  INFINITY, INFINITY, INFINITY, INFINITY, 6.15,     6.15,
                  6.15,     6.15,     6.15,     6.15,     6.15,     6.15,
                  6.15,     6.15,     6.15};
    for (int i = 0; i < 10; ++i) { // NOTE: INFINITY - INFINITY is NaN
      expect_true(g_expected[i] == g[i]);
    }
    for (int i = 10; i < 21; ++i) {
      expect_true(std::abs(g_expected[i] - g[i]) < 1e-10);
    }
  }

  test_that("Grids formed with zero budget work correctly") {
    MinGrid g(1);
    res.minimize(dps, 0.0);
    g.grid(res, 0.0);
    expect_true(g[0] == 0.0);
    res.minimize(dps, -0.5);
    g.grid(res, 0.0);
    expect_true(g[0] == INFINITY);
    res.minimize(dps, 0.5);
    g.grid(res, 0.0);
    expect_true(g[0] == INFINITY);
  }

  test_that("Grids formed when the mean is changed begin with inf.") {
    MinGrid g(21);
    res.minimize(dps, 1e-6);
    g.grid(res, 0.5);
    expect_true(g[0] == INFINITY);
    res.minimize(dps, -1e-6);
    g.grid(res, 0.5);
    expect_true(g[0] == INFINITY);
  }
}

class MinTreeTest {
public:
  MinTreeTest(const MinTree &t) : t(t){};
  const size_t &get_m() const { return t.m; };
  const double &get_gamma() const { return t.gamma; };
  const size_t &get_i_tau() const { return t.i_tau; };
  const std::vector<double> &get_taus() const { return t.taus; };
  const std::vector<MinRes> &get_ress() const { return t.ress; };
  const std::vector<MinRes> &get_combs() const { return t.combs; };
  const MinGrid &get_grid() const { return t.grid; };
  const std::vector<std::vector<double>> &get_beta_min() const {
    return t.beta_min;
  };
  const std::vector<std::vector<double>> &get_beta_max() const {
    return t.beta_max;
  };

private:
  const MinTree &t;
};

context("Minimization trees (MinTree)") {
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
    MinTree t_(data, 11, 0.1);
    MinTreeTest t(t_);

    expect_true(t.get_m() == 11);
    expect_true(t.get_gamma() == 0.1);
    expect_true(t.get_i_tau() == 1);
    expect_true(t.get_taus().size() == 2);
    expect_true(t.get_ress().size() == 4);
    expect_true(t.get_combs().size() == 3);
    expect_true(t.get_grid().size() == 11);
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
    MinTree t_(data, 11, 0.2);
    t_.minimize();
    MinTreeTest t(t_);
    std::vector<MinResTest> ress;
    for (auto &res : t.get_ress())
      ress.push_back(MinResTest(res));
    std::vector<MinResTest> combs;
    for (auto &res : t.get_combs())
      combs.push_back(MinResTest(res));

    expect_true(t.get_taus().size() == 2);
    std::vector<double> expected_taus = {0.0, 0.0};
    for (int i = 0; i < 2; ++i)
      expect_true(std::abs(t.get_taus()[i] - expected_taus[i]) < 1e-10);
    std::vector<double> expected_epsilon_cum = {0.0, 0.2, 0.6, 1.2, 2.0, 3.0};
    std::vector<double> expected_Sigma_cum = {0.0,   -0.18, -0.46,
                                              -0.76, -1.0,  -1.1};
    std::vector<double> expected_delta_cum = {1.0, 0.8, 0.6, 0.4, 0.2, 0.0};
    std::vector<double> expected_kappa_cum = {0.5, 1.0, 1.5, 2.0, 2.5, 0.0};
    std::vector<int> offset = {11, 9, 7, 5, 3, 1};
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 6; ++i) {
        expect_true(std::abs(ress[j].get_epsilon_cum()[33 - offset[i]] -
                             expected_epsilon_cum[i]) < 1e-10);
        expect_true(std::abs(ress[j].get_Sigma_cum()[33 - offset[i]] -
                             expected_Sigma_cum[i]) < 1e-10);
        expect_true(std::abs(ress[j].get_delta_cum()[33 - offset[i]] -
                             expected_delta_cum[i]) < 1e-10);
        expect_true(std::abs(ress[j].get_kappa_cum()[33 - offset[i]] -
                             expected_kappa_cum[i]) < 1e-10);
      }
    }
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 3; ++j) {
        expect_true(std::abs(combs[j].get_epsilon_cum()[i] -
                             (4 - j) * expected_epsilon_cum[i]) < 1e-10);
        expect_true(std::abs(combs[j].get_Sigma_cum()[i] -
                             (4 - j) * expected_Sigma_cum[i]) < 1e-10);
        expect_true(std::abs(combs[j].get_delta_cum()[i] -
                             expected_delta_cum[i]) < 1e-10);
        expect_true(std::abs(combs[j].get_kappa_cum()[i] -
                             (4 - j) * expected_kappa_cum[i]) < 1e-10);
      }
    }
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
    MinTree t_(data, 10, 0.6);
    t_.minimize();
    t_.regress();
    MinTreeTest t(t_);
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
    t_.minimize();
    t_.regress();

    // Recalculate beta, noting that the first stratum has 0.15 added to it,
    // which increases rho0 by 0.15 * 4 = 0.6, and increases the sum of squares
    // by 0.22
    //              sigma_1     - sigma_0
    double beta1 = (10.0 / 14.0 - 10.0 / 14.0)
                   // sigma_0   * rho_0 - r_0 + sigma_1   * rho_1 - r_1
                   + (10.0 / 14.0 * 5.7 - 4.5 + 10.0 / 14.0 * 6.5 - 5.5)
                         // \sum R^2 - rho_0^2 / n_0 - rho_1^2 / n_1
                         / (7.36 - 5.7 * 5.7 / 14.0 - 6.5 * 6.5 / 14.0)
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
    MinTree t_(data, 10, 0.2);
    MinTreeTest t(t_);
    t_.remean(std::vector<double>{0.0, 0.1});
    expect_true(t.get_taus()[0] == 0.0);
    expect_true(t.get_taus()[1] == 0.1);
    expect_true(t.get_i_tau() == 1);
    t_.remean(std::vector<double>{0.1, 0.1});
    expect_true(t.get_taus()[0] == 0.1);
    expect_true(t.get_taus()[1] == 0.1);
    expect_true(t.get_i_tau() == 0);

    t_.remean(std::vector<double>{1e-6, 1e-6});
    t_.minimize();
    t_.regress();

    expect_true(t.get_beta_min()[0][0] == INFINITY);
    expect_true(t.get_beta_max()[0][0] == -INFINITY);
  }
}
