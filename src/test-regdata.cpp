#include "regdata.h"
#include <testthat.h>

context("Regression Data (RegData)"){test_that(
    "Values and dimensions of members are correct"){
    // First stratum has three people, second has two, third has one, fourth has
    // seven.
    auto es = std::vector<std::vector<double>>{
        {0.0, 0.1, 0.2},
        {0.3, 0.4},
        {0.5},
        {0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9}};
auto lwrs = std::vector<std::vector<double>>{
    {0.0, 0.0, 0.0}, {0.0, 0.0}, {0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
auto uprs = std::vector<std::vector<double>>{
    {1.0, 1.0, 1.0}, {1.0, 1.0}, {1.0}, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};
RegData data(es, lwrs, uprs);

expect_true(data.es.size() == 4);
expect_true(data.es[0].size() == 3);
expect_true(data.es[1].size() == 2);
expect_true(data.es[2].size() == 1);
expect_true(data.es[3].size() == 7);
expect_true(data.lwrs.size() == 4);
expect_true(data.lwrs[0].size() == 3);
expect_true(data.lwrs[1].size() == 2);
expect_true(data.lwrs[2].size() == 1);
expect_true(data.lwrs[3].size() == 7);
expect_true(data.uprs.size() == 4);
expect_true(data.uprs[0].size() == 3);
expect_true(data.uprs[1].size() == 2);
expect_true(data.uprs[2].size() == 1);
expect_true(data.uprs[3].size() == 7);
expect_true(data.dps[0].size() == 9);
expect_true(data.dps[1].size() == 6);
expect_true(data.dps[2].size() == 3);
expect_true(data.dps[3].size() == 21);
expect_true(data.n.size() == 2);
expect_true(data.n[0] == 4);
expect_true(data.n[1] == 9);
expect_true(data.sigma.size() == 2);
expect_true(std::abs(data.sigma[0] - 0.25) < 1e-10);
expect_true(std::abs(data.sigma[1] - 7.0 / 9.0) < 1e-10);
expect_true(data.rho.size() == 2);
expect_true(std::abs(data.rho[0] - 0.2) < 1e-10);
expect_true(std::abs(data.rho[1] - 5.95 / 9) < 1e-10);
expect_true(std::abs(data.numerator - (0.25 * 0.5 - 0.75 * 0.1 +
                                       7.0 / 9.0 * 0.75 - 2.0 / 9.0 * 0.35)));
expect_true(std::abs(data.ss - 4.5575) < 1e-10);
}
}
;
