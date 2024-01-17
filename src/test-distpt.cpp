#include "distpt.h"
#include <testthat.h>

context("Distinguished Point (DistPt)") {
  test_that("DistPt initializer list constructor works") {
    DistPt pt = {0.5, PtType::est};
    expect_true(pt.val == 0.5);
    expect_true(pt.type == PtType::est);
  }
}

context("Distinguished Point Vector (PtVector)") {
  auto e = std::vector<double>{1, 0.75, 0.5, 0.25, 0};
  auto l = std::vector<double>(5, 0);
  auto u = std::vector<double>(5, 1);

  PtVector pv(e, l, u);

  std::vector<double> vals = {0,    0, 0, 0, 0, 0, 0.25, 0.5,
                              0.75, 1, 1, 1, 1, 1, 1};
  std::vector<PtType> typ = {PtType::lwr, PtType::lwr, PtType::lwr, PtType::lwr,
                             PtType::lwr, PtType::est, PtType::est, PtType::est,
                             PtType::est, PtType::est, PtType::upr, PtType::upr,
                             PtType::upr, PtType::upr, PtType::upr};

  test_that("DistPt values are correct") {
    for (std::size_t i = 0; i < 3 * e.size(); ++i) {
      expect_true(pv[i].val == vals[i]);
    }
  }
  test_that("DistPt types are correct") {
    for (std::size_t i = 0; i < 3 * e.size(); ++i) {
      expect_true(pv[i].type == typ[i]);
    }
  }
}
