#include "iter.h"
#include <testthat.h>

/*
 * @brief Test fixture for the $\ell^1$ iterator.
 */
class EllOneIterTest {
public:
  EllOneIterTest(const EllOneIter &iter) : iter(iter) {}

  const bool get_valid() const { return iter.valid; }
  const size_t get_g() const { return iter.g; }
  const double get_delta() const { return iter.delta; }
  const std::vector<double> &get_d() const { return iter.d; }
  const std::vector<bool> &get_dir() const { return iter.dir; }
  const std::vector<double> &get_taus() const { return iter.taus; }
  const double &get_eta() const { return iter.eta; }
  const std::vector<double> &get_etas() const { return iter.etas; }
  const std::vector<double> &get_lwr() const { return iter.lwr; }
  const std::vector<double> &get_upr() const { return iter.upr; }
  const double get_epsilon() const { return iter.epsilon; }

private:
  const EllOneIter &iter;
};

context("$\\ell^1$ Iterator (EllOneIter)") {
  std::vector<double> etas = {0.1, 0.2};
  std::vector<double> lwr = {-0.99, -0.99};
  std::vector<double> upr = {0.99, 1.99};
  double epsilon = 0.274;
  EllOneIter iter_(2, 0.025, etas, lwr, upr, epsilon);
  EllOneIterTest iter(iter_);

  test_that("Initialization is correct") {
    expect_true(iter.get_valid());
    expect_true(iter.get_g() == 2);
    expect_true(std::abs(iter.get_delta()) < 1e-10);
    expect_true(std::abs(iter.get_epsilon() - epsilon) < 1e-10);
    expect_true(std::abs(iter.get_eta() - 0.025) < 1e-10);
    for (size_t i = 0; i < 2; ++i) {
      expect_true(iter.get_d()[i] == 0.0);
      expect_true(iter.get_dir()[i]);
      expect_true(iter.get_taus()[i] == 0.0);
      expect_true(std::abs(iter.get_etas()[i] - etas[i]) < 1e-10);
      expect_true(std::abs(iter.get_lwr()[i] - lwr[i]) < 1e-10);
      expect_true(std::abs(iter.get_upr()[i] - upr[i]) < 1e-10);
    }
  }

  test_that("Iteration is correct") {
    iter_.next();
    expect_true(iter.get_valid());
    expect_true(std::abs(iter.get_taus()[0] - 0.1) < 1e-10);
    expect_true(std::abs(iter.get_d()[0] - 0.025) < 1e-10);

    // Iterate nineteen additional times
    for (size_t i = 0; i < 19; ++i) {
      iter_.next();
    }

    expect_true(iter.get_valid());
    expect_true(std::abs(iter.get_taus()[0] - 0.1) < 1e-10);
    expect_true(std::abs(iter.get_taus()[1] - 0.2) < 1e-10);
    expect_true(std::abs(iter.get_d()[0] - 0.025) < 1e-10);
    expect_true(std::abs(iter.get_d()[1] - 0.025) < 1e-10);
    expect_true(std::abs(iter.get_delta() - 0.05) < 1e-10);
    expect_true(iter.get_dir()[0]);
    expect_true(iter.get_dir()[1]);

    // Iterate nine additional times
    for (size_t i = 0; i < 9; ++i) {
      iter_.next();
    }

    expect_true(iter.get_valid());
    expect_true(std::abs(iter.get_taus()[0] + 0.1) < 1e-10);
    expect_true(std::abs(iter.get_taus()[1] - 0.2) < 1e-10);
    expect_true(std::abs(iter.get_d()[0] - 0.025) < 1e-10);
    expect_true(std::abs(iter.get_d()[1] - 0.025) < 1e-10);
    expect_true(std::abs(iter.get_delta() - 0.05) < 1e-10);
    expect_true(!iter.get_dir()[0]);
    expect_true(iter.get_dir()[1]);

    // Iterate 152 additional times
    for (size_t i = 0; i < 152; ++i) {
      iter_.next();
    }

    expect_true(iter.get_valid());
    expect_true(std::abs(iter.get_taus()[0] + 0.6) < 1e-10);
    expect_true(std::abs(iter.get_taus()[1] + 0.8) < 1e-10);
    expect_true(std::abs(iter.get_d()[0] - 0.15) < 1e-10);
    expect_true(std::abs(iter.get_d()[1] - 0.1) < 1e-10);
    expect_true(std::abs(iter.get_delta() - 0.25) < 1e-10);
    expect_true(!iter.get_dir()[0]);
    expect_true(!iter.get_dir()[1]);

    // Iterate 1 additional time
    // This should cause the iterator to become invalid
    iter_.next();

    expect_true(!iter.get_valid());
  }
}

class TauChunkTest {
public:
  TauChunkTest(const TauChunk &chunk) : chunk(chunk) {}

  const size_t get_size() const { return chunk.size; }
  const EllOneIter &get_iter() const { return chunk.iter; }
  const std::vector<std::vector<double>> &get_taus() const {
    return chunk.taus;
  }

private:
  const TauChunk &chunk;
};

context("Iterator chunks (TauChunk)") {
  std::vector<double> etas = {0.1, 0.2};
  std::vector<double> lwr = {-0.99, -0.99};
  std::vector<double> upr = {0.99, 1.99};
  double epsilon = 0.274;
  EllOneIter iter_(2, 0.025, etas, lwr, upr, epsilon);
  std::mutex mtx;
  TauChunk chunk_(10, mtx, iter_);
  TauChunkTest chunk(chunk_);

  test_that("Initialization is correct") {
    expect_true(chunk.get_size() == 10);
    expect_true(chunk.get_iter().g == 2);
    expect_true(chunk.get_iter().valid);
    expect_true(chunk.get_taus().size() == 10);
    for (size_t i = 0; i < 10; ++i) {
      expect_true(chunk.get_taus()[i].size() == 2);
      expect_true(std::abs(chunk.get_taus()[i][0] - 0.1 * i) < 1e-10);
      expect_true(chunk.get_taus()[i][1] == 0.0);
    }
  }

  test_that("Chunk will not grab taus after invalidity") {
    // Refill seventeen times
    for (size_t i = 0; i < 17; ++i) {
      chunk_.refill();
      expect_true(chunk.get_taus().size() == 10);
    }
    chunk_.refill();
    expect_true(chunk.get_taus().size() == 2);
  }

  test_that("Chunk loops through taus correctly.") {
    size_t i = 0;
    for (auto tau : chunk_) {
      expect_true(std::abs(tau[0] - 0.1 * i) < 1e-10);
      expect_true(tau[1] == 0.0);
      ++i;
    }
  }

  test_that("Creating a vector of chunks works correctly.") {
    EllOneIter iter1_(2, 0.025, etas, lwr, upr, epsilon);
    EllOneIter iter2_(2, 0.025, etas, lwr, upr, epsilon);
    std::mutex mtx1;
    std::mutex mtx2;
    TauChunk base_chunk_(1000, mtx1, iter1_);
    TauChunkTest base_chunk(base_chunk_);
    std::vector<TauChunk> chunks;
    for (int i = 0; i < 11; ++i)
      chunks.push_back(TauChunk(10, mtx2, iter2_));
    int i = 0;
    for (auto chunk : chunks) {
      for (auto tau : chunk) {
        expect_true(tau[0] == base_chunk.get_taus()[i][0]);
        expect_true(tau[1] == base_chunk.get_taus()[i][1]);
        ++i;
      }
    }
  }
}
