#pragma once
#include <mutex>
#include <numeric>
#include <vector>

class EllOneIter {
  friend class EllOneIterTest; // NOTE: For unit tests only

public:
  bool valid;     // Whether the iterator is valid
  const size_t g; // Number of groups

  EllOneIter(const size_t g, const double eta, const std::vector<double> &etas,
             const std::vector<double> &lwr, const std::vector<double> &upr,
             const double epsilon);

  void next();
  const std::vector<double> &get_taus() const;

private:
  double delta;             // Current $\ell^1$ norm
  std::vector<double> d;    // Contribution to $\ell^1$ norm at index
  std::vector<bool> dir;    // Direction of the step
  std::vector<double> taus; // Current state
  const double eta;         // Step size cost, in the same units as epsilon
  const std::vector<double> etas; // Step sizes
  const std::vector<double> lwr;  // Lower bounds
  const std::vector<double> upr;  // Upper bounds
  const double epsilon;           // $\ell^1$ bound

  std::vector<double>
  compute_p(const std::vector<std::vector<double>> &es) const;
  void advance(size_t i);
};

class TauChunk {
  friend class TauChunkTest; // NOTE: For unit tests only

public:
  bool &valid;

  TauChunk(const size_t size, std::mutex &mtx, EllOneIter &iter);

  typedef std::vector<std::vector<double>>::iterator iterator;
  typedef std::vector<std::vector<double>>::const_iterator const_iterator;

  std::vector<std::vector<double>>::iterator begin();
  std::vector<std::vector<double>>::iterator end();
  std::vector<std::vector<double>>::const_iterator begin() const;
  std::vector<std::vector<double>>::const_iterator end() const;
  std::vector<std::vector<double>>::const_iterator cbegin() const;
  std::vector<std::vector<double>>::const_iterator cend() const;

  std::vector<double> &operator[](const size_t i);
  const std::vector<double> &operator[](const size_t i) const;

  size_t refill();

private:
  const size_t size;
  EllOneIter &iter;
  std::mutex &mtx;
  std::vector<std::vector<double>> taus;
};
