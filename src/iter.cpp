#include "iter.h"

///////////////////////////// $\ell^1$ Iterator ////////////////////////////////

/*
 * @brief Constructor for the $\ell^1$ iterator.
 * @param g       The number of groups.
 * @param lwr     The lower bound.
 * @param upr     The upper bound.
 * @param eta     The step size cost, in the same units as epsilon.
 * @param eta     The step sizes.
 * @param epsilon The $\ell^1$ bound.
 */
EllOneIter::EllOneIter(const size_t g, const double eta,
                       const std::vector<double> &etas,
                       const std::vector<double> &lwr,
                       const std::vector<double> &upr, const double epsilon)
    : valid(true), g(g), delta(0), d(g, 0.0), dir(g, true), taus(g, 0.0),
      eta(eta), etas(etas), lwr(lwr), upr(upr), epsilon(epsilon) {}

/*
 * @brief Private method for advancing the iterator.
 * @param j The index to advance.
 */
void EllOneIter::advance(size_t j) {
  // If increasing...
  if (dir[j]) {
    // Update mean
    taus[j] += etas[j];
    // Update $\ell^1$ norm
    delta += eta;
    d[j] += eta;
    // Check that we're in the feasible region
    if (taus[j] >= upr[j] || delta > epsilon) {
      // If not, go back to zero
      taus[j] = 0.0;
      dir[j] = false;
      delta -= d[j];
      d[j] = 0.0;
    }
  }

  // If decreasing...
  if (!dir[j]) {
    // Update mean
    taus[j] -= etas[j];
    // Update $\ell^1$ norm
    delta += eta;
    d[j] += eta;
    // Check that we're in the feasible region
    if (taus[j] <= lwr[j] || delta > epsilon) {
      // If not, go back to zero...
      taus[j] = 0.0;
      dir[j] = true;
      delta -= d[j];
      d[j] = 0.0;
      // ...and advance the next index
      if (j < g - 1) {
        advance(j + 1);
        // unless there are no more indices to advance
      } else {
        valid = false;
      }
    }
  }
}

/*
 * @brief Generate the next value of the iterator.
 * @pre The iterator is valid.
 */
void EllOneIter::next() { advance(0); }

/*
 * @brief Return the current value of the iterator.
 * @pre The iterator is valid.
 */
const std::vector<double> &EllOneIter::get_taus() const { return taus; }

///////////////////////////////// Tau Chunk ////////////////////////////////////

/*
 * @brief Constructor for the tau chunk, a container for a chunk of mean
 *        perturbations to chunk work for a thread.
 * @param size The (maximum) size of the chunk.
 * @param iter The iterator.
 */
TauChunk::TauChunk(const size_t size, std::mutex &mtx, EllOneIter &iter)
    : valid(iter.valid), size(size), iter(iter), mtx(mtx), taus(size) {
  // Fill the chunk
  refill();
}

/*
 * @brief Iterator begin method.
 * @return An iterator to the beginning of the chunk.
 */
TauChunk::iterator TauChunk::begin() { return taus.begin(); }

/*
 * @brief Iterator end method.
 * @return An iterator to the end of the chunk.
 */
TauChunk::iterator TauChunk::end() { return taus.end(); }

/*
 * @brief Const iterator begin method.
 * @return A const iterator to the beginning of the chunk.
 */
TauChunk::const_iterator TauChunk::begin() const { return taus.begin(); }

/*
 * @brief Const iterator end method.
 * @return A const iterator to the end of the chunk.
 */
TauChunk::const_iterator TauChunk::end() const { return taus.end(); }

/*
 * @brief Const iterator begin method.
 * @return A const iterator to the beginning of the chunk.
 */
TauChunk::const_iterator TauChunk::cbegin() const { return taus.cbegin(); }

/*
 * @brief Const iterator end method.
 * @return A const iterator to the end of the chunk.
 */
TauChunk::const_iterator TauChunk::cend() const { return taus.cend(); }

/*
 * @brief Indexing operator.
 * @param i The index.
 * @return A reference to the tau at index `i`.
 */
std::vector<double> &TauChunk::operator[](const size_t i) { return taus[i]; }

/*
 * @brief Const indexing operator.
 * @param i The index.
 * @return A const reference to the tau at index `i`.
 */
const std::vector<double> &TauChunk::operator[](const size_t i) const {
  return taus[i];
}

/*
 * @brief Refill the chunk, by drawing at most `size` samples from the iterator.
 */
size_t TauChunk::refill() {
  // Lock the mutex
  std::lock_guard<std::mutex> lock(mtx);

  // Fill the chunk
  for (size_t i = 0; i < size; ++i) {

    // If the iterator is valid...
    if (valid) {
      // Copy the current taus
      taus[i] = iter.get_taus();
      // Advance the iterator.
      iter.next();
      // If the iterator is invalid...
    } else {
      // ... shrink the chunk to the number of samples drawn and break.
      taus.resize(i);
      break;
    }
  }
  return taus.size();
}
