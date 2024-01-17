#include "regdata.h"

//////////////////////////// RegData Implementation ////////////////////////////

/*
 * @brief Private class helper function for converting estimates and upper and
 *        lower bounds into PtVectors.
 * @param es Risk vectors.
 * @param lwrs Lower bounds.
 * @param uprs Upper bounds.
 * @pre The i-th element of `es` is the risk vector for the i-th stratum, and
 *      has the same length as the i-th element of `lwrs` and `uprs`. These
 *      lengths are all greater than or equal to 1.
 * @return A vector of PtVectors generated from the input values.
 */
std::vector<PtVector>
RegData::compute_dps(const std::vector<std::vector<double>> &es,
                     const std::vector<std::vector<double>> &lwrs,
                     const std::vector<std::vector<double>> &uprs) const {
  // Initialize a temporary PtVector vector
  std::vector<PtVector> dps_tmp(es.size());
  for (size_t i = 0; i < es.size(); ++i)
    dps_tmp[i] = PtVector(es[i], lwrs[i], uprs[i]);
  return dps_tmp;
}

/*
 * @brief Private class helper function for converting estimates into group
 *        sizes.
 * @param es Risk vectors.
 */
std::vector<size_t>
RegData::compute_n(const std::vector<std::vector<double>> &es) const {
  std::vector<size_t> n_tmp(g);
  for (size_t j = 0; j < g; ++j)
    n_tmp[j] = es[j].size() + es[j + g].size();
  return n_tmp;
}

/*
 * @brief Private class helper function for converting estimates into group
 *        search rates.
 * @param es Risk vectors.
 */
std::vector<double>
RegData::compute_sigma(const std::vector<std::vector<double>> &es) const {
  std::vector<double> sigma_tmp(g);
  for (size_t j = 0; j < g; ++j)
    sigma_tmp[j] = (double)es[j + g].size() / (double)n[j];
  return sigma_tmp;
}

/*
 * @brief Private class helper function for converting estimates into group
 *        average risk.
 * @param es Risk vectors.
 */
std::vector<double>
RegData::compute_rho(const std::vector<std::vector<double>> &es) const {
  std::vector<double> rho_tmp(g);
  for (size_t j = 0; j < g; ++j)
    rho_tmp[j] = (std::accumulate(es[j].begin(), es[j].end(), 0.0) +
                  std::accumulate(es[j + g].begin(), es[j + g].end(), 0.0)) /
                 n[j];
  return rho_tmp;
}

/*
 * @brief Private class helper function for converting estimates into
 *        numerator for regression coefficients.
 * @param es Risk vectors.
 */
double
RegData::compute_numerator(const std::vector<std::vector<double>> &es) const {
  double numerator_tmp = 0.0;
  for (size_t j = 0; j < g; ++j) {
    numerator_tmp += sigma[j] * rho[j] * n[j] -
                     std::accumulate(es[j + g].begin(), es[j + g].end(), 0.0);
  }
  return numerator_tmp;
}

/*
 * @brief Private class helper function for converting estimates into
 *        sum of squares of risk.
 * @param es Risk vectors.
 */
double RegData::compute_ss(const std::vector<std::vector<double>> &es) const {
  double ss_tmp = 0.0;
  for (size_t j = 0; j < es.size(); ++j)
    for (size_t i = 0; i < es[j].size(); ++i)
      ss_tmp += std::pow(es[j][i], 2);
  return ss_tmp;
}

/*
 * @brief Constructor for RegData (regression data) class.
 * @param es Risk vectors
 * @param lwrs Lower bounds
 * @param uprs Upper bounds
 * @pre `es`, `lwrs`, and `uprs` have the same length, which is even.
 * @pre The i-th element of `es` is the risk vector for the i-th stratum, and
 *      has the same length as the i-th element of `lwrs` and `uprs`. These
 *      lengths are all greater than or equal to 1. The unobserved strata come
 *      first, with each corresponding observed stratum occurring es.size() / 2
 *      indices later.
 * @pre All of the unobserved strata occur before all of the observed strata.
 * @pre Each element of `es`, `lwrs`, and `uprs` is sorted in ascending order.
 * @pre Each element of `lwrs` is less than or equal to the corresponding risk
 *      vector element, which is less than or equal to the corresponding element
 *      of `uprs`.
 * @pre Each element of `es`, `lwrs`, and `uprs` contains probabilities.
 * @pre There are at least two groups, i.e., `es.size()` is greater than or
 *      equal to four.
 */
RegData::RegData(const std::vector<std::vector<double>> &es,
                 const std::vector<std::vector<double>> &lwrs,
                 const std::vector<std::vector<double>> &uprs)
    : g(es.size() / 2), es(std::move(es)), lwrs(std::move(lwrs)),
      uprs(std::move(uprs)), dps(compute_dps(es, lwrs, uprs)), n(compute_n(es)),
      sigma(compute_sigma(es)), rho(compute_rho(es)),
      numerator(compute_numerator(es)), ss(compute_ss(es)) {}
