#include "max.h"

/////////////////////////// MaxRes Implementation //////////////////////////////

/*
 * @brief Constructor for MaxRes object.
 * @param n Length of the vectors.
 * @pre `n` is greater than or equal to 1.
 * @return A MaxRes object with n-length vectors allocated.
 */
MaxRes::MaxRes(size_t n)
    : n(n), capacity(n), epsilon_cum(n), Sigma_cum(n), delta_cum(n){};

/*
 * @brief Resize the vectors in a MaxRes object.
 * @param n New length of the vectors.
 * @pre `n` is greater than or equal to 1.
 */
void MaxRes::resize(size_t n) {
  this->n = n;
  if (n > capacity) {
    capacity = n;
    epsilon_cum.resize(n);
    Sigma_cum.resize(n);
    delta_cum.resize(n);
  }
}

/*
 * @brief Maximize the sum of squares of a risk vector.
 * @param e Risk vector.
 * @param lwr Lower bound vector.
 * @param upr Upper bound vector.
 * @param tau Average change in risk.
 * @pre `e`, `upr`, and `lwr` have the same length.
 * @pre `e`, `upr`, and `lwr` are sorted in ascending order.
 * @pre `lwr[i]` <= `e[i]` <= `upr[1]`.
 * @pre `e`, `upr`, and `lwr` all contain probabilities.
 * @pre The length is greater than or equal to 1.
 * @return A MaxRes object representing all solutions with mean tau.
 */
void MaxRes::maximize(const std::vector<double> &e,
                      const std::vector<double> &lwr,
                      const std::vector<double> &upr, const double tau) {
  // Store vector length
  size_t n = e.size();

  // Initialize "active" values for remeaning and maximization
  size_t i_lwr = 0;        // Lower "active" index
  size_t i_upr = n - 1;    // Upper "active" index
  double a_lwr = e[0];     // Lower "active" value
  double a_upr = e[n - 1]; // Upper "active" value

  // REMEAN
  // Find the minimum L_1 distance perturbation having the target mean
  double D = n * tau; // Difference between required and actual mean
  double epsilon_min = std::abs(D); // Minimum budget expenditure to "fix" mean
  double Sigma = 0;                 // Change in sum of squares

  if (D < 0) { // If mean is to be decreased...
    D = std::abs(D);
    double d = e[i_lwr] - lwr[i_lwr];
    while (D > d && i_lwr < n - 1) { // ... decrese to lwr on lower end.
      D -= d;
      Sigma -= d * (e[i_lwr] + lwr[i_lwr]);
      ++i_lwr;
      d = e[i_lwr] - lwr[i_lwr];
    }
    Sigma +=
        D *
        (D -
         2 * e[i_lwr]); // Final reduction does not go all the way down to lwr.
    a_lwr = e[i_lwr] - D;
  } else if (D > 0) { // ... increase to upr on the upper end.
    double d = upr[i_upr] - e[i_upr];
    while (D > d && i_upr > 0) {
      D -= d;
      Sigma += d * (e[i_upr] + upr[i_upr]);
      --i_upr;
      d = upr[i_upr] - e[i_upr];
    }
    Sigma +=
        D *
        (D + 2 * e[i_upr]); // Final increase does not go all the way up to upr.
    a_upr = e[i_upr] + D;
  } else {
    epsilon_min = 0; // No change in mean, so no budget is spent.
    Sigma = 0;
  }

  // MAXIMIZE
  // Find the maximum sum of squares (SS) among all solutions within given L_1
  // distance
  n = i_upr - i_lwr + 1;
  resize(n);
  double delta = a_upr - a_lwr;
  double Delta_lwr = a_lwr - lwr[i_lwr];
  double Delta_upr = upr[i_upr] - a_upr;
  double Delta = std::min(Delta_lwr, Delta_upr);

  // Store results after zero iterations
  epsilon_cum[0] = epsilon_min;
  Sigma_cum[0] = Sigma;
  delta_cum[0] = delta;

  // Loop over remaining iterations
  for (size_t i = 1; i < n; ++i) {
    // Advance one iteration
    if (Delta_lwr == Delta) {
      ++i_lwr;                        // Increment lower "active" index
      a_lwr = e[i_lwr];               // Calcualte new lower "active" value
      Delta_lwr = a_lwr - lwr[i_lwr]; // Calculate new lower "active" gap
      a_upr = a_upr + Delta;          // Calculate new upper "active" value
      Delta_upr = Delta_upr - Delta;  // Calculate new upper "active" gap
    } else {
      --i_upr;                        // Decrement upper "active" index
      a_upr = e[i_upr];               // Calculate new upper "active" value
      Delta_upr = upr[i_upr] - a_upr; // Calculate new upper "active" gap
      a_lwr = a_lwr - Delta;          // Calculate new lower "active" value
      Delta_lwr = Delta_lwr - Delta;  // Calculate new lower "active" gap
    }

    epsilon_cum[i] = 2 * Delta; // Budget spent in this iteration
    Sigma_cum[i] = 2 * Delta * (delta + Delta); // Increase in sum of squares
    Delta = std::min(Delta_upr,
                     Delta_lwr); // Next gap is minimum of upper and lower gaps
    delta =
        a_upr - a_lwr; // Next difference is upper minus lower "active" values
    delta_cum[i] = delta; // Store next difference
  }
  delta_cum[n - 1] = 0.0; // Last difference is always zero

  // Calculate cumulative values
  for (size_t i = 1; i < n; ++i) {
    epsilon_cum[i] += epsilon_cum[i - 1];
    Sigma_cum[i] += Sigma_cum[i - 1];
  }
}

/////////////////////////// MaxGrid Implementation /////////////////////////////

/*
 * @brief Constructor for MaxGrid object.
 * @param m Length of the grid.
 * @pre `m` is greater than or equal to 1.
 * @return A MaxGrid object with length `m`.
 */
MaxGrid::MaxGrid(size_t m) : m(m), g(m){};

/*
 * @brief Constructor for MaxGrid object.
 * @param init_list Initializer list for the grid.
 * @pre The length of `init_list` is greater than or equal to 1.
 * @return A MaxGrid object with length `init_list.size()`.
 */
MaxGrid::MaxGrid(std::initializer_list<double> init_list)
    : m(init_list.size()), g(init_list){};

/*
 * @brief Size of the grid.
 * @return The length of the grid.
 */
size_t MaxGrid::size() const { return m; }

/*
 * @brief Subscript operator for MaxGrid object.
 * @param index Index of the grid point.
 * @pre `index` is greater than or equal to 0.
 * @pre `index` is less than or equal to the length of the grid.
 * @return The value of the grid at the given index.
 */
double &MaxGrid::operator[](std::size_t index) { return g[index]; }

/*
 * @brief Subscript operator for MaxGrid object.
 * @param index Index of the grid point.
 * @pre `index` is greater than or equal to 0.
 * @pre `index` is less than or equal to the length of the grid.
 * @return The value of the grid at the given index.
 */
const double &MaxGrid::operator[](std::size_t index) const { return g[index]; }

/*
 * @brief Evaluate a MaxGrid object on a grid of its points.
 * @param res MaxRes object.
 * @param m Length of the grid.
 * @param gamma Step size.
 * @pre `epsilon` is greater than or equal to 0.
 * @pre m is greater than or equal to 1.
 * @pre res.n has been maximized.
 * @return A MaxGrid object representing the maximum sum of squares on a grid.
 */
void MaxGrid::grid(const MaxRes &res, const double gamma) {
  size_t i = 0;     // Pointer to position in res
  size_t n = res.n; // Length of result of maximization routine
  double epsilon_max = res.epsilon_cum[n - 1]; // Maximum budget expendable
  double epsilon_cur = res.epsilon_cum[i]; // Budget expended at "current" step
  double epsilon_itr = 0.0; // Budget exhausted up to current iteration
  // Budget expended at "next" step, difference between active values, change
  // in sum of squares, change in difference between active values
  double epsilon_nxt, delta, Sigma, Delta;

  // Early return if the MaxRes is length one
  if (n == 1) {
    Sigma = res.Sigma_cum[n - 1];
    for (size_t k = 0; k < m; ++k) {
      g[k] = (epsilon_itr <= epsilon_max)
                 ? -std::numeric_limits<double>::infinity()
                 : Sigma;
      epsilon_itr += gamma;
    }
    return;
  }

  epsilon_nxt = res.epsilon_cum[i + 1]; // Budget expended at "next" step
  delta = res.delta_cum[i];             // Difference between active values
  // NOTE: If the current point is feasible, then Sigma is highest value;
  // otherwise it takes the sentinel value of -Inf. Due to floating point
  // errors, equality is checked with a tolerance.
  Sigma = (epsilon_itr >= epsilon_cur - 1e-10)
              ? res.Sigma_cum[i]
              : -std::numeric_limits<double>::infinity();

  for (size_t k = 0; k < m; ++k) { // Loop through grid points
    while (epsilon_itr >= epsilon_nxt && i < n - 1) { // Advance to next point
      ++i;
      epsilon_cur = res.epsilon_cum[i];
      delta = res.delta_cum[i];
      Sigma = res.Sigma_cum[i];
      if (i < n - 1) {
        epsilon_nxt = res.epsilon_cum[i + 1];
      } else {
        epsilon_nxt = epsilon_max;
      }
    }
    // NOTE: Due to floating point errors, equality is checked with a tolerance.
    if (epsilon_itr >= epsilon_cur - 1e-10) {
      if (epsilon_itr < epsilon_max) {
        // Change in difference between active values
        Delta = epsilon_itr - epsilon_cur;
        // Increase in sum of squares from extra budget
        g[k] = Sigma + Delta * (delta + Delta / 2); // Increase
      } else {
        g[k] = Sigma; // epsilon_max exceeded, sum of squares cannot increase
      }
    } else {
      g[k] = -std::numeric_limits<double>::infinity(); // Infeasible
    }
    epsilon_itr += gamma;
  }
}

/*
 * @brief Combine values of grids of points to determine the maximum sum of
 *       squares across multiple groups.
 * @param g1 A grid of sum of square values.
 * @param g2 A grid of sum of square values.j
 * @pre The lengths of g1 and g2 match the calling MaxGrid object.
 * @return The (approximate) maximum sum of squares for each given budget across
 *         both strata.
 */
void MaxGrid::combine(const MaxGrid &g1, const MaxGrid &g2) {
  for (size_t i = 0; i < m; ++i) {
    double M = -std::numeric_limits<double>::infinity(); // Initialized to
                                                         // negative infinity
    for (size_t j = 0; j <= i; ++j) {
      size_t k = i - j;
      M = std::max(M, g1[j] + g2[k]);
    }
    g[i] = M;
  }
}

//////////////////////////// MaxTree Implementation ////////////////////////////

/*
 * @brief Constructor for MaxTree objects.
 * @param data Regression data.
 * @param m Length of the grids.
 * @param gamma Step size.
 * @pre `m` is greater than or equal to 1.
 * @pre `gamma` is greater than 0.
 * @return A MaxTree object with `ns.size()` MaxRes objects, with the i-th
 *         object initiliazed to length `ns[i]`, `ns.size()` MaxGrid objects to
 *         hold the results of the maximization routine, along with their
 *         corresponding grids; and `ns.size()` MaxGrid objects to hold the
 *         combined results. (The final combined MaxGrid object is a dummy, and
 *         points to the final result MaxGrid object.) The maximization routine
 *         is run on the observed strata, and the results are combined
 *         appropriately.
 */
MaxTree::MaxTree(const RegData &data, const size_t m, const double gamma)
    : data(data), g(data.g), m(m), gamma(gamma), i_tau(g - 1), taus(g, 0),
      beta_min(g - 1,
               std::vector<double>(m, std::numeric_limits<double>::infinity())),
      beta_max(g - 1, std::vector<double>(
                          m, -std::numeric_limits<double>::infinity())) {
  // Initialize MaxRes objects
  ress.reserve(2 * g);
  for (size_t i = 0; i < 2 * g; ++i) {
    ress.emplace_back(data.es[i].size());
  }

  // Initialize the result grids
  res_grids.reserve(2 * g);
  for (size_t i = 0; i < 2 * g; ++i) {
    res_grids.emplace_back(m);
  }

  // Initialize the combining grids
  comb_grids.reserve(2 * g - 1);
  for (size_t i = 0; i < 2 * g - 1; ++i) {
    comb_grids.emplace_back(m);
  }

  // Maximize the sum of squares in the observed strata
  for (size_t i = i_tau; i < 2 * g; ++i) {
    ress[i].maximize(data.es[i], data.lwrs[i], data.uprs[i], 0.0);
    res_grids[i].grid(ress[i], gamma);
  }

  // Combine the results of the maximization routine
  // NOTE: The final combined MaxGrid object is combined from the final *two*
  // result MaxGrid objects.
  comb_grids[2 * g - 2].combine(res_grids[2 * g - 2], res_grids[2 * g - 1]);
  // NOTE: Since size_t is unsigned, the loop must be written in this way.
  for (size_t i = 2 * g - 2; i > i_tau; --i)
    comb_grids[i - 1].combine(res_grids[i - 1], comb_grids[i]);
}

/*
 * @brief Accessor function for the minimum regression coefficients.
 * @return The minimum regression coefficients observed so far.
 */
const std::vector<std::vector<double>> &MaxTree::get_beta_min() const {
  return beta_min;
}

/*
 * @brief Accessor function for the maximum regression coefficients.
 * @return The maximum regression coefficients observed so far.
 */
const std::vector<std::vector<double>> &MaxTree::get_beta_max() const {
  return beta_max;
}

/*
 * @brief Remean the risk vectors.
 * @param taus Mean changes in risk.
 * @pre   The length of `taus` is `g`.
 * @return The calling object with the risk vectors remeaned, and the pointer
 *         i_tau updated to the largest index that changed in the remeaned
 *         vectors.
 */
void MaxTree::remean(const std::vector<double> &taus) {
  for (size_t i = 0; i < g; ++i) {
    if (this->taus[i] != taus[i]) {
      this->taus[i] = taus[i];
      i_tau = i;
    }
  }
}

/*
 * @brief Maximize the sum of squares of a risk vector.
 * @return The (approximate) maximum sum of squares of risk vectors across
 *         strata, subject to the given budget, upper and lower bounds, and mean
 *         changes in risk.
 */
void MaxTree::maximize() {
  // Maximize the sum of squares in the strata that need to be updated.
  for (size_t i = 0; i <= i_tau; ++i) {
    ress[i].maximize(data.es[i], data.lwrs[i], data.uprs[i], taus[i]);
    res_grids[i].grid(ress[i], gamma);
  }

  // Combine the results of the maximization routine.
  // NOTE: Since size_t is unsigned, the loop must be written in this way.
  for (size_t i = i_tau + 1; i > 0; --i) {
    comb_grids[i - 1].combine(res_grids[i - 1], comb_grids[i]);
  }
}

/*
 * @brief Calculate the regression coefficients and compare to existing records.
 * @return The largest and smallest regression coefficients observed so far,
 *         stored in `beta_min` and `beta_max` as appropriate.
 */
void MaxTree::regress() {
  data.regress(comb_grids[0], taus, beta_min, beta_max);
}
