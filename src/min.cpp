#include "min.h"

/////////////////////////// MinRes Implementation //////////////////////////////

/*
 * @brief Constructor for MinRes.
 * @param n Length of the vectors to be allocated.
 * @pre `n` is greater than or equal to 1.
 * @return A MinRes object with n-length vectors allocated.
 */
MinRes::MinRes(size_t n)
    : n(n), capacity(n), epsilon_cum(n), Sigma_cum(n), delta_cum(n),
      kappa_cum(n){};

/*
 * @brief Resize the vectors in MinRes object.
 * @param n New length of the vectors.
 * @pre `n` is greater than or equal to 1.
 */
void MinRes::resize(size_t n) {
  this->n = n;
  if (n > capacity) {
    capacity = n;
    epsilon_cum.resize(n);
    Sigma_cum.resize(n);
    delta_cum.resize(n);
    kappa_cum.resize(n);
  }
}

/*
 * @brief Minimize the sum of squares of a risk vector.
 * @param pts Vector of DistPts.
 * @param tau Average change in risk.
 * @pre The length of pts is greater than or equal to 1.
 * @pre The DistPts in pts are sorted in increasing order of value, and
 *      correspond to a valid risk vector.
 * @return A MinRes object representing all solutions with mean tau.
 */
void MinRes::minimize(const PtVector &dps, const double tau) {
  // Store vector length
  size_t n = dps.size();

  // Initialize "active" values for remeaning and minimization
  size_t i_lwr = 0;              // Lower "active" index
  size_t i_upr = n - 1;          // Upper "active" index
  double t_lwr = dps[i_lwr].val; // Lower thresholds
  double t_upr = dps[i_upr].val; // Upper thresholds
  double k_lwr = 0;              // Lower number of active indices
  double k_upr = 0;              // Upper number of active indices

  // MINIMIZE
  // Find the minimum sum of squares (SS) among all solutions with given L_1
  // distance
  double D = n * tau / 3; // Difference between required and actual mean
  double epsilon_min = std::abs(D); // Minimum budget expenditure to "fix" mean
  double Sigma = 0;                 // Change in sum of squares
  DistPt dp;                        // Distinguished point
  if (D > 0) {                      // If mean is to be increased
    double t_nxt = i_lwr == n - 1 ? t_lwr : dps[i_lwr + 1].val;
    double d = 0; // (Relying on fact that first DistPt is lower bound)
    while (D > d && i_lwr < n - 1) { // ...increase threshold on lower end
      D -= d;
      Sigma += d * (t_nxt + t_lwr);
      ++i_lwr;
      t_lwr = t_nxt;
      dp = dps[i_lwr];
      t_nxt = i_lwr == n - 1 ? t_lwr : dps[i_lwr + 1].val;
      if (dp.type == PtType::est) {
        ++k_lwr;
      } else if (dp.type == PtType::upr) {
        --k_upr;
      }
      d = k_lwr * (t_nxt - t_lwr);
    }
    // NOTE: The final iteration cannot have k_lwr == 0 unless the entire
    // budget has been exhausted (i.e., e[i] <- upr[i] for all i).
    if (i_lwr < n - 1) {
      Sigma = D == 0 ? Sigma : Sigma + D * (2 * t_lwr + D / k_lwr);
      t_lwr = D == 0 ? t_lwr : t_lwr + D / k_lwr;
    }
  } else if (D < 0) { // If mean is to be decreased...
    D = std::abs(D);
    double t_nxt = i_upr == 0 ? t_upr : dps[i_upr - 1].val;
    double d = 0; // (Relying on fact that last DistPt is upper bound)
    while (D > d && i_upr > 0) { // ... decrease threshold on upper end.
      D -= d;
      Sigma -= d * (t_nxt + t_upr);
      --i_upr;
      t_upr = t_nxt;
      dp = dps[i_upr];
      t_nxt = i_upr == 0 ? t_upr : dps[i_upr - 1].val;
      if (dp.type == PtType::est) {
        ++k_upr;
      } else if (dp.type == PtType::lwr) {
        --k_upr;
      }
      d = k_upr * (t_upr - t_nxt);
    }
    // NOTE: The final iteration cannot have k_upr == 0 unless the entire
    // budget has been exhausted (i.e., e[i] <- lwr[i] for all i).
    if (i_upr > 0) {
      Sigma = D == 0 ? Sigma : Sigma - D * (2 * t_upr - D / k_upr);
      t_upr = D == 0 ? t_upr : t_upr - D / k_upr;
    }
  } else {
    epsilon_min = 0.0;
    Sigma = 0.0;
  }

  // MINIMIZE
  // Find the minimum sum of squares (SS) among all solutions within given L_1
  // distance from the input
  // Initialize some intermediate variables
  n = i_upr - i_lwr + 1; // Number of iterations required
  resize(n);
  double delta = t_upr - t_lwr; // Difference between thresholds
  double K = 0.0;               // Coefficient of SS growth in delta

  // Store results after zero iterations
  epsilon_cum[0] = epsilon_min;
  Sigma_cum[0] = Sigma;
  delta_cum[0] = delta;
  kappa_cum[0] = K;

  // If i_lwr and i_upr are already equal, we're done
  if (n == 1) {
    return;
  }

  // Initialize the active gaps
  double Delta_lwr = k_lwr * (dps[i_lwr + 1].val - t_lwr); // Lower active gap
  double Delta_upr = k_upr * (t_upr - dps[i_upr - 1].val); // Upper active gap
  double Delta = std::min(Delta_lwr, Delta_upr); // Smaller of the gaps

  // Loop over remaining iterations until i_lwr and i_upr are adjacent
  if (n > 2) {
    for (size_t i = 1; i < n - 1; ++i) {
      // Advance one iteration
      if (Delta_lwr == Delta) {
        ++i_lwr;         // Increment lower active index
        dp = dps[i_lwr]; // Fetch the distinguished point
        t_lwr = dp.val;  // Update the lower threshold
        if (dp.type == PtType::est) {
          ++k_lwr; // Increment k_lwr if dp is est
        } else if (dp.type == PtType::upr) {
          --k_lwr; // Decrement k_lwr if dp is upr
        }
        Delta_upr = Delta_upr - Delta_lwr; // Calculate new upper active gap
        // Calculate new upper threshold
        t_upr = k_upr == 0 ? t_upr : t_upr - Delta_lwr / k_upr;
        // Calculate new lower active gap
        Delta_lwr = k_lwr * (dps[i_lwr + 1].val - t_lwr);
      } else {
        --i_upr;         // Decrement upper active index
        dp = dps[i_upr]; // Fetch the distinguished point
        t_upr = dp.val;  // Update the upper threshold
        if (dp.type == PtType::est) {
          ++k_upr; // Increment k_upr if dp is est
        } else if (dp.type == PtType::lwr) {
          --k_upr; // Decrement k_upr if dp is lwr
        }
        Delta_lwr = Delta_lwr - Delta_upr; // Calculate new lower active gap
        // Calculate new lower threshold
        t_lwr = k_lwr == 0 ? t_lwr : t_lwr + Delta_upr / k_lwr;
        // Calculate new upper active gap
        Delta_upr = k_upr * (t_upr - dps[i_upr - 1].val);
      }
      // Record results
      epsilon_cum[i] = 2 * Delta;
      Sigma_cum[i] = Delta * (K * Delta - 2 * delta);
      Delta = std::min(Delta_upr, Delta_lwr); // Next gap is minimum of gaps
      delta = t_upr - t_lwr;                  // Next diff between thresholds
      // Next coefficient for change in sum of squares
      K = k_lwr == 0 || k_upr == 0 ? 0.0 : 1 / k_lwr + 1 / k_upr;
      delta_cum[i] = delta;
      kappa_cum[i] = K == 0 ? 0.0 : 1 / K;
    }
  }
  // Final iteration, when i_lwr and i_upr are adjacent
  Delta = K == 0 ? 0.0 : (t_upr - t_lwr) / K;
  epsilon_cum[n - 1] = 2 * Delta;
  Sigma_cum[n - 1] = Delta * (K * Delta - 2 * delta);
  delta_cum[n - 1] = 0;
  kappa_cum[n - 1] = 0;
  // Take cumulative sums of epsilon and Sigma
  for (size_t i = 1; i < n; ++i) {
    epsilon_cum[i] += epsilon_cum[i - 1];
    Sigma_cum[i] += Sigma_cum[i - 1];
  }
}

/*
 * @brief Combine the output of the minimization route across levels.
 * @param res_1 Output of the minimization routine.
 * @param res_2 Output of the minimization routine.
 * @return All possible solutions for all values of delta in res1 and res2.
 */
void MinRes::combine(const MinRes &res1, const MinRes &res2) {
  // Initialize intermediate variables
  size_t n1 = res1.n;                      // Length of res1 vectors
  size_t n2 = res2.n;                      // Length of res2 vectors
  size_t i0 = 0;                           // Active index for combined
  size_t i1 = 0;                           // Active index for res1
  size_t i2 = 0;                           // Active index for res2
  double epsilon1 = res1.epsilon_cum[i1];  // Budget exhausted in res1
  double epsilon2 = res2.epsilon_cum[i2];  // Budget exhausted in res2
  double Sigma1 = res1.Sigma_cum[i1];      // Change in SS in res1
  double Sigma2 = res2.Sigma_cum[i2];      // Change in SS in res2
  double delta1 = res1.delta_cum[i1];      // Active difference in res1
  double delta2 = res2.delta_cum[i2];      // Active difference in res2
  double kappa1 = res1.kappa_cum[i1];      // Coefficient of SS growth in res1
  double kappa2 = res2.kappa_cum[i2];      // Coefficient of SS growth in res2
  double delta = std::max(delta1, delta2); // Active difference in combined
  double df_delta = 0;                     // Change in difference

  // Resize vectors
  resize(n1 + n2);

  // Iterate through both containers
  for (size_t i = 0; i < n1 + n2; i++) {
    if (delta == delta1) {               // If delta is from res1...
      kappa1 = res1.kappa_cum[i1];       // Update kappa1 with res1
      epsilon1 = res1.epsilon_cum[i1];   // Update epsilon1 with res1
      Sigma1 = res1.Sigma_cum[i1];       // Update Sigma1 with res1
      epsilon2 += 2 * kappa2 * df_delta; // Update epsilon2 with formula
      // Sigma2 is changed according to the formula
      Sigma2 += kappa2 * df_delta * (df_delta - 2 * (delta + df_delta));
      if (i1 < n1 - 1) {
        i1++;                        // Increment active index in res1
        delta1 = res1.delta_cum[i1]; // Get next value of delta1...
      } else {
        delta1 = -std::numeric_limits<double>::infinity(); // Unless it doesn't
                                                           // exist.
      }
      df_delta =
          delta - std::max(delta1, delta2); // Calculate jump for next iteration
    } else {                                // Otherwise, if from res2...
      kappa2 = res2.kappa_cum[i2];          // Update kappa2 with res2
      epsilon2 = res2.epsilon_cum[i2];      // Update epsilon2 with res2
      Sigma2 = res2.Sigma_cum[i2];          // Update Sigma2 with res2
      epsilon1 += 2 * kappa1 * df_delta;    // Update epsilon1 with formula
      // Sigma1 is changed according to the formula
      Sigma1 += kappa1 * df_delta * (df_delta - 2 * (delta + df_delta));
      if (i2 < n2 - 1) {
        i2++;                        // Increment active index in res2
        delta2 = res2.delta_cum[i2]; // Get next value of delta2...
      } else {
        delta2 = -std::numeric_limits<double>::infinity(); // Unless it doesn't
                                                           // exist.
      }
      df_delta =
          delta - std::max(delta1, delta2); // Calculate jump for next iteration
    }
    // Record the current values of the parameters
    if (!(std::abs(df_delta) < 1e-10) && std::isfinite(df_delta)) {
      delta_cum[i0] = delta;                 // Record delta jumping to
      epsilon_cum[i0] = epsilon1 + epsilon2; // Record the sum of epsilons
      Sigma_cum[i0] = Sigma1 + Sigma2;       // Record the sum of Sigmas
      kappa_cum[i0] = kappa1 + kappa2;       // Record the sum of kappas
      ++i0;
    }
    // The final value that is recorded is when both containers are empty.
    if (df_delta == std::numeric_limits<double>::infinity()) {
      delta_cum[i0] = delta;                 // Record delta jumping to
      epsilon_cum[i0] = epsilon1 + epsilon2; // Record the sum of epsilons
      Sigma_cum[i0] = Sigma1 + Sigma2;       // Record the sum of Sigmas
      kappa_cum[i0] = kappa1 + kappa2;       // Record the sum of kappas
    }
    delta = std::max(delta1, delta2); // Calculate next delta to jump to
  }
  // Correct the final value of kappa_cum to zero
  kappa_cum[i0] = 0.0;

  // Update the length of the result container
  resize(i0 + 1);
}

//////////////////////////// MinGrid Implementation ////////////////////////////

/*
 * @brief Constructor for MinGrid.
 * @param m Length of the grid.
 * @pre `m` is greater than or equal to 1.
 * @return A MinGrid object with length `m`.
 */
MinGrid::MinGrid(size_t m) : m(m), g(m){};

/*
 * @brief Constructor for MinGrid.
 * @param init_list Initializer list for the grid.
 * @pre `init_list.size()` is greater than or equal to 1.
 * @return A MinGrid object with length `init_list.size()`.
 */
MinGrid::MinGrid(std::initializer_list<double> init_list)
    : m(init_list.size()), g(init_list){};

/*
 * @brief Size of the grid.
 * @return The length of the grid.
 */
size_t MinGrid::size() const { return m; }

/*
 * @brief Subscript operator for MinGrid object.
 * @param index Index of the grid point.
 * @pre `index` is greater than or equal to 0.
 * @pre `index` is less than or equal to the length of the grid.
 * @return The value of the grid at the given index.
 */
double &MinGrid::operator[](std::size_t index) { return g[index]; }

/*
 * @brief Subscript operator for MinGrid object.
 * @param index Index of the grid point.
 * @pre `index` is greater than or equal to 0.
 * @pre `index` is less than or equal to the length of the grid.
 * @return The value of the grid at the given index.
 */
const double &MinGrid::operator[](std::size_t index) const { return g[index]; }

/*
 * @brief Evaluate a MinGrid object on a grid of its points.
 * @param res MinRes object.
 * @param m Length of the grid.
 * @param gamma Step size.
 * @pre `epsilon` is greater than or equal to 0.
 * @pre m is greater than or equal to 1.
 * @pre res has been minimized.
 * @return A MinGrid object representing the minimum sum of squares on a grid.
 */
void MinGrid::grid(const MinRes &res, const double gamma) {
  size_t i = 0;     // Pointer to position in res
  size_t n = res.n; // Length of result of minimization routine
  double epsilon_max = res.epsilon_cum[n - 1]; // Maximum budget expendable
  double epsilon_cur = res.epsilon_cum[i]; // Budget expended at "current" step
  double epsilon_itr = 0.0; // Budget exhausted up to current iteration
  // Budget expended at "next" step, difference between active values, change
  // in sum of squares, change in difference between active values, growth
  // coefficient of SS in difference between active values.
  double epsilon_nxt, delta, Sigma, Delta, kappa;

  // Early return if the MinRes is length one
  if (n == 1) {
    Sigma = res.Sigma_cum[n - 1];
    for (size_t k = 0; k < m; ++k) {
      g[k] = (epsilon_itr <= epsilon_max)
                 ? std::numeric_limits<double>::infinity()
                 : Sigma;
      epsilon_itr += gamma;
    }
    return;
  }

  epsilon_nxt = res.epsilon_cum[i + 1]; // Budget expended at "next" step
  delta = res.delta_cum[i];             // Difference between active values
  // NOTE: If the current point is feasible, then Sigma is highest value;
  // otherwise it takes the sentinel value of Inf. Due to floating point
  // errors, equality is checked with a tolerance.
  Sigma = (epsilon_itr >= epsilon_cur - 1e-10)
              ? res.Sigma_cum[i]
              : std::numeric_limits<double>::infinity();
  kappa = res.kappa_cum[i]; // Coefficient of SS growth in difference

  for (size_t k = 0; k < m; ++k) { // Loop through grid points
    while (epsilon_itr >= epsilon_nxt && i < n - 1) { // Advance to next point
      ++i;
      epsilon_cur = res.epsilon_cum[i];
      delta = res.delta_cum[i];
      Sigma = res.Sigma_cum[i];
      kappa = res.kappa_cum[i];
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
        g[k] = Sigma + Delta * (Delta / (4 * kappa) - delta);
      } else {
        g[k] = Sigma; // epsilon_max exceeded, sum of squares cannot increase
      }
    } else {
      g[k] = std::numeric_limits<double>::infinity(); // Infeasible
    }
    epsilon_itr += gamma;
  }
}

//////////////////////////// MinTree Implementation ////////////////////////////

/*
 * @brief Constructor for MinTree objects.
 * @param data Regression data.
 * @param m Length of the grid.
 * @param gamma Step size.
 * @pre `m` is greater than or equal to 1.
 * @pre `gamma` is greater than 0.
 * @return A MinTree object with 2 * `ns.size()` - 1 MinRes objects stored in
 *         the `ress` and `combs` members. The i-th MinRes object in `ress` is
 *         initialized to be the result of the minimization routine for the i-th
 *         stratum, and the i-th MinRes object in `combs` is the result of
 *         combining the i-th and (i + 1)-th MinRes objects in `ress`. (The
 *         final MinRes object in `combs` is a dummy and points to the final
 *         MinRes object in `ress`.) The minimization routine is run on the
 *         observed strata, and the results are combined appropriately.
 */
MinTree::MinTree(const RegData &data, const size_t m, const double gamma)
    : data(data), g(data.g), m(m), gamma(gamma), i_tau(g - 1), taus(g, 0),
      grid(m),
      beta_min(g - 1,
               std::vector<double>(m, std::numeric_limits<double>::infinity())),
      beta_max(g - 1, std::vector<double>(
                          m, -std::numeric_limits<double>::infinity())) {
  // Initialize MinRes result objects
  ress.reserve(2 * g);
  for (size_t i = 0; i < 2 * g; ++i)
    ress.emplace_back(data.es[i].size());

  // Initialize the MinRes combination objects
  combs.reserve(2 * g - 1);
  for (size_t i = 0; i < 2 * g - 1; ++i)
    combs.emplace_back(data.es[i].size());

  // INITIALIZE
  // Minimize the sum of squares in the observed strata
  for (size_t i = i_tau; i < 2 * g; ++i)
    ress[i].minimize(data.dps[i], 0.0);

  // Combine the results of the minimization routine
  // NOTE: The final combined MinRes object is combined from the final *two*
  // result MaxRes objects.
  combs[2 * g - 2].combine(ress[2 * g - 2], ress[2 * g - 1]);
  // NOTE: Since size_t is unsigned, the loop must be written in this way.
  for (size_t i = 2 * g - 2; i > i_tau; --i)
    combs[i - 1].combine(ress[i - 1], combs[i]);
}

/*
 * @brief Accessor function for the minimum regression coefficients.
 * @return The minimum regression coefficients observed so far.
 */
const std::vector<std::vector<double>> &MinTree::get_beta_min() const {
  return beta_min;
}

/*
 * @brief Accessor function for the maximum regression coefficients.
 * @return The maximum regression coefficients observed so far.
 */
const std::vector<std::vector<double>> &MinTree::get_beta_max() const {
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
void MinTree::remean(const std::vector<double> &taus) {
  for (size_t i = 0; i < g; ++i) {
    if (this->taus[i] != taus[i]) {
      this->taus[i] = taus[i];
      i_tau = i;
    }
  }
}

/*
 * @brief Minimize the sum of squares of a risk vector.
 * @return The minimum sum of squares of risk vectors across strata, subject to
 *         the given budget, upper and lower bounds, and mean changes in risk.
 *         (Stored in the inital MinRes object in `combs`.)
 */
void MinTree::minimize() {
  // Minimize the sum of squares in the strata that need to be updated.
  for (size_t i = 0; i <= i_tau; ++i)
    ress[i].minimize(data.dps[i], taus[i]);

  // Combine the results of the minimization routine.
  // NOTE: Since size_t is unsigned, the loop must be written in this way.
  for (size_t i = i_tau + 1; i > 0; --i)
    combs[i - 1].combine(ress[i - 1], combs[i]);

  // Store the results in the grid.
  grid.grid(combs[0], gamma);
}

/*
 * @brief Calculate the regression coefficients and compare to existing records.
 * @return The largest and smallest regression coefficients observed so far,
 *         stored in `beta_min` and `beta_max` as appropriate.
 */
void MinTree::regress() { data.regress(grid, taus, beta_min, beta_max); }
