#include "sens.h"

/*
 * @brief Run the risk-adjusted regression sensitivity analysis on the data.
 * @param e The estimated risks.
 * @param lwr Lower bounds for risks.
 * @param upr Upper bounds for risks.
 * @param lwr_iter The lower bound on the mean risk of unobserved individuals in
 * each group.
 * @param upr_iter The upper bound on the mean risk of unobserved individuals in
 * each group.
 * @param grp The group to which each observation belongs.
 * @param a Whether an individual is observed.
 * @param epsilon L_1 bound on the differences. (Aggregate.)
 * @param eta The distance between gridpoints for the mean probability for
 * unobserved individuals. (The algorithm scales like O(1/eta^|g|).)
 * @param m Size of grid for maximum approximation. (The
 * algorithm scales like O(m^2).)
 * @param chunk_size The number of observations to process in each thread.
 * @param num_threads The number of threads to use.
 * @pre e, lwr, upr, g, and a must all have the same length.
 * @pre lwr .≤ e .≤ upr.
 * @pre lwr, p, and upr must be sorted in increasing order.
 * @pre The set of all values of g must take the form 1:G for some G.
 * @pre lwr_iter and upr_iter must have length G.
 * @pre epsilon must be non-negative.
 * @pre eta must be strictly positive unless epsilon == 0.
 * @pre a must be true for at least one member of each group and false for at
 * least one member of each group.
 * @pre The data should be laid out as follows. (1) `a` should be ascending. (2)
 * within each value of `a`, `g` should be ascending. (3) within each value of
 * `a` and `g`, `e`, `lwr`, and `upr` should be ascending.
 * @pre There should be at least two groups.
 * @pre chunk_size must be at least one.
 * @pre num_threads must be at least one.
 * @return The maximum and minimum values of $\beta$ subject to the constraints.
 */
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
sens(const std::vector<double> &e, const std::vector<double> &lwr,
     const std::vector<double> &upr, const std::vector<double> &lwr_iter,
     const std::vector<double> &upr_iter, const std::vector<int> &grp,
     const std::vector<bool> &a, const double epsilon, const double eta,
     const size_t m, const size_t chunk_size, const int num_threads) {
  // Number of observations
  size_t N = e.size();
  // Number of groups
  size_t g = *std::max_element(grp.begin(), grp.end());
  // Number of gridpoints
  double gamma = epsilon / (m - 1);

  // Group observations
  std::vector<std::vector<double>> es(2 * g);
  std::vector<std::vector<double>> lwrs(2 * g);
  std::vector<std::vector<double>> uprs(2 * g);

  // Create the groups
  for (size_t i = 0; i < e.size(); ++i) {
    if (!a[i]) {
      es[grp[i] - 1].push_back(e[i]);
      lwrs[grp[i] - 1].push_back(lwr[i]);
      uprs[grp[i] - 1].push_back(upr[i]);
    } else {
      es[g + grp[i] - 1].push_back(e[i]);
      lwrs[g + grp[i] - 1].push_back(lwr[i]);
      uprs[g + grp[i] - 1].push_back(upr[i]);
    }
  }

  // Create the regression data
  RegData data(es, lwrs, uprs);

  // Create the eta vector
  std::vector<double> etas(g);
  std::transform(es.begin(), es.begin() + g, etas.begin(),
                 [N, &eta](auto &e) { return eta * N / e.size(); });

  // Create the iterator and mutex
  EllOneIter iter(g, eta, etas, lwr_iter, upr_iter, epsilon);
  std::mutex mtx;

  // For each thread, create a chunk of observations and a min and max tree
  std::vector<std::thread> threads;
  std::vector<TauChunk> chunks;
  for (int i = 0; i < num_threads; ++i)
    chunks.emplace_back(chunk_size, mtx, iter);
  std::vector<MinTree> min_trees(num_threads, MinTree(data, m, N * gamma));
  std::vector<MaxTree> max_trees(num_threads, MaxTree(data, m, N * gamma));

  // Anonymous function called within each thread
  auto compute = [](TauChunk &chunk, MinTree &min_tree, MaxTree &max_tree) {
    // Check if the iter is valid.
    do {
      // For each set of means in the chunk...
      for (const auto &tau : chunk) {
        // ... compute the minimum and maximum
        min_tree.remean(tau);
        min_tree.minimize();
        min_tree.regress();
        max_tree.remean(tau);
        max_tree.maximize();
        max_tree.regress();
      }
    } while (chunk.refill() > 0);
  };

  // Create the threads and then join them
  for (int i = 0; i < num_threads; ++i) {
    threads.emplace_back(compute, std::ref(chunks[i]), std::ref(min_trees[i]),
                         std::ref(max_trees[i]));
  }
  for (std::thread &thread : threads) {
    thread.join();
  }

  // Return the results
  std::vector<std::vector<double>> beta_min(
      g - 1, std::vector<double>(m, std::numeric_limits<double>::infinity()));
  std::vector<std::vector<double>> beta_max(
      g - 1, std::vector<double>(m, -std::numeric_limits<double>::infinity()));
  for (const auto &min_tree : min_trees) {
    for (size_t j = 0; j < g - 1; ++j) {
      for (size_t i = 0; i < m; ++i) {
        beta_min[j][i] =
            std::min(beta_min[j][i], min_tree.get_beta_min()[j][i]);
        beta_max[j][i] =
            std::max(beta_max[j][i], min_tree.get_beta_max()[j][i]);
      }
    }
  }
  for (const auto &max_tree : max_trees) {
    for (size_t j = 0; j < g - 1; ++j) {
      for (size_t i = 0; i < m; ++i) {
        beta_min[j][i] =
            std::min(beta_min[j][i], max_tree.get_beta_min()[j][i]);
        beta_max[j][i] =
            std::max(beta_max[j][i], max_tree.get_beta_max()[j][i]);
      }
    }
  }

  return std::tuple<decltype(beta_min), decltype(beta_max)>(beta_min, beta_max);
}

/*
 * @brief Run the risk-adjusted regression sensitivity analysis on the data. (R
 * wrapper.)
 * @param e The estimated risks.
 * @param lwr Lower bounds for risks.
 * @param upr Upper bounds for risks.
 * @param lwr_iter The lower bound on the mean risk of unobserved individuals in
 * each group.
 * @param upr_iter The upper bound on the mean risk of unobserved individuals in
 * each group.
 * @param grp The group to which each observation belongs.
 * @param a Whether an individual is observed.
 * @param epsilon L_1 bound on the differences. (Aggregate.)
 * @param eta The distance between gridpoints for the mean probability for
 * unobserved individuals. (The algorithm scales like O(1/eta^|g|).)
 * @param m Size of grid for maximum approximation. (The
 * algorithm scales like O(m^2).)
 * @param chunk_size The number of observations to process in each thread.
 * @param num_threads The number of threads to use.
 * @pre e, lwr, upr, g, and a must all have the same length.
 * @pre lwr .≤ e .≤ upr.
 * @pre lwr, p, and upr must be sorted in increasing order.
 * @pre The set of all values of g must take the form 1:G for some G.
 * @pre lwr_iter and upr_iter must have length G.
 * @pre epsilon must be non-negative.
 * @pre eta must be strictly positive unless epsilon == 0.
 * @pre a must be true for at least one member of each group and false for at
 * least one member of each group.
 * @pre The data should be laid out as follows. (1) `a` should be ascending. (2)
 * within each value of `a`, `g` should be ascending. (3) within each value of
 * `a` and `g`, `e`, `lwr`, and `upr` should be ascending.
 * @pre There should be at least two groups.
 * @pre chunk_size must be at least one.
 * @pre num_threads must be at least one.
 * @return The maximum and minimum values of $\beta$ subject to the constraints.
 */
[[cpp11::register]] cpp11::list
sens_(const cpp11::doubles &e, const cpp11::doubles &lwr,
      const cpp11::doubles &upr, const cpp11::doubles &lwr_iter,
      const cpp11::doubles &upr_iter, const cpp11::integers &grp,
      const cpp11::logicals &a, const double epsilon, const double eta,
      const int m, const int chunk_size, const int num_threads) {
  using namespace cpp11::literals;

  // Convert the vectors to `std::vector`
  std::vector<double> e_vec(e.begin(), e.end());
  std::vector<double> lwr_vec(lwr.begin(), lwr.end());
  std::vector<double> upr_vec(upr.begin(), upr.end());
  std::vector<double> lwr_iter_vec(lwr_iter.begin(), lwr_iter.end());
  std::vector<double> upr_iter_vec(upr_iter.begin(), upr_iter.end());
  std::vector<int> grp_vec(grp.begin(), grp.end());
  std::vector<bool> a_vec(a.begin(), a.end());

  std::vector<std::vector<double>> beta_min, beta_max;
  std::tie<std::vector<std::vector<double>>, std::vector<std::vector<double>>>(
      beta_min, beta_max) =
      sens(e_vec, lwr_vec, upr_vec, lwr_iter_vec, upr_iter_vec, grp_vec, a_vec,
           epsilon, eta, m, chunk_size, num_threads);

  // Return the results as an `R` list object
  cpp11::writable::list beta_min_list(beta_min.size());
  cpp11::writable::list beta_max_list(beta_max.size());
  for (size_t i = 0; i < beta_min.size(); ++i) {
    beta_min_list[i] = (cpp11::doubles)beta_min[i];
    beta_max_list[i] = (cpp11::doubles)beta_max[i];
  }
  cpp11::writable::list out(
      {"beta_min"_nm = beta_min_list, "beta_max"_nm = beta_max_list});

  return out;
}
