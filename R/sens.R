#' Perform sensitivity analysis on a risk-adjusted regression
#'
#' `sens()` performs sensitivity analysis on a risk-adjusted regression by
#' computing the maximum and minimum regression coefficients consistent with the
#' data and the analyst's prior knowledge, expressed through `epsilon`, the
#' bound on the mean absolute difference between the true and estimated risks.
#' It additionally can provide bootstrapped pointwise confidence intervals for
#' the regression coefficients.
#'
#' @section Details:
#'
#' The sensitivity analysis assumes that every group contains at least one
#' observed and one unobserved individual, and that the estimated risks and
#' upper and lower bounds are "sortable," i.e., that there exists a permutation
#' of the rows such that the estimated risks and upper and lower bounds are all
#' non-decreasing within each group and observation status. If these conditions
#' are not met, the function will throw an error.
#'
#' To ensure that these conditions continue to hold, the bootstrap resamples are
#' stratified by group and observation status. As a result, in small samples,
#' the confidence intervals may be slightly narrowed, since they do not account
#' for uncertainty in the number of individuals in each group, and the number of
#' observed and unobserved individuals within each group.
#'
#' @param df The data frame containing the data.
#' @param group_col The name of the column containing the group labels. This
#' column should be a factor or coercible to a factor.
#' @param obs_col The name of the column containing whether or not the outcome
#' was observed. This column should be a logical or coercible to a logical.
#' @param p_col The name of the column containing the estimated risks. These
#' risks should be expressed on the probability scale, i.e., be between 0 and 1.
#' @param base_group The name of the base group. This group will be used as the
#' reference group in the regression.
#' @param epsilon The bound on the mean absolute difference between the true and
#' estimated risks.
#' @param lwr_col The name of the column containing the lower bounds on the true
#' risk. (Defaults to 0 for all observations.)
#' @param upr_col The name of the column containing the upper bounds on the true
#' risk. (Defaults to 1 for all observations.)
#' @param eta The step size for the grid search.  Note that while steps are
#' taken at the group level, the step size is expressed at the level of change
#' in average risk *across the entire population*. In other words, smaller
#' groups will have proportionally larger steps. (Defaults to 0.01.)
#' @param m The grid size for the maximization approximation. (Defaults to
#' `101`.)
#' @param N The number of bootstrap resamples to use to compute pointwise
#' confidence intervals. (Defaults to 0, which performs no bootstrap.)
#' @param alpha The confidence level for the pointwise confidence intervals.
#' (Defaults to 0.05.)
#' @param chunk_size The number of repetitions to perform in each chunk when run
#' in parallel. Larger chunk sizes make it less likely that separate threads
#' will block on each other, but also make it more likely that the threads will
#' finish at different times. (Defaults to 100.)
#' @param n_threads The number of threads to use when running in parallel.
#' (Defaults to 1, i.e., serial execution.)
#' @return A data frame containing the following columns:
#' * `epsilon`: Values of epsilon ranging from 0 to the input value of `epsilon`
#' in `m` steps.
#' * `beta_min_{group}`: The minimum value of the regression coefficient for the
#' group `group`. (Note that the base group is not included in this list.)
#' * `beta_max_{group}`: The maximum value of the regression coefficient for the
#' group `group`. (Note that the base group is not included in this list.)
#' * (**If `N > 0`**) `beta_min_{group}_{alpha/2}`: The `alpha/2` quantile of
#' the bootstrap distribution of the minimum value of the regression coefficient
#' for group `group`. (Note that the base group is not included in this list.)
#' * (**If `N > 0`**) `beta_min_{group}_{1 - alpha/2}`: The `1 - alpha/2`
#' quantile of the bootstrap distribution of the minimum value of the regression
#' coefficient for group `group`. (Note that the base group is not included in
#' this list.)
#' * (**If `N > 0`**) `beta_max_{group}_{alpha/2}`: The `alpha/2` quantile of
#' the bootstrap distribution of the maximum value of the regression coefficient
#' for group `group`. (Note that the base group is not included in this list.)
#' * (**If `N > 0`**) `beta_max_{group}_{1 - alpha/2}`: The `1 - alpha/2`
#' quantile of the bootstrap distribution of the maximum value of the regression
#' coefficient for group `group`. (Note that the base group is not included in
#' this list.)
#' @examples
#' # Generate some data
#' set.seed(1)
#' df <- tibble::tibble(
#'   group = factor(
#'     sample(c("a", "b"), 1000, replace = TRUE),
#'     levels = c("a", "b")
#'   ),
#'   p = runif(1000)^2,
#'   frisked = runif(1000) < p + 0.1 * (group != "a")
#' )
#'
#' # Compute the sensitivity analysis
#' sens(df, group, frisked, p, "a", 0.1)
#'
#' # Search over a finer grid
#' sens(df, group, frisked, p, "a", 0.1, eta = 0.001)
#'
#' # Increase the accuracy of the maximization approximation
#' sens(df, group, frisked, p, "a", 0.1, m = 1001)
#'
#' \donttest{
#' # Calculate 90% pointwise confidence intervals
#' sens(df, group, frisked, p, "a", 0.1, N = 1000, alpha = 0.1)
#'
#' # Run in parallel, adjusting the chunk size to avoid blocking
#' sens(df, group, frisked, p, "a", 0.1, n_threads = 2, eta = 0.0001,
#'      chunk_size = 1000)
#' }
#' @importFrom dplyr rename mutate select arrange n_distinct anti_join summarize
#' group_by filter bind_cols slice_sample pull
#' @importFrom rlang enquo as_name abort quo_is_null quo_text warn set_names
#' @importFrom glue glue
#' @importFrom purrr map_dfr
#' @importFrom tidyr expand_grid
#' @importFrom vctrs vec_is
#' @importFrom tibble tibble
#' @export
sens <- function(df, group_col, obs_col, p_col, base_group, epsilon,
                 lwr_col = NULL, upr_col = NULL, eta = 0.01,
                 m = 101L, N = 0L, alpha = 0.05, chunk_size = 100L,
                 n_threads = 1L) {
  # Get the user-provided columns.
  group_col <- rlang::enquo(group_col)
  obs_col <- rlang::enquo(obs_col)
  p_col <- rlang::enquo(p_col)
  lwr_col <- rlang::enquo(lwr_col)
  upr_col <- rlang::enquo(upr_col)

  # Check that the user-provided columns can be coerced to strings, if they are
  # not null.
  tryCatch(
    {
      rlang::as_name(group_col)
    },
    error = function(e) {
      rlang::abort(glue::glue(
        "Invalid column name in `group_col`.",
        .sep = " "
      ))
    }
  )
  tryCatch(
    {
      rlang::as_name(obs_col)
    },
    error = function(e) {
      rlang::abort(glue::glue(
        "Invalid column name in `obs_col`.",
        .sep = " "
      ))
    }
  )
  tryCatch(
    {
      rlang::as_name(p_col)
    },
    error = function(e) {
      rlang::abort(glue::glue(
        "Invalid column name in `p_col`.",
        .sep = " "
      ))
    }
  )
  if (!rlang::quo_is_null(lwr_col)) {
    tryCatch(
      {
        rlang::as_name(lwr_col)
      },
      error = function(e) {
        rlang::abort(glue::glue(
          "Invalid column name in `lwr_col`.",
          .sep = " "
        ))
      }
    )
  }
  if (!rlang::quo_is_null(upr_col)) {
    tryCatch(
      {
        rlang::as_name(upr_col)
      },
      error = function(e) {
        rlang::abort(glue::glue(
          "Invalid column name in `upr_col`.",
          .sep = " "
        ))
      }
    )
  }

  # Check that `df` is a data frame
  if (!is.data.frame(df)) {
    rlang::abort("Argument `df` must be a data frame.")
  }

  # Check that all of the columns are in the dataframe
  if (any(!c(
    rlang::as_name(group_col),
    rlang::as_name(obs_col),
    rlang::as_name(p_col)
  ) %in% colnames(df))) {
    rlang::abort(glue::glue(
      "Columns `group_col` ({ rlang::as_name(group_col)}),",
      "`obs_col` ({ rlang::as_name(obs_col)}),",
      "and `p_col` ({ rlang::as_name(p_col)})",
      "must be present in `df`",
      "(cols: { paste(colnames(df), collapse = ', ') }).",
      .sep = " "
    ))
  }

  # If specified, lwr_col and upr_col must be in the dataframe
  if (
    !rlang::quo_is_null(lwr_col) && !rlang::as_name(lwr_col) %in% colnames(df)
  ) {
    rlang::abort(glue::glue(
      "Column `lwr_col` ({ rlang::as_name(lwr_col) }) must be present in `df`",
      "(cols: { paste(colnames(df), collapse = ', ') }) if specified.",
      .sep = " "
    ))
  }
  if (
    !rlang::quo_is_null(upr_col) && !rlang::as_name(upr_col) %in% colnames(df)
  ) {
    rlang::abort(glue::glue(
      "Column `upr_col` ({ rlang::as_name(upr_col) }) must be present in `df`",
      "(cols: { paste(colnames(df), collapse = ', ') }) if specified.",
      .sep = " "
    ))
  }

  df <- df %>%
    dplyr::rename(
      group = !!group_col,
      obs = !!obs_col,
      p = !!p_col
    )

  # Add in lwr and upper if not present
  if (!rlang::quo_is_null(lwr_col)) {
    df <- df %>%
      dplyr::rename(lwr = !!lwr_col)
  } else {
    df <- df %>%
      dplyr::mutate(lwr = rep(0, nrow(df)))
  }
  if (!rlang::quo_is_null(upr_col)) {
    df <- df %>%
      dplyr::rename(upr = !!upr_col)
  } else {
    df <- df %>%
      dplyr::mutate(upr = rep(1, nrow(df)))
  }

  df <- df %>%
    dplyr::select(group, obs, p, lwr, upr)

  ################################# TYPE CHECKS ################################

  # Check that there are no missing values.
  if (any(is.na(df))) {
    rlang::abort(c(
      "There are missing values in the data frame.",
      "x" = glue::glue(
        "No risk estimate, group, or observation status, or lower or upper",
        "bound (if provided) can be missing.",
        .sep = " "
      )
    ))
  }

  # Check that group is coercible to factor (i.e., factor, character, numeric,
  # logical)
  if (with(
    df,
    !is.factor(group) &&
      !vctrs::vec_is(group, character()) &&
      !vctrs::vec_is(group, numeric()) &&
      !vctrs::vec_is(group, integer()) &&
      !vctrs::vec_is(group, logical())
  )) {
    rlang::abort(c(
      glue::glue("Invalid group column '{ rlang::quo_text(group_col) }'."),
      "x" = "Column  must be a factor, character, numeric, or logical vector."
    ))
  }

  # Check that base group occurs in the group column.
  if (length(base_group) != 1 || with(df, !base_group %in% group)) {
    rlang::abort(glue::glue(
      "Base group must occur in the group column."
    ))
  } else {
    base_group <- as.character(base_group)
  }

  # Check that there are at least two groups.
  if (with(df, length(unique(group)) <= 1)) {
    rlang::abort(glue::glue(
      "There must be at least two distinct groups in the group column."
    ))
  }

  # Check that obs is logical or coercible
  if (with(
    df,
    !vctrs::vec_is(obs, logical()) && (
      !vctrs::vec_is(obs, numeric()) || !all(obs %in% c(0, 1))
    )
  )) {
    rlang::abort(c(
      glue::glue(
        "Invalid observation column '{ rlang::quo_text(obs_col) }'."
      ),
      "x" = paste(
        "Column must be a logical vector, or a numeric or character",
        "vector taking only the values 0 and 1."
      )
    ))
  }

  # Check that p is numeric and between 0 and 1.
  if (with(df, !vctrs::vec_is(p, numeric()) || any(p < 0 | p > 1))) {
    rlang::abort(c(
      glue::glue("Invalid risk column '{ rlang::quo_text(p_col) }'."),
      "x" = paste(
        "column must be a numeric vector taking only values between 0 and",
        "1, inclusive."
      )
    ))
  }

  # Check that lwr is numeric and between 0 and 1.
  if (with(df, !vctrs::vec_is(lwr, numeric()) || any(lwr < 0 | lwr > 1))) {
    rlang::abort(c(
      glue::glue("Invalid lower bound column '{ rlang::quo_text(lwr_col) }'."),
      "x" = paste(
        "column must be a numeric vector taking only values between 0 and",
        "1, inclusive."
      )
    ))
  }

  # Check that upr is numeric and between 0 and 1.
  if (with(df, !vctrs::vec_is(upr, numeric()) || any(upr < 0 | upr > 1))) {
    rlang::abort(c(
      glue::glue("Invalid upper bound column '{ rlang::quo_text(upr_col) }'."),
      "x" = paste(
        "column must be a numeric vector taking only values between 0 and",
        "1, inclusive."
      )
    ))
  }

  # Check that epsilon is numeric and positive.
  if (
    !vctrs::vec_is(epsilon, numeric()) || length(epsilon) != 1 || epsilon <= 0
  ) {
    rlang::abort(c(
      glue::glue("Invalid epsilon ({ epsilon })."),
      "x" = "Argument `epsilon` must be a single, positive numeric constant."
    ))
  } else {
    if (epsilon > 1) {
      rlang::warn(c(
        glue::glue("Argument `epsilon` ({ epsilon }) was greater than 1."),
        "i" = "`epsilon` will be set to 1, its maximum meaningful value."
      ))
      epsilon <- 1
    }
    epsilon <- as.double(epsilon)
  }

  # Check that eta is numeric and positive, and less than epsilon.
  if (!vctrs::vec_is(eta, numeric()) || length(eta) != 1 || eta <= 0) {
    rlang::abort(c(
      glue::glue("Invalid eta ({ eta })."),
      "x" = paste(
        "Argument `eta` must be a single, positive numeric constant less than",
        "`epsilon`."
      )
    ))
  } else {
    if (eta > epsilon / 3) {
      rlang::warn(c(
        glue::glue(
          "Argument `eta` ({ eta }) is large relative to `epsilon`",
          "({ epsilon }).",
          .sep = " "
        ),
        "i" = "Results may be unreliable."
      ))
      eta <- epsilon
    }
    eta <- as.double(eta)
  }

  # Check that `m` is a positive whole number and at least two.
  if (
    !(
      vctrs::vec_is(m, numeric())
      || vctrs::vec_is(m, integer())
    )
    || length(m) != 1
    || m < 2
  ) {
    rlang::abort(c(
      glue::glue("Invalid `m` ({ m })."),
      "x" = "Argument `m` must be a single positive integer, at least two."
    ))
  } else {
    if (m != floor(m)) {
      rlang::warn(c(
        glue::glue("Argument `m` ({ m }) was not a whole number."),
        "i" = "Argument `m` will be rounded down to the nearest whole number."
      ))
    }
    m <- as.integer(m)
  }

  # Check that `n_threads` is a positive whole number.
  if (
    !(
      vctrs::vec_is(n_threads, numeric())
      || vctrs::vec_is(n_threads, integer())
    )
    || length(n_threads) != 1
    || n_threads < 1
  ) {
    rlang::abort(c(
      glue::glue("Invalid `n_threads` ({ n_threads })."),
      "x" = "Argument `n_threads` must be a single positive integer."
    ))
  } else {
    if (n_threads != floor(n_threads)) {
      rlang::warn(c(
        glue::glue(
          "Argument `n_threads` ({ n_threads }) was not a whole number."
        ),
        "i" = paste(
          "Argument `n_threads` will be rounded down to the nearest whole",
          "number."
        )
      ))
    }
    if (n_threads > parallel::detectCores()) {
      rlang::warn(c(
        glue::glue(
          "Argument `n_threads` ({ n_threads }) was greater than the number",
          "of cores.",
          .sep = " "
        ),
        "i" = "Argument `n_threads` will be set to the number of cores."
      ))
      n_threads <- parallel::detectCores()
    }
    n_threads <- as.integer(n_threads)
  }

  # Check that `chunk_size` is a positive whole number.
  if (
    !(
      vctrs::vec_is(chunk_size, numeric())
      || vctrs::vec_is(chunk_size, integer())
    )
    || length(chunk_size) != 1
    || chunk_size < 1
  ) {
    rlang::abort(c(
      glue::glue("Invalid `chunk_size` ({ chunk_size })."),
      "x" = "Argument `chunk_size` must be a single positive integer."
    ))
  } else {
    if (chunk_size != floor(chunk_size)) {
      rlang::warn(c(
        glue::glue(
          "Argument `chunk_size` ({ chunk_size }) was not a whole number."
        ),
        "i" = paste(
          "Argument `chunk_size` will be rounded down to the nearest whole",
          "number."
        )
      ))
    }
    chunk_size <- as.integer(chunk_size)
  }

  # Check that `N` is a non-negative whole number.
  if (
    !(vctrs::vec_is(N, numeric()) || vctrs::vec_is(N, integer()))
    || length(N) != 1
    || N < 0
  ) {
    rlang::abort(c(
      glue::glue("Invalid `N` ({ N })."),
      "x" = "Argument `N` must be a single non-negative integer."
    ))
  } else {
    if (N != floor(N)) {
      rlang::warn(c(
        glue::glue("Argument `N` ({ N }) was not a whole number."),
        "i" = "Argument `N` will be rounded down to the nearest whole number."
      ))
    }
    N <- as.integer(N)
  }

  # Check that `alpha  is numeric and between 0 and 1.
  if (
    !vctrs::vec_is(alpha, numeric())
    || length(alpha) != 1
    || alpha <= 0
    || alpha >= 1
  ) {
    rlang::abort(c(
      glue::glue("Invalid alpha ({ alpha })."),
      "x" = paste(
        "Argument `alpha` must be a single numeric constant between zero and",
        "one."
      )
    ))
  } else {
    if (N > 0 && alpha * N < 1) {
      rlang::abort(c(
        glue::glue(
          "Argument `alpha` ({ alpha }) is too small to produce any confidence",
          "intervals.",
          .sep = " "
        ),
        "i" = "Argument `alpha` should be at least 1 / N."
      ))
    }
    if (N > 0 && alpha * N < 5) {
      rlang::warn(c(
        glue::glue(
          "Argument `alpha` ({ alpha }) is small relative to `N` ({ N })."
        ),
        "i" = "Confidence intervals may be unreliable."
      ))
    }
    alpha <- as.double(alpha)
  }

  ########################### DATA INTEGRITY CHECKS ############################

  df <- df %>%
    dplyr::mutate(
      group = forcats::fct_relevel(forcats::as_factor(group), {{ base_group }}),
      obs = as.logical(obs),
      p = as.numeric(p),
    ) %>%
    dplyr::arrange(obs, group, p, upr, lwr)

  # Check if there are any unused levels and drop them.
  if (any(!levels(df$group) %in% levels(forcats::fct_drop(df$group)))) {
    rlang::warn(glue::glue(
      "Some groups had no observations. Those groups will be dropped from",
      "the analysis.",
      .sep = " "
    ))
    df <- df %>%
      dplyr::mutate(group = forcats::fct_drop(group))
  }

  # Record the number of groups
  G <- with(df, dplyr::n_distinct(group))

  # Check if every group has at least one observed and one unobserved individual
  strata <- with(
    df,
    tidyr::expand_grid(group = levels(group), obs = c(TRUE, FALSE))
  )
  if (nrow(dplyr::anti_join(strata, df, by = c("group", "obs"))) > 0) {
    rlang::abort(glue::glue(
      "Every group must have at least one observed and one unobserved",
      "individual.",
      .sep = " "
    ))
  }

  # Check if the risk estimates are between zero and one
  if (with(df, !all(lwr <= p & p <= upr))) {
    rlang::abort(c(
      glue::glue("The estimated risks (p) are invalid.",
        "*" = "Some values were not between lwr and upr."
      )
    ))
  }

  # Check if the bounds are correctly sorted
  valid <- df %>%
    dplyr::group_by(group, obs) %>%
    dplyr::summarize(
      valid = c(
        all(lwr == cummax(lwr)) & all(p == cummax(p)) & all(upr == cummax(upr))
      )
    ) %>%
    dplyr::pull(valid)
  if (!all(valid)) {
    rlang::abort(c(
      glue::glue("The bounds (lwr, p, upr) are not sortable."),
      "* The bounds and risks cannot be sorted in the same order."
    ))
  }

  ############################### CALL C++ CODE ################################

  # Calculate the upper and lower bounds for the iterator in each group.
  iter_bounds <- df %>%
    dplyr::filter(!obs) %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(
      lwr = mean(lwr) - mean(p),
      upr = max(upr) - mean(p)
    )

  # Generate the estimated betas.
  raw_betas <- with(
    df,
    sens_(
      p, lwr, upr, iter_bounds$lwr, iter_bounds$upr, group, obs, epsilon,
      eta, m, chunk_size, n_threads
    )
  )

  # Restructure the betas as a data frame.
  groups <- with(df, levels(group))
  res <- dplyr::bind_cols(
    epsilon = seq(0, epsilon, length.out = m),
    rlang::set_names(raw_betas[[1]], paste0("beta_min_", groups[2:G])),
    rlang::set_names(raw_betas[[2]], paste0("beta_max_", groups[2:G]))
  )

  # If N > 0, generate the bootstrap resampled confidence intervals.
  if (N > 0) {
    boot_res <- purrr::map_dfr(1:N, .progress = "Resamples", ~ {
      # Generate the bootstrap resample.
      boot_df <- df %>%
        dplyr::group_by(group, obs) %>%
        dplyr::slice_sample(prop = 1, replace = TRUE) %>%
        dplyr::ungroup()

      # Calculate the upper and lower bounds for the iterator in each group.
      boot_iter_bounds <- boot_df %>%
        dplyr::filter(!obs) %>%
        dplyr::group_by(group) %>%
        dplyr::summarize(
          lwr = mean(lwr) - mean(p),
          upr = max(upr) - mean(p)
        )

      # Generate the estimated betas.
      boot_raw_betas <- with(
        boot_df,
        sens_(
          p, lwr, upr, boot_iter_bounds$lwr, boot_iter_bounds$upr, group, obs,
          epsilon, eta, m, chunk_size, n_threads
        )
      )

      # Restructure the betas as a data frame.
      dplyr::bind_cols(
        i = .x,
        epsilon = seq(0, epsilon, length.out = m),
        rlang::set_names(boot_raw_betas[[1]], paste0("beta_min_", groups[2:G])),
        rlang::set_names(boot_raw_betas[[2]], paste0("beta_max_", groups[2:G]))
      )
    })

    # Take the quantiles of the bootstrap resamples and recombine with the
    # original results.
    fns <- list(
      ~ quantile(., probs = alpha / 2),
      ~ quantile(., probs = 1 - alpha / 2)
    )
    fns <- rlang::set_names(fns, c(
      sprintf("%04.1f", 100 * alpha / 2),
      sprintf("%04.1f", 100 * (1 - alpha / 2))
    ))
    res <- boot_res %>%
      dplyr::group_by(epsilon) %>%
      dplyr::reframe(dplyr::across(tidyselect::starts_with("beta_"), fns)) %>%
      dplyr::left_join(res, by = "epsilon")
  }

  # Return the results.
  res
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("group", "obs", "p", "lwr", "upr"))
}
