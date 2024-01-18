test_that("OLS matches sens", {
  n <- 1000
  p <- runif(n)^2
  group <- sample(1:3, n, replace = TRUE)
  obs <- runif(n) <= p + (group == 2) * 0.1
  df <- dplyr::tibble(p, group, obs) %>%
    dplyr::mutate(group = factor(group, levels = 1:3))
  lm_fit <- lm(obs ~ p + group, data = df) %>%
    broom::tidy() %>%
    dplyr::filter(stringr::str_detect(term, "group")) %>%
    dplyr::pull(estimate) %>%
    rep(2)
  sens_fit <- sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200) %>%
    dplyr::filter(epsilon == 0) %>%
    dplyr::select(tidyselect::starts_with("beta")) %>%
    tidyr::pivot_longer(cols = tidyselect::starts_with("beta"), names_to = "term", values_to = "estimate") %>%
    dplyr::pull(estimate)
  expect_equal(lm_fit, sens_fit)
})

test_that("Error if columns not in dataframe", {
  n <- 1000
  p <- runif(n)^2
  group <- sample(1:3, n, replace = TRUE)
  obs <- runif(n) <= p + (group == 2) * 0.1
  df <- dplyr::tibble(p, group, obs) %>%
    dplyr::mutate(group = factor(group, levels = 1:3))
  expect_error(
    sens(df, group_, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "must be present"
  )
  expect_error(
    sens(df, group, obs_, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "must be present"
  )
  expect_error(
    sens(df, group, obs, p_, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "must be present"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, lwr_col = lwr_, n_threads = 1, chunk_size = 200),
    "must be present"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, upr_col = upr_, n_threads = 1, chunk_size = 200),
    "must be present"
  )

  # Pass non-coercible column names
  expect_error(
    sens(df, 2.0, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "Invalid column name"
  )
  expect_error(
    sens(df, group, 2.0, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "Invalid column name"
  )
  expect_error(
    sens(df, group, obs, 2.0, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "Invalid column name"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, lwr_col = 2.0, n_threads = 1, chunk_size = 200),
    "Invalid column name"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, upr_col = 2.0, n_threads = 1, chunk_size = 200),
    "Invalid column name"
  )
})

test_that("Arguments have the right type", {
  n <- 1000
  p <- runif(n)^2
  group <- sample(1:3, n, replace = TRUE)
  obs <- runif(n) <= p + (group == 2) * 0.1
  df <- dplyr::tibble(p, group, obs) %>%
    dplyr::mutate(group = factor(group, levels = 1:3))

  df0 <- as.list(df)
  expect_error(
    sens(df0, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "data frame"
  )

  # Change group to a list column
  df1 <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(group = list(group)) %>%
    dplyr::ungroup()
  expect_error(
    sens(df1, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "must be a factor"
  )

  # Obs is not logical or coercible to logical
  df3 <- df %>%
    dplyr::mutate(obs = as.character(obs))
  expect_error(
    sens(df3, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "must be a logical"
  )

  # risk is not numeric
  df4 <- df %>%
    dplyr::mutate(p = as.character(p))
  expect_error(
    sens(df4, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "must be a numeric"
  )

  # Lwr and upr are not numeric
  df6 <- df %>%
    dplyr::mutate(lwr_ = as.character(0))
  expect_error(
    sens(df6, group, obs, p, 1, 0.5, eta = 0.01, lwr_col = lwr_, n_threads = 1, chunk_size = 200),
    "must be a numeric"
  )
  df7 <- df %>%
    dplyr::mutate(upr_ = as.character(1))
  expect_error(
    sens(df7, group, obs, p, 1, 0.5, eta = 0.01, upr_col = upr_, n_threads = 1, chunk_size = 200),
    "must be a numeric"
  )

  # Check that epsilon is numeric
  expect_error(
    sens(df, group, obs, p, 1, "a", eta = 0.01, n_threads = 1, chunk_size = 200),
    "must be a single, positive numeric"
  )

  # Check that eta is numeric
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = "a", n_threads = 1, chunk_size = 200),
    "must be a single, positive numeric"
  )

  # Check that m is numeric
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, m = "a", n_threads = 1, chunk_size = 200),
    "must be a single positive integer"
  )

  # Check that n_threads is numeric
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = "a", chunk_size = 200),
    "single positive integer"
  )

  # Check that chunk size is numeric
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = "a"),
    "single positive integer"
  )

  # Check that N is numeric
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200, N = "a"),
    "single non-negative integer"
  )

  # Check that alpha is numeric
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200, alpha = "a"),
    "must be a single numeric"
  )
})

test_that("Errors when argument integrity violated", {
  n <- 1000
  p <- runif(n)^2
  group <- sample(1:3, n, replace = TRUE)
  obs <- runif(n) <= p + (group == 2) * 0.1
  df <- dplyr::tibble(p, group, obs) %>%
    dplyr::mutate(group = factor(group, levels = 1:3))

  # Use a basegroup that is not a vector
  expect_error(
    sens(df, group, obs, p, list(), 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "must occur in the group column"
  )

  # Use a base group that is not in the data
  expect_error(
    sens(df, group, obs, p, 4, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "Base group"
  )

  # Base group is not a single value
  expect_error(
    sens(df, group, obs, p, c(1, 2), 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "Base group"
  )

  # Only one group
  df2 <- df %>%
    dplyr::filter(group == 1)
  expect_error(
    sens(df2, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "two distinct groups"
  )

  # Risk is not between 0 and 1
  df5 <- df %>%
    dplyr::mutate(p = p + 2)
  expect_error(
    sens(df5, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "between 0 and 1"
  )

  # Lwr and upr are not between 0 and 1
  df8 <- df %>%
    dplyr::mutate(lwr_ = -1)
  expect_error(
    sens(df8, group, obs, p, 1, 0.5, eta = 0.01, lwr_col = lwr_, n_threads = 1, chunk_size = 200),
    "Invalid lower bound"
  )
  df9 <- df %>%
    dplyr::mutate(upr_ = 2)
  expect_error(
    sens(df9, group, obs, p, 1, 0.5, eta = 0.01, upr_col = upr_, n_threads = 1, chunk_size = 200),
    "Invalid upper bound"
  )

  # Check that epsilon is positive, length one, and less than one.
  expect_error(
    sens(df, group, obs, p, 1, -1, eta = 0.01, n_threads = 1, chunk_size = 200),
    "must be a single, positive numeric"
  )
  expect_error(
    sens(df, group, obs, p, 1, c(1, 2), eta = 0.01, n_threads = 1, chunk_size = 200),
    "must be a single, positive numeric"
  )
  expect_warning(
    sens(df, group, obs, p, 1, 1.1, eta = 0.01, n_threads = 1, chunk_size = 200),
    "greater than 1"
  )

  # Check that eta is strictly positive, length one, and not too large
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0, n_threads = 1, chunk_size = 200),
    "must be a single, positive numeric"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = c(1, 2), n_threads = 1, chunk_size = 200),
    "must be a single, positive numeric"
  )
  expect_warning(
    sens(df, group, obs, p, 1, 0.5, eta = 2, n_threads = 1, chunk_size = 200),
    "large relative to",
  )

  # Check that m is strictly positive and length one
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, m = 0, n_threads = 1, chunk_size = 200),
    "must be a single positive integer"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, m = c(1, 2), n_threads = 1, chunk_size = 200),
    "must be a single positive integer"
  )

  # Check that n_threads is strictly positive, length one, and roundable to
  # integer, and not larger than the number of cores
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 0, chunk_size = 200),
    "single positive integer"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = c(1, 2), chunk_size = 200),
    "single positive integer"
  )
  expect_warning(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1.1, chunk_size = 200),
    "not a whole number"
  )
  # NOTE: To avoid using too many cores on CRAN, we only check that the number
  # of threads is not larger than the number of cores when not on CRAN.
  if (!testthat:::on_cran()) {
    expect_warning(
      sens(df, group, obs, p, 1, 0.5,
        eta = 0.01, n_threads =
          parallel::detectCores() + 1, chunk_size = 200
      ),
      "greater than the number of cores"
    )
  }

  # Check that chunk size is strictly positive, length one, and roundable to
  # integer
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 0),
    "single positive integer"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = c(1, 2)),
    "single positive integer"
  )
  expect_warning(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 1.1),
    "not a whole number"
  )

  # Check that `N` is non-negative, length one, and roundable to integer
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200, N = -1),
    "single non-negative integer"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200, N = c(1, 2)),
    "single non-negative integer"
  )
  expect_warning(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200, N = 0.1),
    "not a whole number"
  )

  # Check that alpha is between 0 and 1, length one, greater than one when multiplied by N, and
  # greater than 5 when multiplied by N
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200, alpha = -1),
    "must be a single numeric"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200, alpha = c(1, 2)),
    "must be a single numeric"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200, alpha = 1.1),
    "between zero and one"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200, alpha = 0),
    "between zero and one"
  )
  expect_error(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200, N = 1, alpha = 1 / 2),
    "too small"
  )
  expect_warning(
    sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200, N = 4, alpha = 1 / 2),
    "small relative to"
  )
})

test_that("Errors when data integrity violated", {
  n <- 1000
  p <- runif(n)^2
  group <- sample(1:3, n, replace = TRUE)
  obs <- runif(n) <= p + (group == 2) * 0.1
  # Upr and lwr formed from logit of p
  logit <- function(x) log(x / (1 - x))
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  lwr <- inv_logit(logit(p) - 0.1)
  upr <- inv_logit(logit(p) + 0.1)
  df <- dplyr::tibble(p, group, obs, lwr, upr) %>%
    dplyr::mutate(group = factor(group, levels = 1:3))

  # Error if less than two groups
  df2 <- df %>%
    dplyr::mutate(group = 1)
  expect_error(
    sens(df2, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "two distinct groups"
  )

  # Warn if unused levels in group
  df3 <- df %>%
    dplyr::mutate(group = factor(group, levels = 1:4))
  expect_warning(
    sens(df3, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "Some groups had no observations"
  )

  # Error if no observations in a group or observation status
  df3 <- df %>%
    dplyr::filter(group != 1 | obs)
  expect_error(
    sens(df3, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "at least one observed and one unobserved"
  )
  df4 <- df %>%
    dplyr::filter(obs == FALSE)
  expect_error(
    sens(df4, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "at least one observed and one unobserved"
  )

  # Error if there are any missing values
  df5 <- df %>%
    dplyr::mutate(p = ifelse(runif(n) < 0.1, NA, p))
  expect_error(
    sens(df5, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 200),
    "can be missing"
  )

  # Error if not lwr <= p <= upr
  df6 <- df %>%
    dplyr::mutate(p = dplyr::if_else(runif(n) < 0.1, inv_logit(logit(lwr) - 0.1), p))
  expect_error(
    sens(df6, group, obs, p, 1, 0.5, eta = 0.01, lwr_col = lwr, upr_col = upr, n_threads = 1, chunk_size = 200),
    "risks .* are invalid"
  )
  df7 <- df %>%
    dplyr::mutate(p = dplyr::if_else(runif(n) < 0.1, inv_logit(logit(upr) + 0.1), p))
  expect_error(
    sens(df7, group, obs, p, 1, 0.5, eta = 0.01, lwr_col = lwr, upr_col = upr, n_threads = 1, chunk_size = 200),
    "risks .* are invalid"
  )

  # Error if the data are not sortable
  df8 <- df %>%
    dplyr::mutate(
      lwr = pmin(p, 1 - p),
      upr = pmax(p, 1 - p)
    )
  expect_error(
    sens(df8, group, obs, p, 1, 0.5, eta = 0.01, lwr_col = lwr, upr_col = upr, n_threads = 1, chunk_size = 200),
    "sorted"
  )
})

test_that("Bootstrapped confidence intervals are valid", {
  # This test is slow, so don't run it on CRAN
  skip_on_cran()
  # It's also probabilistic, and can occasionally fail, so don't run it on CI
  skip_on_ci()

  n <- 1000
  p <- runif(n)^2
  group <- sample(1:2, n, replace = TRUE)
  obs <- runif(n) <= p + (group == 2) * 0.1
  # Upr and lwr formed from logit of p
  logit <- function(x) log(x / (1 - x))
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  lwr <- inv_logit(logit(p) - 0.1)
  upr <- inv_logit(logit(p) + 0.1)
  df <- dplyr::tibble(p, group, obs, lwr, upr) %>%
    dplyr::mutate(group = factor(group, levels = 1:2))

  # Check that the confidence intervals are valid
  sens_fit <- sens(df, group, obs, p, 1, 0.5, eta = 0.01, n_threads = 1, chunk_size = 100, N = 1000)
  expect_true(all(sens_fit$beta_min_2_02.5 <= sens_fit$beta_min_2))
  expect_true(all(sens_fit$beta_min_2 <= sens_fit$beta_min_2_97.5))
  expect_true(all(sens_fit$beta_max_2_02.5 <= sens_fit$beta_max_2))
  expect_true(all(sens_fit$beta_max_2 <= sens_fit$beta_max_2_97.5))
})
