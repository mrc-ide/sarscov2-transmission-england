transform_beta_date <- function(sample, i, days_earlier, i_betas) {
  p <- sample$pars[i, ]
  ret <- sample$predict$transform(p)

  beta_date <- sircovid::sircovid_date(sample$info$beta_date)
  beta_date <- sircovid::sircovid_date(
    c("2020-03-16", "2020-03-23", "2020-03-25",
      "2020-05-11", "2020-06-15", "2020-07-04",
      "2020-08-03", "2020-09-01", "2020-09-14",
      "2020-10-14", "2020-11-05"))
  beta_date[i_betas] <- beta_date[i_betas] - days_earlier

  beta_value <- p[sprintf("beta%s", seq_along(beta_date))]

  ret$beta_step <-
    sircovid::sircovid_parameters_beta(beta_date, beta_value, ret$dt)

  ret
}


transform_par <- function(sample, i, par, sens) {
  p <- sample$pars[i, ]
  p[par] <- p[par] * sens
  sample$predict$transform(p)
}


## Ordinarily we could run this in parallel with dust::dust_simulate,
## a list of parameters and a matrix of starting positions. However,
## that needs a little update to cope with the staggered starting
## date; here we just run the samples in series instead.
run_counterfactual <- function(sample, transform) {
  p <- sample$pars
  res_it <- lapply(seq_len(nrow(p)), run_counterfactual1, sample, transform)
  trajectories <- list(
    state = join_arrays(res_it),
    date = sample$trajectories$date[!sample$trajectories$predicted])
  trajectories$predicted <- rep(FALSE, length(trajectories$date))
  list(trajectories = trajectories)
}


run_counterfactual1 <- function(i, sample, transform) {
  p_i <- transform(sample, i)
  model <- sample$predict$filter$model$new(p_i, 1, 1)
  initial <- sample$predict$filter$initial(model$info(), 1, p_i)
  model$reset(p_i, initial$step)
  model$set_state(initial$state, initial$step)
  index <- c(date = 1, sample$predict$filter$index(model$info())$state)
  step <- c(initial$step, sample$predict$filter$data$step_end)
  traj <- dust::dust_iterate(model, step, index)
  traj[!grepl("^(S|cum_admit)_", rownames(traj)), , ]
}


join_arrays <- function(x) {
  ret <- array(0, dim = c(nrow(x[[1]]), length(x), ncol(x[[1]])))

  for (i in seq_along(x)) {
    ret[, i, ] <- x[[i]]
  }
  rownames(ret) <- rownames(x[[1]])
  ret
}


thin_sample <- function(sample, n_par) {
  n_par <- min(n_par, nrow(sample$pars))
  i <- seq_len(n_par)

  sample$pars <- sample$pars[i, ]
  sample$probabilities <- sample$probabilities[i, ]
  sample$state <- sample$state[, i]
  sample$chain <- sample$chain[i]
  sample$iteration <- sample$iteration[i]
  sample$trajectories$state <- sample$trajectories$state[, i, ]

  sample
}


version_check <- function(package, version) {
  if (packageVersion(package) < version) {
    stop(sprintf(
      paste("Please update %s with: drat:::add('ncov-ic');",
            "install.packages('%s')"),
      package, package))
  }
}
