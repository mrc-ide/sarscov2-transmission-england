carehomes_spim_thin <- function(samples, control) {
  thin <- ceiling(control$n_chains * (control$n_mcmc - control$burnin) /
                    control$n_sample)

  ret <- mcstate::pmcmc_thin(samples, control$burnin, thin)

  incidence_states <- c("deaths", "deaths_hosp", "deaths_comm",
                        "admitted", "new", "infections",
                        "sympt_cases", "sympt_cases_over25")

  ret$trajectories$date <- ret$trajectories$step / ret$trajectories$rate

  ret$trajectories <- sircovid:::add_trajectory_incidence(
    ret$trajectories, incidence_states)

  ret
}


calculate_spim_Rt <- function(samples) {
  step <- samples$trajectories$step

  index_S <- grep("^S_", names(samples$predict$index))
  S <- samples$trajectories$state[index_S, , , drop = FALSE]

  pars <- lapply(seq_len(nrow(samples$pars)), function(i)
    samples$predict$transform(samples$pars[i, ]))

  ret <- sircovid::carehomes_Rt_trajectories(
    step, S, pars,
    initial_step_from_parameters = TRUE,
    shared_parameters = FALSE)

  ret
}


check_rows_data_rtm <- function(data_rtm) {
  ## Check the number of rows for each date
  rows_out <- data_rtm %>% group_by(date) %>% summarise(rows = n())
  all_regions <- unique(data_rtm$region)

  dates_incomplete <- which(rows_out$rows < length(all_regions))

  if (length(dates_incomplete > 0)) {
    for (d in rows_out$date[dates_incomplete]) {
      i <- !all_regions %in% data_rtm[data_rtm$date == d, "region"]
      missing_regions <- all_regions[i]

      tmp <-  data_rtm %>% filter(date == d)
      tmp_add <- tmp[seq_along(missing_regions), ]
      tmp_add[3:ncol(tmp_add)] <- NA
      tmp_add$region <- missing_regions
      tmp_add$date <- d
      data_rtm <- data_rtm %>% bind_rows(tmp_add)
    }
  }

  data_rtm
}


sircovid_info <- function(samples, region, beta_date) {
  data <- samples$predict$transform(samples$pars[1, ])
  info <- samples$predict$filter$model$new(data, 0, 1)$info()
  list(version = packageVersion("sircovid"),
       info = info,
       data = data,
       region = region,
       beta_date = beta_date)
}
