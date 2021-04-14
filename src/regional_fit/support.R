calc_ifr_t <- function(sample) {
  step <- sample$trajectories$step

  ifr_t <- function(i) {
    p <- sample$predict$transform(sample$pars[i, ])
    p_C <- p$p_C

    expand <- function(x) {
      sircovid::sircovid_parameters_beta_expand(step, x)
    }

    p_H <- outer(p$psi_H, expand(p$p_H_step))
    p_ICU <- outer(p$psi_ICU, expand(p$p_ICU_step))
    p_ICU_D <- outer(p$psi_ICU_D, expand(p$p_ICU_D_step))
    p_H_D <- outer(p$psi_H_D, expand(p$p_H_D_step))
    p_G_D <- outer(p$psi_G_D, expand(p$p_G_D_step))
    p_W_D <- outer(p$psi_W_D, expand(p$p_W_D_step))

    IHR_by_age <- p_C * p_H * (1 - p_G_D) * 100
    IFR_by_age <- p_C * p_H * p_G_D * 100 +
      IHR_by_age * (p_ICU * (p_ICU_D + (1 - p_ICU_D) * p_W_D) +
                      (1 - p_ICU) * p_H_D)

    pop <- p$N_tot
    nms <- names(sample$predict$index)
    index_S <- grep("^S_", nms)

    j <- length(which(!sample$trajectories$predicted))
    S <- sample$trajectories$state[index_S, i, j]
    S_mat <- matrix(S, nrow = p$n_groups)
    S <- rowSums(S_mat) # sum across all vaccination states
    infected <- unname(pop - S)
    gen_pop <- c(pop[1:p$n_age_groups], 0, 0)
    gen_pop_infected <- c(infected[1:p$n_age_groups], 0, 0)

    IFR_age <- split(IFR_by_age, seq_len(nrow(IFR_by_age)))
    names(IFR_age) <- gsub("S", "IFR", nms[index_S][seq_len(p$n_groups)])
    IHR_age <- split(IHR_by_age, seq_len(nrow(IHR_by_age)))
    names(IHR_age) <- gsub("S", "IHR", nms[index_S][seq_len(p$n_groups)])

    ret <-
      list(IFR_pop_all = apply(IFR_by_age, 2, weighted.mean, pop),
           IFR_pop_gen = apply(IFR_by_age, 2, weighted.mean, gen_pop),
           IFR_AR_all  = apply(IFR_by_age, 2, weighted.mean, infected),
           IFR_AR_gen  = apply(IFR_by_age, 2, weighted.mean, gen_pop_infected),
           IHR_pop_all = apply(IHR_by_age, 2, weighted.mean, pop),
           IHR_pop_gen = apply(IHR_by_age, 2, weighted.mean, gen_pop),
           IHR_AR_all  = apply(IHR_by_age, 2, weighted.mean, infected),
           IHR_AR_gen  = apply(IHR_by_age, 2, weighted.mean, gen_pop_infected))

    c(ret, IFR_age, IHR_age)
  }

  res <- lapply(seq_len(dim(sample$pars)[1]), ifr_t)

  ## These are stored in a list-of-lists and we convert to a
  ## list-of-matrices here
  collect <- function(nm) {
    matrix(unlist(lapply(res, "[[", nm)), length(step), length(res))
  }

  nms <- names(res[[1]])
  ret <- lapply(nms, collect)
  names(ret) <- nms
  ret$date <- sample$trajectories$date

  ret
}


carehomes_spim_thin <- function(samples, control) {

  thin <- ceiling(control$pmcmc$n_chains *
                    (control$pmcmc$n_steps - control$burnin) /
                    control$n_sample)

  ## Note we add one to the burnin to account for the initial values
  ret <- mcstate::pmcmc_thin(samples, control$burnin + 1L, thin)

  incidence_states <- c("deaths", "deaths_hosp", "deaths_comm",
                        "deaths_carehomes", "admitted", "diagnoses",
                        "infections", "sympt_cases", "sympt_cases_over25")

  ret$trajectories$date <- ret$trajectories$step / ret$trajectories$rate

  ret$trajectories <- sircovid:::add_trajectory_incidence(
    ret$trajectories, incidence_states)
  keep <- !grepl("^(I_weighted|D_hosp|prob_strain|cum_n_vaccinated)_",
                 rownames(ret$trajectories$state))
  ret$trajectories$state <- ret$trajectories$state[keep, , ]

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


extract_outputs_by_age <- function(sample, what) {

  trajectories <- sample$trajectories$state
  cum_output <- trajectories[grep(paste0("^", what),
                                  rownames(trajectories)), , ]

  total_output <- cum_output[, , dim(cum_output)[3]]
  prop_output <- t(total_output) / colSums(total_output) * 100

  output <- apply(cum_output, 1:2, diff)
  mean_output <- apply(output, 1:2, mean)

  out <- list(prop_total_output = prop_output,
              mean_prop_total_output = colMeans(prop_output),
              output_t = mean_output,
              lower_bound = apply(output, 1:2, quantile, 0.025),
              upper_bound = apply(output, 1:2, quantile, 0.975))

  out <- aggregate_outputs_by_age(out, what)

  out
}

aggregate_outputs_by_age <- function(object, what) {

  which_df <- c("output_t", "lower_bound", "upper_bound")
  out <- NULL

  # Model outputs have 19 age groups: 17 5-year age groups up to 80 and 80+,
  # plus CHW and CHR
  # CHW are evenly distributed across 8 age bands between 25 and 65
  # CHR are split 0.05, 0.05, 0.15 and 0.75 between the 65-69, 70-74, 75-79 and
  # 80+ bands.

  # Let's split modelled trajectories to match age bands in data
  # Younger age groups will get aggregated, given overall low numbers
  if (what == "cum_admit") {

    for (df in which_df) {

      # Under 18s' admissions include all of bands 1-3, and 3/5th of 15-19s
      # 18-64s the remaining 2/5th, all of bands 5-13 and all CHW (band 18)
      # 65-84s bands 14-16, 40% (assumed) of 80+, and all CHR under 80 plus 40%
      # (assumed) of CHR 80+
      adm_0 <- rowSums(object[[df]][, 1:3], na.rm = TRUE) +
        (object[[df]][, 4] * 0.6)
      adm_18 <- rowSums(object[[df]][, 5:13], na.rm = TRUE) +
        (object[[df]][, 4] * 0.4) + object[[df]][, 18]
      adm_65 <- rowSums(object[[df]][, 14:16], na.rm = TRUE) +
        (object[[df]][, 17] * 0.4) +
        (object[[df]][, 19] * 0.55)
      adm_85 <- (object[[df]][, 17] * 0.6) + (object[[df]][, 19] * 0.45)

      out[[df]] <- cbind(adm_0, adm_18, adm_65, adm_85)
    }

  } else if (what == "D_hosp") {

    for (df in which_df) {

      # CHW and CHR deaths are assumed uniformly distributed as per their age
      # composition. CHW deaths are thus 6/8th in under 55s and 2/8th in 55-64s,
      # and CHR deaths are 10% in 65-74s and 90% in 75+s
      death_0 <- rowSums(object[[df]][, 1:11], na.rm = TRUE) +
        (object[[df]][, 18] * 0.75)
      death_55 <- rowSums(object[[df]][, 12:13], na.rm = TRUE) +
        (object[[df]][, 18] * 0.25)
      death_65 <- rowSums(object[[df]][, 14:15], na.rm = TRUE) +
        (object[[df]][, 19] * 0.1)
      death_75 <- rowSums(object[[df]][, 16:17], na.rm = TRUE) +
        (object[[df]][, 19] * 0.9)

      out[[df]] <- cbind(death_0, death_55, death_65, death_75)

    }
  }
  out
}

sort_data_admissions <- function(df, r) {

  df <- admissions_nhse
  r <- region

  vector_age_bands <- unique(df$age)
  age_bands <- c("date", "region", "adm_0", "adm_18", "adm_65", "adm_85")
  df <- df[df$region == r, ]

  df <- df %>%
    dplyr::group_by(date, age) %>%
    dplyr::mutate(age = paste0("adm_", sub("\\_.*", "", age))) %>%
    tidyr::pivot_wider(., names_from = age, values_from = count) %>%
    mutate(adm_0 = adm_0 + adm_6) %>%
    select(all_of(age_bands))

  df
}
