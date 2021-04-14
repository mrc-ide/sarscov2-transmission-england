
carehomes_spim_data <- function(date, rtm, serology, region,
                                trim_deaths, steps_per_day,
                                full_data = FALSE) {
  rtm <- carehomes_spim_data_rtm(date, rtm, region, full_data)
  serology <- carehomes_spim_data_serology(date, serology, region)

  ## Merge the two datasets on date
  stopifnot(all(serology$date %in% rtm$date))
  i <- match(rtm$date, serology$date)
  serology <- serology[i, ]
  rownames(serology) <- NULL
  data <- cbind(rtm, serology[setdiff(names(serology), "date")])

  ## We need a dummy row here so that the sampler stops before the
  ## first day in order to get the differences in cumulative deaths
  ## correct.
  data <- rbind(data[1, ], data)
  data[1, ] <- NA
  data$date[[1]] <- data$date[[2]] - 1

  ## At this point we'll save our "real" date into date_string and
  ## work with "sircovid" dates which are just integer dates into 2020
  data$date_string <- data$date
  data$date <- sircovid::sircovid_date(data$date)

  ## Set last 'trim_deaths' days with deaths reported to NA, as these
  ## are too likely to be back-filled to be reliable
  i <- seq(to = nrow(data), length.out = trim_deaths)
  data[i, c("deaths", "deaths_non_hosp", "deaths_hosp",
            "deaths_comm", "deaths_carehomes")] <- NA

  mcstate::particle_filter_data(data, "date", steps_per_day, 1)
}

carehomes_spim_data_rtm <- function(date, data, region, full_data = FALSE) {
  vars <- c("phe_patients", "phe_occupied_mv_beds", "phe_admissions",
            "ons_death_hospital", "ons_death_carehome", "ons_death_community",
            "react_positive", "react_samples",
            "pillar2_positives_over25",
            "pillar2_negatives_total_pcr_over25",
            "pillar2_positives_lft_over25")
  data <- data[c("region", "date", vars)]

  ## Remove any data after the date of analysis parameter
  data <- data[as.Date(data$date) <= as.Date(date), ]

  ## Make sure the dates for each region match up
  data <- check_rows_data_rtm(data)

  data <- data[data$region == region, ]
  data$deaths <- NA_integer_
  data$deaths_hosp <- data$ons_death_hospital
  data$deaths_carehomes <- data$ons_death_carehome
  data$deaths_comm <- data$ons_death_community
  data$deaths_non_hosp <- NA_integer_

  ## Use non lateral flow pillar 2 negatives where we have them; data
  ## is not available in the same format for all regions, hence the
  ## need for these checks
  neg_over25_non_lfts_available <-
    !all(is.na(data$pillar2_negatives_total_pcr_over25))
  pos_over25_lfts_available <- !all(is.na(data$pillar2_positives_lft_over25))

  if (neg_over25_non_lfts_available && pos_over25_lfts_available) {
    data$pillar2_negatives_over25 <- data$pillar2_negatives_total_pcr_over25
    data$pillar2_positives_over25 <-
      ifelse(is.na(data$pillar2_positives_over25), 0,
             data$pillar2_positives_over25) -
      ifelse(is.na(data$pillar2_positives_lft_over25), 0,
             data$pillar2_positives_lft_over25)
  }

  ## Ignore pillar 2 testing before 2020-06-01
  drop <- c("pillar2_positives_over25", "pillar2_negatives_over25")
  data[which(data$date < "2020-06-01"), drop] <- NA_integer_

  ## Remove implausible value for MV beds occupancy in east_of_england
  ## on 2020-09-11
  if (region == "east_of_england") {
    data[which(data$phe_patients - data$phe_occupied_mv_beds < 0),
         "phe_occupied_mv_beds"] <- NA_integer_
  }

  ## Check no negative integers have been introduced into pillar 2 timeseries
  stopifnot(all(data$pillar2_negatives_over25 >= 0, na.rm = TRUE),
            all(data$pillar2_positives_over25 >= 0, na.rm = TRUE))

  ## Retain variables in names to be used by sircovid
  ret <- data.frame(
    date = sircovid::as_date(data$date),
    deaths_hosp = data$deaths_hosp,
    deaths_carehomes = data$deaths_carehomes,
    deaths_comm = data$deaths_comm,
    deaths_non_hosp = data$deaths_non_hosp,
    icu = data$phe_occupied_mv_beds,
    general = data$phe_patients - data$phe_occupied_mv_beds,
    hosp = data$phe_patients,
    deaths = data$deaths,
    admitted = NA_integer_,
    diagnoses = NA_integer_,
    all_admission = data$phe_admissions,
    pillar2_tot = NA_integer_,
    pillar2_pos = NA_integer_,
    pillar2_cases = NA_integer_,
    pillar2_over25_tot = data$pillar2_positives_over25 +
      data$pillar2_negatives_over25,
    pillar2_over25_pos = data$pillar2_positives_over25,
    pillar2_over25_cases = data$pillar2_positives_over25,
    react_pos = data$react_positive,
    react_tot = data$react_samples,
    strain_non_variant = NA_integer_,
    strain_tot = NA_integer_,
    stringsAsFactors = FALSE)

  if (!full_data) {
    omit <- c("deaths", "deaths_non_hosp", "hosp", "admitted",
              "diagnoses", "pillar2_tot", "pillar2_pos",
              "pillar2_cases", "pillar2_over25_cases")
    for (i in omit) {
      ret[[i]] <- NA_integer_
    }
  }

  ret[which(ret$general < 0)] <- 0

  stopifnot(min(ret$general, na.rm = TRUE) >= 0)

  ret
}


carehomes_spim_data_serology <- function(date, data, region) {
  data <- data[data$region == region & data$age_group == "15_64", ]

  ## Remove any data after the date parameter
  data <- data[as.Date(data$date) <= as.Date(date), ]

  data.frame(date = sircovid::as_date(data$date),
             npos_15_64 = data$n_positive,
             ntot_15_64 = data$total_samples,
             stringsAsFactors = FALSE)
}


carehomes_spim_parameters <- function(info, prior, proposal, region) {
  info <- info[info$region == region & info$include, ]
  info <- info[setdiff(names(info), c("region", "include"))]
  rownames(info) <- NULL

  prior <- prior[prior$region == region & prior$name %in% info$name, ]
  prior <- prior[setdiff(names(prior), "region")]
  rownames(prior) <- NULL


  proposal <- proposal[proposal$region == region &
                         proposal$name %in% info$name, ]
  stopifnot(!any(duplicated(proposal$name)))
  stopifnot(all(info$name %in% proposal$name),
            all(info$name %in% names(proposal)))
  proposal <- as.matrix(proposal[info$name])
  rownames(proposal) <- colnames(proposal) <- info$name

  list(info = info, prior = prior, proposal = proposal, region = region)
}


carehomes_spim_transform <- function(beta_date, region,
                                     support_severity,
                                     support_progression) {
  force(region)

  beta_date <- sircovid::sircovid_date(beta_date)

  severity <- sircovid::sircovid_parameters_severity(support_severity)

  progression <- sircovid::carehomes_parameters_progression()
  progression[support_progression$parameter] <- support_progression$value

  function(pars) {
    start_date <- pars[["start_date"]]
    p_G_D <- pars[["p_G_D"]] # death in *community*
    p_G_D_CHR <- pars[["p_G_D_CHR"]] # death in *care home*
    eps <- pars[["eps"]]
    m_CHW <- pars[["m_CHW"]]
    m_CHR <- pars[["m_CHR"]]
    p_ICU <- pars[["p_ICU"]]
    p_ICU_D <- pars[["p_ICU_D"]]
    p_H_D <- pars[["p_H_D"]]
    p_W_D <- pars[["p_W_D"]]
    mu_D <- pars[["mu_D"]]
    mu_ICU <- pars[["mu_ICU"]]
    p_H <- pars[["p_H"]]
    p_H_CHR <- pars[["p_H_CHR"]]
    beta_value <- unname(pars[c("beta1", "beta2", "beta3", "beta4", "beta5",
                                "beta6", "beta7", "beta8", "beta9", "beta10",
                                "beta11", "beta12")])
    p_NC <- pars[["p_NC"]]
    rho_pillar2_tests <- pars[["rho_pillar2_tests"]]
    alpha_D <- pars[["alpha_D"]]
    alpha_H <- pars[["alpha_H"]]
    ## Total: 30 parameters

    ## Set age dependent severity based on Bob's analysis.
    ##
    ## For the parameters that we are varying (p_ICU,
    ## p_ICU_D, p_H_D and p_H) we scale the given
    ## parameters by their maximum value, then multiply by the
    ## parameter being varied.
    normalise_and_scale <- function(x, v) {
      v * x / max(x)
    }

    severity$p_ICU <-
      normalise_and_scale(severity$p_ICU, p_ICU)
    severity$p_H <-
      normalise_and_scale(severity$p_H, p_H)
    severity$p_sero_pos[] <- 0.85

    progression$gamma_sero_pre_1 <- 1 / 13
    progression$gamma_sero_pre_2 <- 1 / 13

    ret <- sircovid::carehomes_parameters(
      ## Core
      start_date = start_date,
      region = region,
      beta_date = beta_date,
      beta_value = beta_value,
      initial_I = 30,
      ## Severity
      severity = severity,
      sero_specificity = 0.99,
      sero_sensitivity = 1,
      p_death_carehome = p_G_D_CHR,
      ## Progression
      progression = progression,
      ## Transmission
      eps = eps,
      m_CHW = m_CHW,
      m_CHR = m_CHR,
      pillar2_specificity = 1,
      pillar2_sensitivity = 1,
      react_specificity = 1,
      react_sensitivity = 1,
      p_NC = p_NC)

    ret$I_A_transmission <- 0.223

    date_mu_change <- c("2020-04-01", "2020-06-01")
    p_ICU_D_dates <- sircovid::sircovid_date(date_mu_change)
    p_ICU_D_value <- c(p_ICU_D, p_ICU_D * mu_D)
    ret$p_ICU_D_step <- sircovid::sircovid_parameters_beta(
      p_ICU_D_dates, p_ICU_D_value, ret$dt)
    p_H_D_dates <- sircovid::sircovid_date(date_mu_change)
    p_H_D_value <- c(p_H_D, p_H_D * mu_D)
    ret$p_H_D_step <- sircovid::sircovid_parameters_beta(
      p_H_D_dates, p_H_D_value, ret$dt)
    p_W_D_dates <- sircovid::sircovid_date(date_mu_change)
    p_W_D_value <- c(p_W_D, p_W_D * mu_D)
    ret$p_W_D_step <- sircovid::sircovid_parameters_beta(
      p_W_D_dates, p_W_D_value, ret$dt)
    p_ICU_dates <- sircovid::sircovid_date(date_mu_change)
    p_ICU_value <- c(p_ICU, p_ICU * mu_ICU)
    ret$p_ICU_step <- sircovid::sircovid_parameters_beta(
      p_ICU_dates, p_ICU_value, ret$dt)
    p_star_dates <- sircovid::sircovid_date(c("2020-03-15", "2020-07-01",
                                              "2020-09-20", "2020-12-02"))
    p_star_value <- c(0.1, 0.42, 0.2, 0.25)
    ret$p_star_step <- sircovid::sircovid_parameters_beta(
      p_star_dates, p_star_value, ret$dt)

    ret$p_G_D_step <- p_G_D
    ret$psi_G_D[1:18] <- 1

    ## work out scaling for CHR for p_H and p_G_D
    ret$psi_G_D[19] <- p_G_D_CHR / p_G_D
    ret$psi_H[19] <- p_H_CHR / p_H

    ## Gamma values from fitted Erlang to PCR positivity
    ret$k_PCR_pre <- 1
    ret$gamma_PCR_pre <- 0.1922243
    ret$k_PCR_pos <- 1
    ret$gamma_PCR_pos <- 0.0597321
    ret$k_sero_pos <- 1
    ret$gamma_sero_pos <- 1 / 180

    ## Time to diagnosis if admitted without test
    ret$gamma_U <- 1 / 3

    ret$rho_pillar2_tests <- rho_pillar2_tests
    ret$phi_pillar2_cases <- 0.5

    ## kappa for hospital data streams (not all will actually be used)
    ret$kappa_ICU <- 1 / alpha_H
    ret$kappa_general <- 1 / alpha_H
    ret$kappa_hosp <- 1 / alpha_H
    ret$kappa_admitted <- 1 / alpha_H
    ret$kappa_diagnoses <- 1 / alpha_H
    ret$kappa_all_admission <- 1 / alpha_H
    ret$kappa_death_hosp <- 1 / alpha_H

    ## kappa for other death data streams (not all will actually be used)
    ret$kappa_death_carehomes <- 1 / alpha_D
    ret$kappa_death_comm <- 1 / alpha_D
    ret$kappa_death_non_hosp <- 1 / alpha_D
    ret$kappa_death <- 1 / alpha_D

    ret
  }
}


carehomes_spim_pmcmc_parameters <- function(parameters, beta_date,
                                            support_severity,
                                            support_progression) {
  stopifnot(
    identical(parameters$info$name, parameters$prior$name),
    identical(parameters$info$name, rownames(parameters$proposal)))

  prior <- lapply(split(parameters$prior, parameters$prior$name),
                  make_prior)

  pars <- Map(
    mcstate::pmcmc_parameter,
    name = parameters$info$name,
    initial = parameters$info$initial,
    min = parameters$info$min,
    max = parameters$info$max,
    discrete = parameters$info$discrete,
    prior = prior)

  transform <- carehomes_spim_transform(beta_date, region,
                                        support_severity,
                                        support_progression)

  mcstate::pmcmc_parameters$new(pars, parameters$proposal, transform)
}


make_prior <- function(d) {
  if (d$type == "gamma") {
    shape <- d$gamma_shape
    scale <- d$gamma_scale
    function(p) {
      dgamma(p, shape = shape, scale = scale, log = TRUE)
    }
  } else if (d$type == "beta") {
    shape1 <- d$beta_shape1
    shape2 <- d$beta_shape2
    function(p) {
      dbeta(p, shape1 = shape1, shape2 = shape2, log = TRUE)
    }
  } else if (d$type == "null") {
    NULL
  } else {
    stop("Unknown prior type")
  }
}


calculate_proposal_kernel <- function(samples, region) {
  m <- cov(samples$pars)
  rownames(m) <- NULL
  data.frame(region = region,
             name = colnames(m),
             m,
             stringsAsFactors = FALSE)
}


rtm_parallel_control <- function(n_chains) {
  nt <- as.integer(Sys.getenv("CONTEXT_CORES", Sys.getenv("MC_CORES", 1)))
  pos <- 1:4
  nw <- max(pos[nt %% pos == 0 & pos <= n_chains])
  message(sprintf("Running on %d workers with %d threads", nw, nt))
  list(n_threads_total = nt, n_workers = nw)
}


carehomes_pmcmc_control <- function(short_run, n_chains) {
  if (short_run) {
    n_mcmc <- 20
  } else {
    n_mcmc <- 21000
  }
  parallel <- rtm_parallel_control(n_chains)
  rerun_every <- 100
  mcstate::pmcmc_control(n_mcmc,
                         n_chains = n_chains,
                         n_threads_total = parallel$n_threads_total,
                         n_workers = parallel$n_workers,
                         rerun_every = rerun_every,
                         save_state = TRUE,
                         save_trajectories = TRUE,
                         progress = TRUE)
}


carehomes_control <- function(short_run, n_chains) {
  if (short_run) {
    n_particles <- 10
    n_sample <- 10
    burnin <- 1
  } else {
    n_particles <- 384
    n_sample <- 1000
    burnin <- 1000
  }
  list(pmcmc = carehomes_pmcmc_control(short_run, n_chains),
       n_particles = n_particles,
       n_sample = n_sample,
       burnin = burnin,
       trim_deaths = 0,
       forecast_days = 0)
}
