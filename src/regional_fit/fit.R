
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
  data[i, c("deaths", "deaths_hosp", "deaths_comm")] <- NA

  mcstate::particle_filter_data(data, "date", steps_per_day, 1)
}

carehomes_spim_data_rtm <- function(date, data, region, full_data = FALSE) {
  vars <- c("phe_patients", "phe_occupied_mv_beds",
            "admitted", "new", "phe_admissions",
            "death2", "death3", "death2_long", "death3_long",
            "pillar2_positives", "pillar2_negatives",
            "positives", "negatives", "react_positive", "react_samples",
            "pillar2_positives_over25", "pillar2_negatives_over25",
            "positives_over25", "pillar2_negatives_non_lft",
            "pillar2_negatives_non_lft_over25",
            "pillar2_positives_lft_over25", "pillar2_positives_lft")
  data <- data[c("region", "date", vars)]

  ## Remove any data after the date of analysis parameter
  data <- data[as.Date(data$date) <= as.Date(date), ]

  ## Make sure the dates for each region match up
  data <- check_rows_data_rtm(data)

  ## Set NA deaths to 0
  data[which(is.na(data$death2)), "death2"] <- 0
  data[which(is.na(data$death2_long)), "death2_long"] <- 0
  data[which(is.na(data$death3)), "death3"] <- 0
  data[which(is.na(data$death3_long)), "death3_long"] <- 0

  data <- data[data$region == region, ]
  data$deaths <- NA_integer_
  data$deaths_hosp <- data$death3
  data$deaths_comm <- data$death2 - data$death3

  ## Use non lateral flow pillar 2 negatives where we have them; data
  ## is not available in the same format for all regions, hence the
  ## need for these checks
  negatives_non_lfts_available <- !all(is.na(data$pillar2_negatives_non_lft))
  positives_lfts_available <- !all(is.na(data$pillar2_positives_lft))
  neg_over25_non_lfts_available <-
    !all(is.na(data$pillar2_negatives_non_lft_over25))
  pos_over25_lfts_available <- !all(is.na(data$pillar2_positives_lft_over25))

  if (negatives_non_lfts_available && positives_lfts_available) {
    data$pillar2_negatives <- data$pillar2_negatives_non_lft
    data$pillar2_positives <- ifelse(
      is.na(data$pillar2_positives), 0, data$pillar2_positives) -
      ifelse(is.na(data$pillar2_positives_lft), 0, data$pillar2_positives_lft)
  }
  if (neg_over25_non_lfts_available && pos_over25_lfts_available) {
    data$pillar2_negatives_over25 <- data$pillar2_negatives_non_lft_over25
    data$pillar2_positives_over25 <-
      ifelse(is.na(data$pillar2_positives_over25), 0,
             data$pillar2_positives_over25) -
      ifelse(is.na(data$pillar2_positives_lft_over25), 0,
             data$pillar2_positives_lft_over25)
  }

  ## Ignore pillar 2 testing before 2020-06-01
  drop <- c("pillar2_positives", "pillar2_negatives",
            "pillar2_positives_over25", "pillar2_negatives_over25")
  data[which(data$date < "2020-06-01"), drop] <- NA_integer_
  
  ## Remove implausible value for MV beds occupancy in east_of_england
  ## on 2020-09-11
  if (region == "east_of_england") {
    data[which(data$phe_patients - data$phe_occupied_mv_beds < 0),
         "phe_occupied_mv_beds"] <- NA_integer_
  }

  ## Check no negative integers have been introduced into pillar 2 timeseries
  stopifnot(all(data$pillar2_negatives >= 0, na.rm = TRUE),
            all(data$pillar2_positives >= 0, na.rm = TRUE),
            all(data$pillar2_negatives_over25 >= 0, na.rm = TRUE),
            all(data$pillar2_positives_over25 >= 0, na.rm = TRUE))

  ## Retain variables in names to be used by sircovid
  ret <- data.frame(
    date = sircovid::as_date(data$date),
    deaths_hosp = data$deaths_hosp,
    deaths_comm = data$deaths_comm,
    icu = data$phe_occupied_mv_beds,
    general = data$phe_patients - data$phe_occupied_mv_beds,
    hosp = data$phe_patients,
    deaths = data$deaths,
    admitted = data$admitted,
    new = data$new,
    new_admitted = data$phe_admissions,
    pillar2_tot = data$pillar2_positives + data$pillar2_negatives,
    pillar2_pos = data$pillar2_positives,
    pillar2_cases = data$pillar2_positives,
    pillar2_over25_tot = data$pillar2_positives_over25 +
      data$pillar2_negatives_over25,
    pillar2_over25_pos = data$pillar2_positives_over25,
    pillar2_over25_cases = data$pillar2_positives_over25,
    react_pos = data$react_positive,
    react_tot = data$react_samples,
    stringsAsFactors = FALSE)

  if (!full_data) {
    omit <- c("hosp", "admitted", "new", "pillar2_tot", "pillar2_pos",
              "pillar2_cases", "pillar2_over25_cases")
    for (i in omit) {
      ret[[i]] <- NA_integer_
    }
    if (all(is.na(ret$pillar2_over25_tot))) {
      ret$pillar2_tot <- data$pillar2_positives + data$pillar2_negatives
      ret$pillar2_pos <- data$pillar2_positives
      ret$pillar2_over25_tot <- NA_integer_
      ret$pillar2_over25_pos <- NA_integer_
    }
  }

  ret[which(ret$general < 0)] <- 0

  stopifnot(min(ret$deaths_comm, na.rm = TRUE) >= 0)
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
    p_death_comm <- pars[["p_death_comm"]] # death in *carehome*
    eps <- pars[["eps"]]
    C_1 <- pars[["C_1"]]
    C_2 <- pars[["C_2"]]
    p_ICU_hosp <- pars[["p_ICU_hosp"]]
    p_death_ICU <- pars[["p_death_ICU"]]
    p_death_hosp_D <- pars[["p_death_hosp_D"]]
    p_death_stepdown <- pars[["p_death_stepdown"]]
    mu_death_hosp <- pars[["mu_death_hosp"]]
    mu_ICU_hosp <- pars[["mu_ICU_hosp"]]
    p_hosp_sympt <- pars[["p_hosp_sympt"]]
    beta_value <- unname(pars[c("beta1", "beta2", "beta3", "beta4", "beta5",
                                "beta6", "beta7", "beta8", "beta9", "beta10",
                                "beta11", "beta12")])
    prop_noncovid_sympt <- pars[["prop_noncovid_sympt"]]
    rho_pillar2_tests <- pars[["rho_pillar2_tests"]]
    ## Total: 26 parameters

    ## Set age dependent severity based on Bob's analysis.
    ##
    ## For the parameters that we are varying (p_ICU_hosp,
    ## p_death_ICU, p_death_hosp_D and p_hosp_sympt) we scale the given
    ## parameters by their maximum value, then multiply by the
    ## parameter being varied.
    normalise_and_scale <- function(x, v) {
      v * x / max(x)
    }

    severity$p_ICU_hosp <-
      normalise_and_scale(severity$p_ICU_hosp, p_ICU_hosp)
    severity$p_hosp_sympt <-
      normalise_and_scale(severity$p_hosp_sympt, p_hosp_sympt)
    severity$p_seroconversion[] <- 0.85

    progression$gamma_R_pre_1 <- 1 / 13
    progression$gamma_R_pre_2 <- 1 / 13

    ret <- sircovid::carehomes_parameters(
      ## Core
      start_date = start_date,
      region = region,
      beta_date = beta_date,
      beta_value = beta_value,
      ## Severity
      severity = severity,
      sero_specificity = 0.99,
      sero_sensitivity = 1,
      p_death_carehome = p_death_comm,
      ## Progression
      progression = progression,
      ## Transmission
      eps = eps,
      C_1 = C_1,
      C_2 = C_2,
      pillar2_specificity = 1,
      pillar2_sensitivity = 1,
      react_specificity = 1,
      react_sensitivity = 1,
      prop_noncovid_sympt = prop_noncovid_sympt)

    date_mu_change <- c("2020-04-01", "2020-06-01")
    p_death_ICU_dates <- sircovid::sircovid_date(date_mu_change)
    p_death_ICU_value <- c(p_death_ICU, p_death_ICU * mu_death_hosp)
    ret$p_death_ICU_step <- sircovid::sircovid_parameters_beta(
      p_death_ICU_dates, p_death_ICU_value, ret$dt)
    p_death_hosp_D_dates <- sircovid::sircovid_date(date_mu_change)
    p_death_hosp_D_value <- c(p_death_hosp_D, p_death_hosp_D * mu_death_hosp)
    ret$p_death_hosp_D_step <- sircovid::sircovid_parameters_beta(
      p_death_hosp_D_dates, p_death_hosp_D_value, ret$dt)
    p_death_stepdown_dates <- sircovid::sircovid_date(date_mu_change)
    p_death_stepdown_value <- c(p_death_stepdown,
                                p_death_stepdown * mu_death_hosp)
    ret$p_death_stepdown_step <- sircovid::sircovid_parameters_beta(
      p_death_stepdown_dates, p_death_stepdown_value, ret$dt)
    p_ICU_hosp_dates <- sircovid::sircovid_date(date_mu_change)
    p_ICU_hosp_value <- c(p_ICU_hosp, p_ICU_hosp * mu_ICU_hosp)
    ret$p_ICU_hosp_step <- sircovid::sircovid_parameters_beta(
      p_ICU_hosp_dates, p_ICU_hosp_value, ret$dt)
    p_admit_conf_dates <- sircovid::sircovid_date(c("2020-03-15", "2020-05-01"))
    p_admit_conf_value <- c(0.1, 0.2)
    ret$p_admit_conf_step <- sircovid::sircovid_parameters_beta(
      p_admit_conf_dates, p_admit_conf_value, ret$dt)

    ## Gamma values from fitted Erlang to PCR positivity
    ret$s_PCR_pre <- 1
    ret$gamma_PCR_pre <- 0.1922243
    ret$s_PCR_pos <- 1
    ret$gamma_PCR_pos <- 0.0597321
    ret$s_R_pos <- 1
    ret$gamma_R_pos <- 0

    ## Setting infectiousness of hospitalised and care home cases to zero
    ret$ICU_transmission <- 0
    ret$hosp_transmission <- 0
    ret$comm_D_transmission <- 0

    ## Time to diagnosis if admitted without test
    ret$gamma_test <- 1 / 3

    ret$observation$rho_pillar2_tests <- rho_pillar2_tests
    ret$observation$phi_pillar2_cases <- 0.5
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


calc_serology_prior_means <- function(prior, region) {

  prior <- prior[prior$region == region, ]

  gamma_mean <- function(par) {
    prior[prior$name == par, "gamma_scale"] *
      prior[prior$name == par, "gamma_shape"]
  }
  beta_mean <- function(par) {
    prior[prior$name == par, "beta_shape1"] /
      (prior[prior$name == par, "beta_shape1"] +
         prior[prior$name == par, "beta_shape2"])
  }

  list(gamma_R_pre_1 = gamma_mean("gamma_R_pre_1"),
       gamma_R_pre_2 = gamma_mean("gamma_R_pre_2"),
       p_seroconversion = beta_mean("p_seroconversion"),
       sero_specificity = beta_mean("sero_specificity")
  )
}


carehomes_spim_control <- function(short_run, n_chains) {
  control <- list(
    n_particles = 96, # multiple of 32, so 32, 64, 96, 128, ...
    n_mcmc = 11000,
    n_sample = 1000,
    n_chains = n_chains,
    n_threads = get_n_threads(),
    burnin = 1000,
    forecast_days = 0)

  if (short_run) {
    control$n_particles <- 10
    control$n_mcmc <- 20
    control$n_sample <- 10
    control$burnin <- 1
  }

  control
}
