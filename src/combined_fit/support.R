switch_levels <- function(x) {
  nms <- names(x[[1]])
  y <- lapply(nms, function(z) lapply(x, "[[", z))
  names(y) <- nms
  y
}


calc_cum_inf_age <- function(sample, data, a, date = NULL) {

  if (is.null(date)) {
    date <- tail(data$fitted$date_string, 1)
  }
  p <- sample$predict$transform(sample$pars[1, ])
  S0 <- p$N_tot[a]

  index_S <- grep("^S_", names(sample$predict$index))
  n_groups <- p$n_groups

  St <- apply(sample$trajectories$state[index_S[seq(a, length(index_S),
                                                    n_groups)], , ,
                                        drop = FALSE],
              c(2, 3), sum)
  cum_inf <- (S0 - St) / S0
  sircovid_date <- as.character(sircovid::sircovid_date(date))
  cum_inf[, which(sample$trajectories$date == sircovid_date)]

}


calc_inits_mle <- function(pmcmcs) {

  find_mle <- function(pmcmc, region_name) {
    w <- which.max(pmcmc$probabilities[, "log_posterior"])
    mle <- pmcmc$pars[w, ]
    par_names <- names(mle)
    data.frame(region = region_name,
               name = names(mle),
               initial = unname(mle),
               stringsAsFactors = FALSE)
  }

  res <- do.call(rbind, mapply(FUN = find_mle, pmcmc = pmcmcs,
                               region_name = names(pmcmcs), SIMPLIFY = FALSE))
  rownames(res) <- NULL

  res
}

update_tuned_inputs <- function(pmcmcs, tuned_inputs) {

  inits <- calc_inits_mle(pmcmcs)
  k <- sapply(seq_len(nrow(inits)),
              function(i) which(inits[i, "region"] == tuned_inputs$region
                                & inits[i, "name"] == tuned_inputs$name))
  tuned_inputs$initial[k] <- inits$initial
  tuned_inputs

}


read_and_combine_covmats <- function(regions) {
  path <- sprintf("regional_results/parameters_covariance_%s.csv", regions$key)
  res <- Map(read.csv, path)
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  res
}


aggregate_data <- function(data, region_names, type = "full") {

  d <- lapply(data[region_names], "[[", type)
  x <- d[[1]]

  date_names <- grep("date", names(x), value = TRUE)
  data_names <- setdiff(names(x), date_names)

  dim_t <- min(vnapply(d, nrow))

  ret <- Reduce("+", lapply(d, function(x) x[seq_len(dim_t), data_names]))

  death_nms <- c("deaths_hosp", "deaths_carehomes", "deaths_comm")
  ret$deaths <- pmax(ret$deaths,
                     rowSums(ret[, death_nms], na.rm = TRUE),
                     na.rm = TRUE)

  data.frame(ret, x[, date_names])
}

extract_admissions_by_age <- function(sample) {
  trajectories <- sample$trajectories$state
  cum_admissions <- trajectories[grep("^cum_admit", rownames(trajectories)), , ]

  total_admissions <- cum_admissions[, , dim(cum_admissions)[3]]
  prop_admissions <- t(total_admissions) / colSums(total_admissions) * 100

  admissions <- apply(cum_admissions, 1:2, diff)
  mean_admissions <- apply(admissions, 1:2, mean)

  list(prop_total_admissions = prop_admissions,
       mean_prop_total_admissions = colMeans(prop_admissions),
       mean_admissions_t = mean_admissions)
}

prepare_model_data <- function(admissions, admissions_data) {
  model <- admissions$england$mean_prop_total_admissions

  ## Add CHW and CHR admissions into age groups based on their age weightings
  model[6:13] <- model[6:13] + model[18] / 8
  model[14:17] <- model[14:17] + c(0.05, 0.05, 0.15, 0.75) * model[19]

  model <- data.frame(model[1:17] / 100)
  age <- c("0 to 4", "5 to 9", "10 to 14", "15 to 19", "20 to 24", "25 to 29",
           "30 to 34", "35 to 39", "40 to 44", "45 to 49", "50 to 54",
           "55 to 59", "60 to 64", "65 to 69", "70 to 74", "75 to 79", "80+")
  model$age <- age
  model$age <- factor(model$age, levels = age)

  model <- melt(model)
  model$variable <- "Model"

  rbind(model, admissions_data)
}

calculate_mean_time_from_infection_to_death <- function(p, death_hosp_only) {

  mean_discr_gamma <- function(k, gamma, dt) {
    k * dt / (1 - exp(-gamma * dt))
  }
  mean_time_to_hosp <- mean_discr_gamma(p$k_E, p$gamma_E, p$dt) +
    mean_discr_gamma(p$k_P, p$gamma_P, p$dt) +
    mean_discr_gamma(p$k_C_1, p$gamma_C_1, p$dt) +
    mean_discr_gamma(p$k_C_2, p$gamma_C_2, p$dt)

  mean_G_D <- mean_discr_gamma(p$k_G_D, p$gamma_G_D, p$dt)
  p_G_D <- p$p_G_D

  mean_H_D <- mean_discr_gamma(p$k_H_D, p$gamma_H_D, p$dt)
  p_H_D <- (1 - p$p_G_D) * (1 - p$p_ICU) * p$p_H_D

  mean_W_D <- mean_discr_gamma(p$k_ICU_pre, p$gamma_ICU_pre, p$dt) +
    mean_discr_gamma(p$k_ICU_W_D, p$gamma_ICU_W_D, p$dt) +
    mean_discr_gamma(p$k_W_D, p$gamma_W_D, p$dt)
  p_W_D <- (1 - p$p_G_D) * p$p_ICU * (1 - p$p_ICU_D) * p$p_W_D

  mean_ICU_D <- mean_discr_gamma(p$k_ICU_pre, p$gamma_ICU_pre, p$dt) +
    mean_discr_gamma(p$k_ICU_D, p$gamma_ICU_D, p$dt)
  p_ICU_D <- (1 - p$p_G_D) * p$p_ICU * p$p_ICU_D

  if (death_hosp_only) {
    mean_days_infection_to_death <- mean_time_to_hosp +
      (mean_H_D * p_H_D +
         mean_W_D * p_W_D +
         mean_ICU_D * p_ICU_D) / (p_H_D + p_W_D + p_ICU_D)
  }else {
    mean_days_infection_to_death <- mean_time_to_hosp +
      (mean_G_D * p_G_D +
         mean_H_D * p_H_D +
         mean_W_D * p_W_D +
         mean_ICU_D * p_ICU_D) / (p_G_D + p_H_D + p_W_D + p_ICU_D)
  }

  ## Issues: p_ICU, p_H_D, p_ICU_D, p_W_D vary over time

  mean_days_infection_to_death
}
