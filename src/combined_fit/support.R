switch_levels <- function(x) {
  nms <- names(x[[1]])
  y <- lapply(nms, function(z) lapply(x, "[[", z))
  names(y) <- nms
  y
}


calc_ifr_t <- function(sample) {
  step <- sample$trajectories$step

  ifr_t <- function(i) {
    p <- sample$predict$transform(sample$pars[i, ])
    p_sympt <- p$p_sympt

    expand <- function(x) {
      sircovid::sircovid_parameters_beta_expand(step, x)
    }

    p_hosp_sympt <- outer(p$psi_hosp_sympt, expand(p$p_hosp_sympt_step))
    p_ICU_hosp <- outer(p$psi_ICU_hosp, expand(p$p_ICU_hosp_step))
    p_death_ICU <- outer(p$psi_death_ICU, expand(p$p_death_ICU_step))
    p_death_hosp_D <- outer(p$psi_death_hosp_D, expand(p$p_death_hosp_D_step))
    p_death_comm <- outer(p$psi_death_comm, expand(p$p_death_comm_step))
    p_death_stepdown <- outer(p$psi_death_stepdown,
                              expand(p$p_death_stepdown_step))

    IHR_by_age <- p_sympt * p_hosp_sympt * (1 - p_death_comm) * 100
    IFR_by_age <- p_sympt * p_hosp_sympt * p_death_comm * 100 +
      IHR_by_age * (p_ICU_hosp * (p_death_ICU + (1 - p_death_ICU) *
                                  p_death_stepdown) +
                    (1 - p_ICU_hosp) * p_death_hosp_D)

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

calc_cum_inf_age <- function(sample, data, a, date = NULL) {

  if (is.null(date)) {
    date <- tail(data$fitted$date_string, 1)
  }
  p <- sample$predict$transform(sample$pars[1, ])
  S0 <- p$N_tot[a]

  index_S <- grep("^S_", names(sample$predict$index))
  n_groups <- p$n_groups

  St <- apply(sample$trajectories$state[index_S[seq(a, length(index_S),
                                              n_groups)], , , drop = FALSE],
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

  ret$deaths <- pmax(ret$deaths,
                     ret$deaths_hosp + ret$deaths_comm,
                     na.rm = TRUE)
  ret[-seq_len(nrow(ret) - 5), ] <- NA
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

prepare_model_admissions <- function(model, dashboard){
  model <- mean_admissions

  age <- c("0 to 5", "6 to 17", "18 to 64", "65 to 84", "85+")

  ## Model outputs only for people in the general population
  model <- data.frame(model[1:17, ] / 100)

  ## Aggregate into 5 age brackets comparable with Dashboard data
  row1 <- model[1, ]
  row2 <- colSums(model[2:4, ])
  row3 <- colSums(model[5:13, ])
  row4 <- colSums(model[14:16, ])
  row5 <- model[17, ]
  model <- rbind(row1, row2, row3, row4, row5)
  model$age <- factor(age, levels = age)

  model <- data.frame(melt(model, id.vars = "age"))
  colnames(model) <- c("age", "region", "admissions_prop")
  model$source <- as.factor("Model")

  ## Prepare the dashboard data
  dashboard <- dashboard %>% dplyr::select(age, region, admissions_prop)
  dashboard$source <- as.factor("PHE Dashboard")

  out <- rbind(model, dashboard)
}

prepare_model_data <- function(admissions, admissions_data) {
  model <- admissions$england$mean_prop_total_admissions
  model <- data.frame(model[1:17] / 100)
  age <- c("0 to 4", "5 to 9", "10 to 14", "15 to 19", "20 to 24", "25 to 29",
           "30 to 34", "35 to 39", "40 to 44", "45 to 49", "50 to 54",
           "55 to 59", "60 to 64", "65 to 69", "70 to 74", "75 to 79", "80+")
  model$age <- age
  model$age <- factor(model$age, levels = age)

  model <- melt(model)
  model$variable <- "Model"
  out <- rbind(model, admissions_data)
}
