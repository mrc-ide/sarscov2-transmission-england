source("util.R")

version_check("sircovid", "0.7.4")
version_check("mcstate", "0.2.15")
version_check("dust", "0.5.3")

control <- carehomes_spim_control(short_run, chains)

## Trim deaths for x number of time steps, if potential backfill issues
trim_deaths <- 0

## The last date of data to use
date <- "2020-12-02"

## Dates are as follows
##  1. 2020-03-16 - PM advises WFH, against non-essential travel etc
##  2. 2020-03-23 - PM announces full lockdown
##  3. 2020-03-25 - Lockdown into full effect
##  4. 2020-05-11 - Initial easing of lockdown
##  5. 2020-06-15 - Non-essential shops can open
##  6. 2020-07-04 - Restaurants, pubs etc can open
##  7. 2020-08-03 - "Eat out to help out" scheme starts
##  8. 2020-09-01 - Schools reopen
##  9. 2020-09-14 - "Rule of six" introduced
## 10. 2020-10-14 - Tiered system introduced
## 11. 2020-10-31 - Lockdown 2 announced
## 12. 2020-11-05 - Lockdown 2 begins
beta_date <- c("2020-03-16", "2020-03-23", "2020-03-25",
               "2020-05-11", "2020-06-15", "2020-07-04",
               "2020-08-03", "2020-09-01", "2020-09-14",
               "2020-10-14", "2020-10-31", "2020-11-05")

message("Reading data")
data_rtm <- read_csv("data_rtm.csv")
data_serology <- read_csv("data_serology.csv")
parameters_info <- read_csv("parameters_info.csv")
parameters_prior <- read_csv("parameters_prior.csv")
parameters_proposal <- read_csv("parameters_proposal.csv")
support_progression <- read_csv("support_progression.csv")
support_severity <- read_csv("support_severity.csv")

message("Preparing parameters")
parameters <- carehomes_spim_parameters(parameters_info, parameters_prior,
                                        parameters_proposal, region)

if (kernel_scaling != 1) {
  ## Used when tuning kernels only
  message("  - scaling proposal kernel by ", kernel_scaling)
  parameters$proposal <- parameters$proposal * kernel_scaling
}

pars <- carehomes_spim_pmcmc_parameters(parameters, beta_date,
                                        support_severity, support_progression)

## Perturb initial parameter values across the chains
initial <- replicate(control$n_chains,
                     pars$propose(pars$initial(), 1))

steps_per_day <- 1 / pars$model(pars$initial())$dt
data <- carehomes_spim_data(date, data_rtm, data_serology, region,
                            trim_deaths, steps_per_day, full_data = FALSE)

filter <- sircovid::carehomes_particle_filter(data, control$n_particles,
                                              control$n_threads)

message("Running chains - this will take a while!")
samples <- mcstate::pmcmc(pars, filter, n_steps = control$n_mcmc,
                          n_chains = control$n_chains,
                          initial = initial,
                          save_trajectories = TRUE, progress = TRUE)
samples$info <- sircovid_info(samples, region, beta_date)

message("Remove burnin and thin PMCMC output")
samples_thin <- carehomes_spim_thin(samples, control)

message("Computing Rt")
rt <- calculate_spim_Rt(samples_thin)

## This is only useful if rerunning later
message("Computing new proposal kernel")
proposal <- calculate_proposal_kernel(samples, region)

message("Writing outputs")
samples[c("state", "trajectories", "predict")] <- list(NULL)
dir.create("outputs", FALSE, TRUE)
saveRDS(samples, "outputs/pmcmc_results.rds")
saveRDS(samples_thin, "outputs/sample_pmcmc_results.rds")
saveRDS(rt, "outputs/Rt.rds")
data <- list(fitted = data,
             full = carehomes_spim_data(date, data_rtm, data_serology, region,
                                        trim_deaths, steps_per_day,
                                        full_data = TRUE))
saveRDS(data, "outputs/data.rds")
write.csv(proposal, "outputs/parameters_covariance.csv", row.names = FALSE)
