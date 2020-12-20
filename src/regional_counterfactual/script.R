version_check("sircovid", "0.7.4")
version_check("mcstate", "0.2.15")
version_check("dust", "0.5.3")

sample <- readRDS("inputs/sample_pmcmc.rds")
sample <- thin_sample(sample, n_par)

lockdown_days <- 7
open_days <- 14
eps <- c(0.25, 0.5, 0.75, 1)
lockdown_betas <- 1:3

message("Running lockdown 1 counterfactuals")
lockdown_earlier <-
  run_counterfactual(sample, transform = function(sample, i) {
    transform_beta_date(sample, i, lockdown_days, lockdown_betas)})

lockdown_later <-
  run_counterfactual(sample, transform = function(sample, i) {
    transform_beta_date(sample, i, -lockdown_days, lockdown_betas)})

open_earlier <-
  run_counterfactual(sample, transform = function(sample, i) {
    transform_beta_date(sample, i, open_days, -lockdown_betas)})

open_later <-
  run_counterfactual(sample, transform = function(sample, i) {
    transform_beta_date(sample, i, -open_days, -lockdown_betas)})

message("Running carehomes counterfactuals")
chr_contact_reduce <- lapply(eps, function(e)
  run_counterfactual(sample, transform = function(sample, i)
    transform_par(sample, i, "eps", 1 - e)))

chr_contact_increase <- lapply(eps, function(e)
  run_counterfactual(sample, transform = function(sample, i)
    transform_par(sample, i, "eps", 1 + e)))

names(chr_contact_reduce) <- sprintf("chr_contact_reduce_%s", eps)
names(chr_contact_increase) <- sprintf("chr_contact_increase_%s", eps)

## The sample objects will be larger than needed if we retain
## information needed to run further scenarios, so remove that here to
## save space.
sample$predict <- NULL
traj_names <- intersect(rownames(lockdown_earlier$trajectories$state),
                        rownames(sample$trajectories$state))
sample$trajectories$state <- sample$trajectories$state[traj_names, , ]

message("Saving results")
results <- c(list(fitted = sample,
                  lockdown_earlier = lockdown_earlier,
                  lockdown_later = lockdown_later,
                  open_earlier = open_earlier,
                  open_later = open_later),
             chr_contact_reduce,
             chr_contact_increase)
dir.create("outputs", FALSE, TRUE)
saveRDS(results, "outputs/results.rds")
