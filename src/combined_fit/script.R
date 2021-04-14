regions <- read.csv("regions.csv", stringsAsFactors = FALSE)
rownames(regions) <- regions$key

dashboard <- read.csv("dashboard_region.csv")
admissions_data <- read.csv("admissions_age_data.csv")

path <- sprintf("regional_results/sample_pmcmc_%s.rds", regions$key)
samples <- Map(readRDS, path)
names(samples) <- regions$key

path <- sprintf("regional_results/pmcmc_%s.rds", regions$key)
pmcmcs <- Map(readRDS, path)
names(pmcmcs) <- regions$key

path <- sprintf("regional_results/Rt_%s.rds", regions$key)
Rt_outputs <- Map(readRDS, path)
names(Rt_outputs) <- regions$key

path <- sprintf("regional_results/data_%s.rds", regions$key)
data <- Map(readRDS, path)
names(data) <- regions$key

path <- sprintf("regional_results/ifr_t_%s.rds", regions$key)
ifr_t <- Map(readRDS, path)
names(ifr_t) <- regions$key

path <- sprintf("regional_results/deaths_%s.rds", regions$key)
deaths <- Map(readRDS, path)
names(deaths) <- regions$key

path <- sprintf("regional_results/admissions_%s.rds", regions$key)
admissions <- Map(readRDS, path)
names(admissions) <- regions$key

## Aggregate region samples
region_names <- sircovid::regions("england")

agg_samples <- samples
agg_samples$england$trajectories <-
  sircovid::combine_trajectories(samples[region_names])

agg_regions <- rbind(regions,
                     england = c("England", "england", NA))
agg_data <- data
agg_data$england$full <- aggregate_data(agg_data, region_names, "full")
agg_data$england$fitted <- aggregate_data(agg_data, region_names, "fitted")
agg_data$england$fitted$deaths <- NA

## EpiEstim Rt calculation at national level directly
inc <- agg_samples$england$trajectories$state["infections_inc", , ]
step <- agg_samples$england$trajectories$step
## remove columns with NA incidence
idx_to_remove <- which(is.na(colSums(inc)))
 
## for EpiEstim we only used parameter values that are fixed across regions
## and not estimated so we can use this simple call to get what we want
p <- agg_samples$london$predict$transform(agg_samples$london$pars[1, ])
## p_C is age-varying, and will be weighted by the population in the GT
## calculation so we need the population for England
p$N_tot <- sircovid::carehomes_parameters(1, "england")$N_tot
set.seed(1) # to fix the GT distribution which is stochastically drawn
gt_distr <- sircovid:::gt_distr(p, n = 10000)
rt_EpiEstim <-
  carehomes_rt_trajectories_epiestim(step[-idx_to_remove],
                                     inc[, -idx_to_remove],
                                     p, gt_distr = gt_distr,
                                     sliding_window_ndays = 1)

i <- seq(0, length(gt_distr) - 1, 1)
mean_gt_distr <- sum(gt_distr * i) # 6.9
std_gt_distr <- sqrt(sum(gt_distr * (i - mean_gt_distr)^ 2)) # 3.7

## EpiEstim from deaths
daily_deaths <- agg_data$england$full$deaths
days_deaths <- agg_data$england$full$date_end
to_keep <- which(!is.na(daily_deaths))
daily_deaths <- daily_deaths[to_keep]
days_deaths <- days_deaths[to_keep]
## calculate this averaged across all regions, age groups, parameter sets
mean_days_infection_to_death <-
  mean(unlist(lapply(regions$key, function(region)
    mean(vapply(seq_len(nrow(agg_samples[[region]]$pars)),
                function(e) calculate_mean_time_from_infection_to_death(
                  agg_samples[[region]]$predict$transform(
                    agg_samples[[region]]$pars[e, ]),
                  death_hosp_only = FALSE), rep(0, 19))))))

## for now the below:
## 1) looks at values at the start of the pandemic only
## 2) takes a brutal mean of the delay across all age groups
## 3) only uses a single parameter set from a single region (1, London)
## this does not account for differential risk of infection
## or differential risk of death by age
## but since these are all quite similar (except for carehome deaths)...

mean_days_infection_to_death <- round(mean(mean_days_infection_to_death))


rt_EpiEstim_deaths_7 <-
  sircovid::carehomes_rt_trajectories_epiestim(days_deaths * 4,
                                               t(as.matrix(daily_deaths)),
                                               p, gt_distr = gt_distr,
                                               sliding_window_ndays = 7)


## Aggregate Rt outputs
agg_Rt_outputs <- Rt_outputs
agg_Rt_outputs$england <- sircovid::combine_rt(Rt_outputs[region_names],
                                               samples[region_names])

## IFR
ifr_t$england <- sircovid::combine_rt(ifr_t[region_names],
                                      samples[region_names])

## Summarise Rt outputs by region for ease of plotting
Rt_by_region <- switch_levels(agg_Rt_outputs)

## Parameter outputs
dir.create("outputs", FALSE, TRUE)
rmarkdown::render("paper_numbers.Rmd",
                  output_file = "outputs/paper_numbers.html")

## New model inputs for next run
dir.create("new_inputs", FALSE, TRUE)

## Covariance matrices
cov_mats <- read_and_combine_covmats(regions)
write.csv(cov_mats, "new_inputs/parameters_proposal.csv", row.names = FALSE)

## Optimal tuning for next run
tuned_inputs <- read.csv("pmcmc_inputs.csv", stringsAsFactors = FALSE)
tuned_inputs <- update_tuned_inputs(pmcmcs, tuned_inputs)
write.csv(tuned_inputs, "new_inputs/parameters_info.csv", row.names = FALSE)
date <- "2020-12-02"
## Create figures
dir.create("figs", FALSE, TRUE)
png("figs/forest_plot.png", width = 17, height = 20, res = 350, units = "cm",
    pointsize = 8)
par(bty = "n", mar = c(3, 0, 0.5, 1.5), mgp = c(2, 0.75, 0),
    oma = c(0, 0, 3, 1))
create_forest_plot(samples, date)
dev.off()

png("figs/data_fits_regional.png", width = 3360, height = 1800, res = 200)
par(mfcol = c(6, 7), oma = c(1, 1, 4, 1), mar = c(3, 3, 1, 1),
    par(mgp = c(1.7, 0.5, 0), bty = "n"))
mapply(plot_all_trajectories,
       sample = samples[region_names],
       data = data[region_names],
       MoreArgs = list(date = date,
                       col_nowcast = paper_cols["nowcast"],
                       col_data_fitted = paper_cols["data_fitted"],
                       col_data_unfitted = paper_cols["data_unfitted"],
                       xlim = as.Date(c("2020-03-16", date)),
                       what = c("deaths_hosp", "deaths_carehomes",
                                "deaths_comm", "icu", "general",
                                "all_admission")))
mtext(side = 3, text = regions[region_names, "name"],
      line = 0.5, outer = TRUE, cex = 0.8,
      at = c(0.08, 0.225, 0.365, 0.515, 0.65, 0.795, 0.94))
dev.off()


png("figs/incidence_per_1000.png", width = 2400, height = 1200, res = 200)
par(mfrow = c(2, 4), oma = c(2, 1, 2, 1), mar = c(3, 3, 3, 1),
    mgp = c(1.7, 0.5, 0), bty = "n")
mapply(plot_incidence_per_1000,
       samples = samples,
       title = regions$name,
       col = paper_cols["nowcast"],
       ymax = 10)
dev.off()

png("figs/Rt_eff_general.png", width = 2400, height = 1200, res = 200)
par(mfrow = c(2, 4), oma = c(2, 1, 2, 1), mar = c(3, 3, 3, 1))
mapply(plot_Rt,
       sample_Rt = Rt_by_region$eff_Rt_general,
       Rt_date = Rt_by_region$date,
       title = agg_regions$name,
       MoreArgs = list(date = date, col = paper_cols["nowcast"],
                       ylab = expression(R[e](t))))
dev.off()

png("figs/Rt_general.png", width = 2400, height = 1200, res = 200)
par(mfrow = c(2, 4), oma = c(2, 1, 2, 1), mar = c(3, 3, 3, 1))
mapply(plot_Rt,
       sample_Rt = Rt_by_region$Rt_general,
       Rt_date = Rt_by_region$date,
       title = agg_regions$name,
       MoreArgs = list(date = date, col = paper_cols["nowcast"],
                       ylab = expression(R[0](t))))
dev.off()

## Admissions by age
admissions_spline <- lapply(agg_samples, extract_admissions_by_age)
mean_admissions <- sapply(admissions_spline, "[[", "mean_prop_total_admissions")
write.csv(mean_admissions, "outputs/mean_admissions_by_age.csv")

plot_admissions_spline <- prepare_model_data(admissions_spline, admissions_data)
png("figs/admissions_spline.png", units = "in",
    width = 10, height = 5, res = 300)
ggplot(plot_admissions_spline, aes(x = age, y = value, colour = variable)) +
  geom_point() + xlab("") + ylab("Proportion") + theme_bw() + ylim(c(0, 0.5)) +
  labs(title = "All admissions by age (England)", colour = "Source") +
  ggsci::scale_color_lancet() +
  theme(axis.text.x = element_text(angle = 90),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()

## Plot outputs by age
png("figs/deaths_by_age.png", width = 2400, height = 1800, res = 200)
par(mfrow = c(7, 4), oma = c(1, 1, 4, 1), mar = c(3, 3, 0.5, 0.5))
mapply(FUN = plot_all_by_age,
       sample = deaths[region_names],
       MoreArgs = list(date = date,
                       what = c("death_0", "death_55", "death_65", "death_75"),
                       yield = "deaths"))
mtext(side = 3, text = c("Deaths under 55s", "Deaths 55 to 64",
                         "Deaths 65 to 74", "Deaths 75+"),
      line = 0.5, outer = TRUE, cex = 0.9, at = c(0.15, 0.395, 0.64, 0.89))
mtext(side = 2, text = regions[region_names, "name"],
      line = -0.5, outer = TRUE, cex = 0.7,
      at = c(0.1, 0.24, 0.38, 0.53, 0.67, 0.81, 0.96))
dev.off()


png("figs/admissions_by_age.png", width = 2400, height = 1800, res = 200)
par(mfrow = c(7, 4), oma = c(1, 1, 4, 1), mar = c(3, 3, 0.5, 0.5))
mapply(FUN = plot_all_by_age,
       sample = admissions[region_names],
       MoreArgs = list(date = date,
                       what = c("adm_0", "adm_18", "adm_65", "adm_85"),
                       yield = "admissions"))
mtext(side = 3, text = c("Admissions under 18s", "Admissions 18 to 64",
                         "Admissions 65 to 84", "Admissions 85+"),
      line = 0.5, outer = TRUE, cex = 0.9, at = c(0.15, 0.395, 0.64, 0.89))
mtext(side = 2, text = regions[region_names, "name"],
      line = -0.5, outer = TRUE, cex = 0.7,
      at = c(0.1, 0.24, 0.38, 0.53, 0.67, 0.81, 0.96))
dev.off()

png("figs/gt_distr.png", width = 1000, height = 1000, res = 200)
par(mar = c(5, 5, .5, .5))
plot(seq(0, length(gt_distr) - 1, 1), gt_distr, type = "h",
     xlab = "Days",
     ylab = "Probability mass function \nof the generation time")
dev.off()

png("figs/epiestim.png", width = 1000, height = 1500, res = 200)

Rt_ref <- agg_Rt_outputs$england
Rt_1 <- rt_EpiEstim
Rt_2 <- rt_EpiEstim_deaths_7
shift_t_2 <- mean_days_infection_to_death


plot_epiestim(Rt_ref, Rt_1, Rt_2, shift_t_2)

dev.off()



## Figures for the manuscript itself
dir.create("paper_figs", FALSE, TRUE)

# assign a unique colour to each region
region_cols <- paper_cols[seq_along(region_names)]
names(region_cols) <- region_names

png("paper_figs/fig1.png", width = 2400, height = 2000, res = 200,
    pointsize = 18)
par(bty = "n", mgp = c(1.7, 0.5, 0), oma = c(0, 1, 2, 0.5),
    mgp = c(1.5, 0.7, 0))
plot_paper_fig1(agg_samples, Rt_by_region$eff_Rt_general,
                agg_data, regions, region_names, region_cols,
                paper_cols[c("hosp", "carehomes", "comm")])
dev.off()

png("paper_figs/fig2.png", width = 2400, height = 2000, res = 200,
    pointsize = 18)
par(oma = c(1, 1, 1, 1), cex = 1, bty = "n")
plot_paper_fig2(agg_samples, agg_data, regions, region_names)
dev.off()

png("paper_figs/fig4.png", width = 17, height = 10, units = "cm", res = 400,
    pointsize = 9)
par(oma = c(1, 1, 1, 1), mar = c(3, 3, 3, 1),  bty = "n")
plot_paper_fig4(agg_samples, ifr_t, regions, region_names, region_cols)
dev.off()

maps <- readRDS("NHS_region_map.rds")
png("paper_figs/fig5.png", width = 2800, height = 1600, pointsize = 45)
par(oma = c(2, 1, 2, 1), bty = "n")
plot_paper_fig5(samples, data, date, region_names, maps, regions)
dev.off()

png("paper_figs/supplement_pillar2_react.png", width = 2400, height = 1000,
    res = 200)
par(mfcol = c(2, 7), oma = c(1, 1, 1, 1), mar = c(2, 3, 2, 0), bty = "n",
    mgp = c(1.5, 0.7, 0))
mapply(plot_supplement_testing,
       sample = samples[region_names],
       data = data[region_names],
       title = regions[region_names, "name"],
       MoreArgs = list(col = paper_cols["nowcast"],
                       col_data = grey(0.2)))
dev.off()
