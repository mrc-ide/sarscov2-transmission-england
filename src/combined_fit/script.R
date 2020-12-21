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

## Calculate IFR over time
ifr_t <- lapply(samples, calc_ifr_t)

## Aggregate region samples
region_names <- c("north_west",
                  "north_east_and_yorkshire",
                  "midlands",
                  "east_of_england",
                  "london",
                  "south_west",
                  "south_east")

agg_samples <- samples
agg_samples$england$trajectories <-
  sircovid::combine_trajectories(samples[region_names])

agg_regions <- rbind(regions,
                     england = c("England", "england", NA))
agg_data <- data
agg_data$england$full <- aggregate_data(agg_data, region_names, "full")
agg_data$england$fitted <- aggregate_data(agg_data, region_names, "fitted")
agg_data$england$fitted$deaths <- NA

## Aggregate Rt outputs
agg_Rt_outputs <- Rt_outputs
agg_Rt_outputs$england <- sircovid::combine_rt(Rt_outputs[region_names],
                                               samples[region_names])

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
png("figs/forest_plot.png", width = 2400, height = 1600, res = 200)
create_forest_plot(samples, date)
dev.off()

png("figs/data_fits_regional.png", width = 3360, height = 1800, res = 200)
par(mfcol = c(6, 7), oma = c(1, 1, 4, 1), mar = c(3, 3, 0.5, 0.5),
    par(mgp = c(1.7, 0.5, 0), bty = "n"))
mapply(plot_all_trajectories,
       sample = samples[region_names],
       data = data[region_names],
       MoreArgs = list(date = date,
                       col_nowcast = paper_cols["nowcast"],
                       col_data_fitted = paper_cols["data_fitted"],
                       col_data_unfitted = paper_cols["data_unfitted"],
                       xlim = as.Date(c("2020-03-16", date)),
                       what = c("deaths_hosp", "deaths_comm", "icu",
                                "general", "hosp", "new_admitted")))
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
admissions <- lapply(agg_samples, extract_admissions_by_age)
mean_admissions <- sapply(admissions, "[[", "mean_prop_total_admissions")

plot_admissions <- prepare_model_admissions(admissions, dashboard)
region_name_map <- setNames(c(regions$key, "england"),
                            c(regions$name, "England"))
plot_admissions$region <-
  names(region_name_map)[match(plot_admissions$region, region_name_map)]
p1 <- ggplot(plot_admissions,
             aes(x = age, y = admissions_prop, fill = region)) +
  geom_bar(position = "dodge", stat = "identity", aes(fill = source)) +
  theme_bw() + ggtitle("A) Admissions by age and region") +
  xlab("") + theme(legend.title = element_blank()) +
  ggsci::scale_fill_lancet() +
  ylab("Proportion") +
  facet_wrap(~region) +
  theme(axis.text.x = element_text(angle = 90))

plot_admissions_2 <- prepare_model_data(admissions, admissions_data)
p2 <- ggplot(plot_admissions_2, aes(x = age, y = value, colour = variable)) +
  geom_point() + xlab("") + ylab("Proportion") + theme_bw() + ylim(c(0, 0.5)) +
  labs(title = "B) All admissions by age (England)", colour = "Source") +
  ggsci::scale_color_lancet() + theme(axis.text.x = element_text(angle = 90))
png("figs/dashboard_model.png", units = "in", width = 10, height = 5,
    res = 300)
gridExtra::grid.arrange(p1, p2, nrow = 1)
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
                paper_cols[c("hosp", "comm")])
dev.off()

png("paper_figs/fig2.png", width = 2400, height = 2000, res = 200,
    pointsize = 18)
par(oma = c(1, 1, 1, 1), cex = 1, bty = "n")
plot_paper_fig2(agg_samples, agg_data, regions, region_names)
dev.off()

png("paper_figs/fig4.png", width = 2400, height = 1600, res = 200,
    pointsize = 18)
par(oma = c(1, 1, 1, 1), mar = c(3, 3, 3, 1),  bty = "n")
plot_paper_fig4(agg_samples, ifr_t, regions, region_names, region_cols)
dev.off()

maps <- readRDS("NHS_region_map.rds")
png("paper_figs/fig5.png", width = 2800, height = 1600, pointsize = 45)
par(oma = c(2, 2, 2, 2), bty = "n")
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
