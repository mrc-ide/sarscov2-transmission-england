regions <- read.csv("regions.csv", stringsAsFactors = FALSE)
rownames(regions) <- regions$key

path <- sprintf("inputs/data_%s.rds", regions$key)
data <- Map(readRDS, path)

path <- sprintf("inputs/results_%s.rds", regions$key)
counterfactuals_by_region <- Map(readRDS, path)

names(counterfactuals_by_region) <- names(data) <- regions$key
counterfactuals <- switch_levels(counterfactuals_by_region)

agg_counterfactuals <- lapply(counterfactuals,  function(x)
  list(trajectories = sircovid::combine_trajectories(x)))

agg_data <- data
agg_data$england <- list(full = aggregate_data(data))

message("Rendering paper numbers")
dir.create("outputs", FALSE, TRUE)
saveRDS(agg_counterfactuals, "outputs/counterfactuals.rds")
rmarkdown::render("paper_numbers.Rmd", output_dir = "outputs")

message("Creating plots")
dir.create("figs", FALSE, TRUE)

n_region <- length(counterfactuals$fitted)
col <- paper_cols["green"]
death_col <- grey(0.4)

png("figs/paper_fig.png", width = 2400, height = 1200, res = 200,
    pointsize = 18)
par(mfcol = c(2, 3), bty = "n", mgp = c(1.7, 0.7, 0),
    oma = c(1, 1, 1, 1), mar = c(2, 2.5, 2, 0))

counterfactual_lockdown <- c(
  lockdown_earlier = "A) Lockdown 1 week earlier",
  lockdown_later = "D) Lockdown 1 week later",
  open_later = "B) Relax lockdown 2 weeks later",
  open_earlier = "E) Relax lockdown 2 weeks earlier")
counterfactual_carehome <- c(
  chr_contact_reduce_0.5 = "C) 50% reduced care home visits",
  chr_contact_increase_0.5 = "F) 50% increased care home visits")

Map(plot_lockdown_counterfactual,
    counterfactual = agg_counterfactuals[names(counterfactual_lockdown)],
    title = unname(counterfactual_lockdown),
    col = col,
    days_earlier = c(7, -7, -14, 14),
    label_y = list(c(1, 0.6), c(0.6, 1), c(0.6, 1), c(0.95, 0.6)),
    label_pos = list(c(4, 3), c(2, 4), c(2, 4), c(4, 3)),
    MoreArgs = list(
        sample = agg_counterfactuals$fitted,
        legend = FALSE, data = agg_data$england,
        death_col = death_col,
        ysf = 2000))
Map(plot_carehome_counterfactual,
    counterfactual = agg_counterfactuals[names(counterfactual_carehome)],
    legend = c(TRUE, FALSE),
    title = unname(counterfactual_carehome),
    MoreArgs = list(
        data = agg_data$england,
        sample = agg_counterfactuals$fitted,
        death_col = death_col,
        col = col))
dev.off()

png("figs/supp_paper_fig.png", width = 2400, height = 2400, res = 200)
par(mfcol = c(7, 6), bty = "n", mgp = c(1.7, 0.5, 0),
    oma = c(1, 1, 1, 1), mar = c(3, 3, 2, 0.5))
Map(plot_lockdown_counterfactual, sample = counterfactuals$fitted,
    counterfactual = counterfactuals$lockdown_earlier,
    legend = FALSE,
    data = data, title = sprintf("A%d) %s", seq_len(n_region), regions$name),
    MoreArgs = list(days_earlier = 7, death_col = grey(0.7),
                    col = col, ysf = 1e2, label_pos = c(4, 3)))

Map(plot_lockdown_counterfactual, sample = counterfactuals$fitted,
    counterfactual = counterfactuals$lockdown_later,
    legend = FALSE,
    data = data, title = sprintf("B%d) %s", seq_len(n_region), regions$name),
    MoreArgs = list(days_earlier = -7, death_col = grey(0.7),
                    col = col, ysf = 1e2, label_pos = c(3, 4)))

Map(plot_lockdown_counterfactual, sample = counterfactuals$fitted,
    counterfactual = counterfactuals$open_earlier,
    legend = FALSE,
    data = data, title = sprintf("C%d) %s", seq_len(n_region), regions$name),
    MoreArgs = list(days_earlier = 14, death_col = grey(0.7),
                    col = col, ysf = 1e2, label_pos = c(4, 3)))
Map(plot_lockdown_counterfactual, sample = counterfactuals$fitted,
    counterfactual = counterfactuals$open_later,
    legend = FALSE,
    data = data, title = sprintf("D%d) %s", seq_len(n_region), regions$name),
    MoreArgs = list(days_earlier = -14, death_col = grey(0.7),
                    col = col, ysf = 1e2, label_pos = c(3, 4)))

Map(plot_carehome_counterfactual, sample = counterfactuals$fitted,
    counterfactual = counterfactuals$chr_contact_reduce_0.5,
    data = data, title = sprintf("E%d) %s", seq_len(n_region), regions$name),
    legend = FALSE,
    MoreArgs = list(death_col = grey(0.7),
                    col = col))
Map(plot_carehome_counterfactual, sample = counterfactuals$fitted,
    counterfactual = counterfactuals$chr_contact_increase_0.5,
    data = data, title = sprintf("F%d) %s", seq_len(n_region), regions$name),
    legend = FALSE,
    MoreArgs = list(death_col = grey(0.7),
                    col = col))
dev.off()
