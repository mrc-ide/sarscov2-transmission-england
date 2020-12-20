prepare_trajectories <- function(sample, counterfactual) {
  ## Prepare data including calculating incidence flows
  calc_incid <- function(x) {
    apply(x, 2, diff)
  }

  i_nowcast <- seq_along(counterfactual$trajectories$date)
  nowcast <- sample$trajectories$state[, , i_nowcast]
  nowcast[is.na(nowcast)] <- 0
  nowcasts <- list(
    deaths_inc = calc_incid(t(nowcast["deaths", , ])),
    deaths_comm_inc = calc_incid(t(nowcast["deaths_comm", , ])),
    infections_inc = calc_incid(t(nowcast["infections", , ])),
    icu = nowcast["icu", -1, ],
    general = nowcast["general", -1, ])

  counterfactual <- counterfactual$trajectories$state
  counterfactual[is.na(counterfactual)] <- 0
  counterfactuals <- list(
    deaths_inc = calc_incid(t(counterfactual["deaths", , ])),
    deaths_comm_inc = calc_incid(t(counterfactual["deaths_comm", , ])),
    infections_inc = calc_incid(t(counterfactual["infections", , ])),
    icu = counterfactual["icu", -1, ],
    general = counterfactual["general", -1, ])

  ## extract quantiles
  extract_qs <- function(trajectories) {
    ps <- c(0.025, 0.25, 0.5, 0.75, 0.975)
    apply(trajectories, 1, quantile, ps, na.rm = TRUE)
  }

  list(nowcast = lapply(nowcasts, extract_qs),
       counterfactual = lapply(counterfactuals, extract_qs))
}


## TODO: remove default args
plot_lockdown_counterfactual <- function(sample, counterfactual, title,
                                         days_earlier, data,
                                         policy_date = as.Date("2020-03-23"),
                                         labels = c("Actual", "Alternative"),
                                         label_y = c(0.6,  0.95),
                                         label_pos = c(2, 3),
                                         death_col = "#5C5992FF",
                                         col,
                                         dcol = "black", alpha = 0.3,
                                         legend = FALSE, ysf = 100) {

  dx <- as.Date(data$full$date_string)

  qs <- prepare_trajectories(counterfactual,  sample = sample)

  policy_dates <- policy_date - c(0, days_earlier)

  ymax_deaths <- max(qs$counterfactual$deaths_inc)
  ymax_deaths <- ceiling(ymax_deaths / ysf) * ysf

  xlim <- c(as.Date("2020-02-01"), as.Date("2020-12-01"))
  ylim <- c(0, ymax_deaths)

  plot(xlim[1], 0, type = "n", xlim = xlim, ylim = ylim,
       yaxp = c(0, ymax_deaths, 2), xlab = "", ylab = "Daily deaths")

  add_trajectories(dx[-1], qs$counterfactual$deaths_inc[, -1], col, alpha)
  add_trajectories(dx[-1], qs$nowcast$deaths_inc[, -1], death_col, alpha)

  points(dx, data$full$deaths_hosp + data$full$deaths_comm, pch = 23,
         bg = death_col, col = dcol, cex = 0.5, lwd = 0.6)

  y1 <- label_y * ymax_deaths
  segments(x0 = policy_dates, y0 = 0, y1 = y1, col = c(death_col, col),
           lty = 2,
           lend = 1, lwd = 1.5)

  text(x = policy_dates, pos = label_pos, offset = 0.2, y = y1, cex = 1,
       col = grey(0.2), font = 3,
       labels = labels)

  title(main = title, adj = 0, font.main = 1, cex = 0.8, line = 1)

  if (legend) {
    plot.new()
    legend(x = "topleft", y.intersp = 1.5,
           legend = c("Fitted", labels),
           fill = c(death_col, cols), bty = "n", cex = 0.8)
  }
}


add_trajectories <- function(x, qs, col, alpha) {
  cols <- add_alpha(rep(col, 2), alpha)
  CI_bands(qs[c("2.5%", "25%", "75%", "97.5%"), ], x, cols = cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50%", ], col = col, lty = 1, lwd = 1.5, lend = 1)
}


plot_carehome_counterfactual <- function(sample, counterfactual, title, data,
                                         death_col = "#5C5992FF", col,
                                         dcol = "black", alpha = 0.3,
                                         legend = FALSE) {
  dx <- as.Date(data$full$date_string)
  dy <- data$full$deaths_comm
  qs <- prepare_trajectories(counterfactual, sample = sample)

  ymax_deaths <- max(qs$counterfactual$deaths_comm_inc)
  ymax_deaths <- ceiling(max(ymax_deaths, dy, na.rm = TRUE) / 50) * 50

  xlim <- c(as.Date("2020-03-01"), as.Date("2020-12-01"))
  ylim <- c(0, ymax_deaths)

  plot(xlim[1], 0, type = "n", xlim = xlim, ylim = ylim,
       yaxp = c(0, ymax_deaths, 2), xlab = "",
       ylab = "Daily care home deaths")

  add_trajectories(dx[-1], qs$counterfactual$deaths_comm_inc[, -1], col, alpha)
  add_trajectories(dx[-1], qs$nowcast$deaths_comm_inc[, -1], death_col, alpha)

  points(dx, dy, pch = 23, bg = death_col, col = dcol, cex = 0.5, lwd = 0.6)

  title(main = title, adj = 0, font.main = 1, cex = 0.8, line = 1)

  if (legend) {
    legend(x = "topright", y.intersp = 1.5,
           legend = c("Fitted", "Counterfactual"),
           fill = c(death_col, col), bty = "n", cex = 1)
  }
}
