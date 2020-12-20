plot_all_trajectories <- function(sample, data, date, what, col_nowcast,
                                  col_data_fitted, col_data_unfitted,
                                  xlim, alpha = 0.3) {
  lapply(what, plot_trajectories, sample, data, date, col_nowcast,
         col_data_fitted, col_data_unfitted, xlim, alpha)
}


plot_trajectories <- function(what,
                              sample,
                              data,
                              date,
                              col_nowcast,
                              col_data_fitted,
                              col_data_unfitted,
                              xlim = as.Date(c("2020-02-01", "2020-12-01")),
                              alpha = 0.3) {

  trajnames <- c(deaths = "deaths_inc", deaths_comm = "deaths_comm_inc",
                 deaths_hosp = "deaths_hosp_inc", icu = "icu", hosp = "hosp",
                 general = "general", new = "new_inc",
                 admitted = "admitted_inc", new_admitted = "new_admitted_inc")

  labs <- c(deaths = "Daily deaths",
            deaths_comm = "Daily care home deaths",
            deaths_hosp = "Daily hospital deaths",
            icu = "ICU beds",
            general = "general beds",
            hosp = "Hospital beds",
            new = "Daily inpatient diagnoses",
            admitted = "Daily admissions",
            new_admitted = "Daily admissions (all)")

  dates <- sample$trajectories$date
  i_date <- dates <= sircovid::sircovid_date(date)
  ## Remove first date as the first time period is always much more
  ## than one day before data is recorded.
  i_date[1] <- FALSE
  x <- sircovid::sircovid_date_as_date(dates[i_date])
  state <- sample$trajectories$state[, , i_date]

  ## Extract trajectory of interest
  if (trajnames[what] == "new_admitted_inc") {
    y <- state["new_inc", , ] + state["admitted_inc", , ]
  } else {
    y <- state[trajnames[what], , ]
  }

  ## Calculate trajectory quantiles
  ps <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  qs <- apply(y, 2, quantile, ps, na.rm = TRUE)

  ## Extract data
  dx <- as.Date(data$fitted$date_string)
  dy <- data$fitted[, what]

  ## Plot data points in green if not used for fitting
  if (all(is.na(dy))) {
    dy <- data$full[, what]
    col_data <- col_data_unfitted
  } else {
    col_data <- col_data_fitted
  }

  ylim <- c(0, max(dy, y, na.rm = TRUE))
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       xlab = "",
       ylab = labs[what])

  CI_bands(qs[c("2.5%", "25%", "75%", "97.5%"), ], x,
           cols = add_alpha(rep(col_nowcast, 2), alpha),
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50%", ], col = col_nowcast, lty = 1,
        lwd = 1.5, lend = 1)

  points(dx, dy, pch = 23, bg = col_data, col = grey(0.2), cex = 0.7,
         lwd = 0.6)
}


plot_incidence_per_1000 <- function(samples, title, col, ymax, alpha = 0.3) {
  ## Remove first date as the first time period is always much more
  ## than one day before data is recorded.
  res <- samples$trajectories$state["infections_inc", , -1L]
  res <- res / sum(samples$predict$transform(samples$pars[1, ])$N_tot) * 1000
  x <- sircovid::sircovid_date_as_date(samples$trajectories$date[-1L])

  ps <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  qs <- apply(res, 2, quantile, ps, na.rm = TRUE)

  xlim <- range(x)
  ylim <- c(0, ymax)
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       main = toupper(title),
       font.main = 1,
       xlab = "", ylab = "Incidence per day per 1000")

  CI_bands(qs[c("2.5%", "25%", "75%", "97.5%"), ], x,
           cols = add_alpha(rep(col, 2), alpha),
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50%", ], col = col, lty = 1, lwd = 1.5)
}


plot_Rt <- function(sample_Rt,
                    Rt_date,
                    col,
                    date,
                    ylab,
                    title = "",
                    alpha = 0.3) {
  ## Drop the first value at the moment as it's likely the first time
  ## period is much more than one day
  sample_Rt <- sample_Rt[-1L, ]
  Rt_date <- Rt_date[-1L, ]

  x <- sircovid::sircovid_date_as_date(Rt_date[, 1])

  ps <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  qs <- apply(sample_Rt, 1, quantile, ps, na.rm = TRUE)

  par(mgp = c(1.7, 0.5, 0), bty = "n")
  xlim <- c(x[1], as.Date(date))
  ylim <- c(0, 4)
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       main = toupper(title),
       font.main = 1,
       xlab = "", ylab = ylab)

  cols <- add_alpha(rep(col, 2), alpha)

  CI_bands(qs[c("2.5%", "25%", "75%", "97.5%"), ], x, cols = cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50%", ], col = col, lty = 1, lwd = 1.5, lend = 1)
  abline(h = 1, lty = 2)

}


plot_paper_fig1 <- function(agg_samples, Rt_general, agg_data, regions,
                            region_names, col_region, col_fit) {
  samples <- agg_samples[region_names]
  Rt_general <- Rt_general[region_names]
  labels <- regions[region_names, "label"]
  region_cols <- region_cols[region_names]
  names(labels) <- region_names

  layout_mat <- matrix(c(1, 2, 3, 4,
                         5, 6, 7, 8,
                         9, 9, 9, 9),
                       nrow = 3, byrow = TRUE)

  layout(mat = layout_mat, heights = c(5, 5, 10), widths = c(1, 1, 1, 1))

  ## Plot A,
  par(mar = c(2.5, 0.5, 0.5, 0.5))
  ordered_regions <- start_date_forest_plot(samples, region_cols, labels)
  ordered_regions <- rev(ordered_regions)
  title("A) Epidemic start date", adj = 0, font.main = 1, line = 0, xpd = NA,
        cex.main = 1)


  ## Plot B - H,
  par(mar = c(2.5, 1.5, 0.5, 0))
  main <- sprintf("%s) %s (%s)",
                  LETTERS[seq_along(region_names) + 1L],
                  regions[region_names, "name"],
                  regions[region_names, "label"])
  plot_carehome_vs_hosp_deaths(agg_samples[region_names],
                               agg_data[region_names],
                               regions[region_names, ],
                               main = main,
                               col = col_fit,
                               ymax = 250,
                               surtitle = "Daily COVID-19 deaths",
                               ylab = "")

  par(mar = c(2.5, 1.5, 1, 0))
  plot_annotated_rt(samples, Rt_general, col = region_cols, regions = regions)
  title(expression("I) Transmission over time:"~italic({R[t]^eff})), adj = 0,
        font.main = 1, line = 0.5, xpd = NA)
}

plot_carehome_vs_hosp_deaths <- function(samples, data, regions, col, ymax,
                                         main, ylab, surtitle) {

  ## Plot legend at the end
  legend <- c(TRUE, rep(FALSE, length(samples) - 1))
  plot_surtitle <- legend

  Map(plot_carehome_vs_hosp_deaths1, samples, data, main, ymax, legend,
      plot_surtitle,
      MoreArgs = list(col = col, surtitle = surtitle, ylab = ylab))
}


start_date_forest_plot <- function(samples, cols, labels) {
  ## Order by start_date
  mean_start_date <- vapply(samples, function(x)
    mean(x$pars[, "start_date"]), numeric(1))

  names(cols) <- names(samples)
  ordered_regions <- names(sort(mean_start_date, decreasing = TRUE))
  cols <- cols[ordered_regions]
  samples <- samples[ordered_regions]
  labels <- labels[ordered_regions]

  ## Extract start_date
  numeric_start_date <- lapply(samples, function(x) {
    as.numeric(sircovid::sircovid_date_as_date(x$pars[, "start_date"]))
  })

  n_regions <- length(samples)
  ylim <- c(0.5, n_regions + 0.5)
  col_line <- grey(0.2)
  at <- seq_len(n_regions)

  ## start_date plot
  xlim <- as.Date(c("2020-01-01", "2020-03-01"))
  plot(x = as.Date("2020-01-01"),
       y = 1,
       type = "n",
       ylab = "",
       xlab = "",
       xlim = xlim,
       ylim = ylim,
       axes = FALSE)
  axis.Date(1, at = as.Date(c("2020-01-01", "2020-02-01", "2020-03-01")),
            format = "%b")

  mapply(FUN = plot_CI_bar, res = numeric_start_date, at = at, width = 0.15,
         col = cols, cex = 1.5)
  text(x = sapply(numeric_start_date, quantile, 0.975) + 5, y = at,
       labels = labels, font = 3)
}


plot_carehome_vs_hosp_deaths1 <- function(sample, data, col, main = "",
                                          ymax = 500, legend = FALSE,
                                          plot_surtitle = FALSE,
                                          surtitle = "", title_col = "black",
                                          max_date = as.Date("2020-12-01"),
                                          alpha = 0.3,
                                          ylab = "Daily COVID-19 deaths") {
  dx <- as.Date(data$fitted$date_string)

  extract_qs <- function(trajectories, what) {
    res <- t(trajectories$state[what, , ])
    x <- sircovid::sircovid_date_as_date(trajectories$date)
    rownames(res) <- as.character(x)
    res <- res[x <= max_date, ]
    ps <- c(0.025, 0.25, 0.5, 0.75, 0.975)
    apply(res,  MARGIN = 1, FUN = quantile, ps, na.rm = TRUE)
  }

  qs_death_comm <- extract_qs(sample$trajectories, what = "deaths_comm_inc")
  qs_death_hosp <- extract_qs(sample$trajectories, what = "deaths_hosp_inc")

  add_trajectories <- function(qs, col) {
    cols <- add_alpha(rep(col, 2), alpha)
    x <- as.Date(colnames(qs))
    CI_bands(qs[c("2.5%", "25%", "75%", "97.5%"), ], x, cols = cols,
             horiz = FALSE, leg = FALSE)
    lines(x, qs["50%", ], col = col, lty = 1, lwd = 1.5, lend = 1)
  }

  xlim <- c(as.Date("2020-03-01"), max_date)
  ylim <- c(0, ymax)
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       las = 1,
       yaxp = c(0, ymax, 2),
       xlab = "", ylab = ylab)
  add_trajectories(qs = qs_death_comm, col = col["comm"])
  add_trajectories(qs = qs_death_hosp, col = col["hosp"])
  points(dx, data$fitted[, "deaths_hosp"], pch = 23, bg = col["hosp"],
         cex = 0.5, lwd = 0.6)
  points(dx, data$fitted[, "deaths_comm"], pch = 23, bg = col["comm"],
         cex = 0.5, lwd = 0.6)

  title(main = main, adj = 0.07, font.main = 3, cex.main = 0.95,
        col.main = title_col, line = -0.5)

  if (plot_surtitle) {
    title(main = surtitle, adj = 0, font.main = 1, xpd = NA,
          line = 1)
  }
  if (legend) {
    legend("topright", inset = c(0.1, 0.1),
           y.intersp = 1.5,
           legend = c("Hospital", "Care home"), text.font = 3,
           fill = col[c("hosp", "comm")],
           bty = "n", cex = 1)
  }
}


plot_annotated_rt <- function(samples, Rt_general, plot_num = "",
                              col, legend = FALSE, regions) {
  x <- sircovid::sircovid_date_as_date(samples$london$trajectories$date)
  y <- sapply(Rt_general, rowMeans)
  policy_dates <- as.Date(c("2020-03-16", "2020-03-23", "2020-03-25",
                            "2020-05-11", "2020-06-15", "2020-07-04",
                            "2020-08-03", "2020-09-01", "2020-09-14",
                            "2020-10-14", "2020-11-05"))
  ## Start from blank plot with appropriate axes
  ylim <- c(0, 4.5)
  xlim <- as.Date(c("2020-03-01", "2020-12-01"))
  x[1] <- xlim[1]
  w <- x <= xlim[2]

  plot(x, rep(0, length(x)), type = "n", ylim = ylim, xlim = xlim, xaxt = "n",
       xlab = "", ylab = expression(italic({R ^ eff}[t])), yaxp = c(0, 4, 4),
       las = 1)
  axis_dates <- as.Date(
    c("2020-03-01", "2020-04-01", "2020-05-01", "2020-06-01",
      "2020-07-01", "2020-08-03", "2020-09-01", "2020-10-01",
      "2020-11-01", "2020-12-01"))
  axis(1, axis_dates, format(axis_dates, "%b"))

  ## Add policy date lines
  y_label <- c(3.9, 4, 2.2, 1.5, 2.3, 3.6, 2.3, 3.4, 4, 2.5, 2)
  segments(policy_dates, y0 = -0.1, y1 = y_label, lty = 3)

  ## Plot effective Rt over time
  matlines(x[w], y[w, ], type = "l", lty = 1, col = col, lwd = 1.5, lend = 1)

  ## Add line for Reff = 1
  segments(xlim[1], 1, xlim[2], lty = 2)

  ## Annotate with important dates
  ## See individual task script.R for details of these dates
  text(x = policy_dates,
       y = y_label,
       labels = c("Mar 16:\nMovement\nrestrictions",
                  "Mar 23:\nLockdown\nannounced",
                  "Mar 25:\nLockdown\nin full effect",
                  "May 11:\nInitial easing",
                  "Jun 15:\nShops\nre-open",
                  "Jul 04:\nPubs and\nrestaurants\nre-open",
                  "Aug 03:\nEat-out-to-\nhelp-out",
                  "Sep 01:\nSchools\nreturn",
                  "Sep 14:\nRule of 6",
                  "Oct 14:\nTiers\nintroduced",
                  "Nov 05:\nSecond\nlockdown"),
       pos = c(2, 4, 4, 2, 2, 2, 3, 2, 3, 3, 4), cex = 1, font = 3)

  ## Add legend underneath plot
  ranked_regions <- colnames(y)[order(y[nrow(y), ], decreasing = TRUE)]
  legend("topright", inset = c(0.03, 0),
         ncol = 1, legend = regions[ranked_regions, "label"],
         fill = col[ranked_regions],
         bty = "n", cex = 0.9, text.font = 3, y.intersp = 1.2)
}


plot_paper_fig2 <- function(agg_samples, agg_data, regions, region_names) {
  par(mar = c(3, 4, 3, 3), mgp = c(1.7, 0.7, 0))
  layout(matrix(c(1, 1, 1, 1,
                  2, 3, 4, 5,
                  6, 7, 8, 9), nrow = 3, byrow = TRUE),
         heights = c(2, 1, 1))

  cols <- c(paper_seq_palette(3)[2], paper_seq_palette2(3)[2])

  ## Plot A
  plot_carehome_vs_comm_incid1(agg_samples$england, max_date = "2020-12-01",
                               exch_col = cols[1],
                               ch_col = cols[2], alpha = 0.4,
                               ymax = 2e5)
  title("A) England daily COVID-19 infections", adj = 0, line = 1,
        font.main = 3)

  legend("topleft", inset = c(0.1, 0.01), fill = cols,
         legend = c("Community", "Care home residents"),
         bty = "n", text.font = 3, y.intersp = 1.5)
  regions <- regions[region_names, ]

  ## Plot B - H
  par(mar = c(2, 3, 1, 0))
  mapply(FUN = plot_pillar2,
         sample = agg_samples[region_names],
         data = agg_data[region_names],
         title = sprintf("%s) %s",
                         LETTERS[seq_along(region_names) + 1L],
                         regions$name),
         MoreArgs = list(xlim = as.Date(c("2020-06-01", "2020-12-01")),
                         ylim = c(0, 25), col = paper_seq_palette3(3)[2],
                         col_data = grey(0.2), ylab = "Pillar-2 positive (%)"))
}


plot_carehome_vs_comm_incid1 <- function(sample, main,
                                         exch_col, ch_col,
                                         max_date,
                                         ymax,
                                         alpha = 0.3,
                                         ysf_exch = 1e4, ysf_ch = 500,
                                         ymax_ch_fac = 100) {
  S_CHR <- sample$trajectories$state["S_CHR", , ]
  incid_CHR <-  rbind(0, apply(S_CHR[1, 1] - S_CHR, 1, diff))
  incid_CHR[1:2, ] <- NA
  incid_exCHR <- t(sample$trajectories$state["infections_inc", , ]) - incid_CHR

  dx <- sircovid::sircovid_date_as_date(sample$trajectories$date)

  extract_qs <- function(res) {
    res <- res[dx <= max_date, ]
    ps <- seq(0.025, 0.975, 0.005)
    apply(res,  MARGIN = 1, FUN = quantile, ps, na.rm = TRUE)
  }

  qs_CHR <- extract_qs(incid_CHR)
  qs_exCHR <- extract_qs(incid_exCHR)

  add_trajectories <- function(qs, col) {
    cols <- add_alpha(rep(col, 2), alpha)
    x <- dx[dx <= max_date]
    CI_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = cols,
             horiz = FALSE, leg = FALSE)
    lines(x, qs["50.0%", ], col = col, lty = 1, lwd = 1.5, lend = 1)
  }

  xlim <- c(as.Date("2020-03-15"), max_date)
  if (is.null(ymax)) {
    ymax <- ceiling(max(incid_exCHR, na.rm = TRUE) / ysf_exch) * ysf_exch
  }
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = c(0, ymax),
       axes = FALSE,
       xlab = "", ylab = "")
  axis.Date(1, at = seq.Date(as.Date("2020-03-01"), xlim[2], by = "month"))
  add_trajectories(qs = qs_exCHR, col = exch_col)
  yaxis <- seq(0, ymax, length.out = 3)
  axis(side = 2, at = yaxis,
       labels = formatC(yaxis, format =  "d", big.mark = ","),
       col = exch_col, las = 1, lwd = 2)
  par(new = TRUE)
  ymax <- ceiling(max(incid_CHR, na.rm = TRUE) / ysf_ch) * ysf_ch
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = c(0, ymax),
       axes = FALSE,
       xlab = "", ylab = "")
  yaxis <- seq(0, ymax, length.out = 3)
  axis(side = 4, at = yaxis, col = ch_col, las = 1, lwd = 2)
  add_trajectories(qs = qs_CHR, col = ch_col)
}


plot_pillar2 <- function(sample, data, col, col_data, title, xlim, ylim, ylab,
                         alpha = 0.3) {
  npos <- data$fitted[, "pillar2_over25_pos"]
  ntot <- data$fitted[, "pillar2_over25_tot"]

  dx <- as.Date(data$fitted$date_string)
  dy <- npos / ntot * 100

  trajectories <- sample$trajectories$state
  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)

  model_params <- sample$predict$transform(sample$pars[1, ])

  prop_noncovid_sympt <- sample$pars[, "prop_noncovid_sympt"]

  pos <- trajectories["sympt_cases_over25_inc", , ]
  neg <- (sum(model_params$N_tot[6:19]) - pos) * prop_noncovid_sympt

  res <- (pos * model_params$pillar2_sensitivity +
            neg * (1 - model_params$pillar2_specificity)) / (pos + neg) * 100


  ps <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       las = 1,
       main = "",
       font.main = 1,
       xlab = "", ylab = ylab)

  CI_bands(qs[c("2.5%", "25%", "75%", "97.5%"), ], x,
           cols = add_alpha(rep(col, 2), alpha),
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50%", ], col = col, lty = 1, lwd = 1.5, lend = 1)
  lines(dx, dy, col = col_data, lend = 1)
  title(title, adj = 0.1, line = 0, font.main = 3, cex.main = 1)
}


plot_paper_fig4 <- function(samples, ifr_t, regions, region_names,
                            region_cols) {
  nms <- rev(region_names)
  labels <- regions[nms, "label"]

  layout(mat = matrix(c(1, 2, 3,
                        4, 5, 5),
                      byrow = TRUE, nrow = 2),
         widths = c(1, 1, 1.8),
         heights = c(5, 5))

  ## Plot IFR/HR gen pop
  col_ifr <- paper_seq_palette(3)[2]
  col_ifr_end <- paper_seq_palette3(3)[2]
  col_ihr <- paper_seq_palette2(3)[2]

  ## Plot age dist
  par(mar = c(3, 2.5, 2, 0), mgp = c(1.7, 0.7, 0), bty = "n")

  ages <- seq(0, 80, 5)
  at_chr <- 87
  qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)

  extract_by_age <- function(region, what) {
    sapply(region[sprintf("%s_%s", what, ages)], function(x) x[nrow(x), ])
  }

  ihr_by_age <- sapply(region_names, function(x)
    apply(extract_by_age(ifr_t[[x]], "IHR"), 2, median))

  ihr_chr <- sapply(region_names, function(x)
    median(tail(ifr_t[[x]]$IHR_CHR)))

  ifr_by_age <- sapply(region_names, function(x)
    colMedians(extract_by_age(ifr_t[[x]], "IFR")))

  plot(0, 0, type = "n", las = 1,
       xlim = c(0, at_chr), ylim = c(0, 50), xlab = "Age", ylab = "(%)",
       xaxt = "n")
  segments(x0 = 0, x1 = at_chr, y0 = seq(0, 50, 10), col = grey(0.9))
  axis(1, seq(0, 80, 20), labels = c(seq(0, 60, 20), "80+"))
  matlines(ages, ihr_by_age, col = region_cols[region_names], lty = 1,
           lend = 1, lwd = 1.5)
  legend("topleft", fill = region_cols, bty = "n",
         xjust = 0, ncol = 2,
         legend = regions[names(region_cols), "label"], text.font = 3)

  title("A) Regional IHR by age", font.main = 1, adj = 0, line = 1)

  plot(0, 0, type = "n", las = 1,
       xlim = c(0, at_chr), ylim = c(0, 15), xlab = "Age",
       xaxt = "n", ylab = "(%)")
  segments(x0 = 0, x1 = at_chr, y0 = seq(0, 50, 5), col = grey(0.9))
  axis(1, seq(0, 80, 20), labels = c(seq(0, 60, 20), "80+"))
  matlines(ages, ifr_by_age, col = region_cols[region_names], lty = 1,
           lend = 1, lwd = 1.5)


  title("B) Regional IFR by age", font.main = 1, adj = 0, line = 1)

  ## Plot IFR/HR by age
  ifr_by_age <- sapply(ifr_t$england[sprintf("IFR_%s", ages)],
                       function(x) x[nrow(x), ])
  q_ifr_by_age <- apply(ifr_by_age, 2, quantile, qs)
  ihr_by_age <- sapply(ifr_t$england[sprintf("IHR_%s", ages)],
                       function(x) x[nrow(x), ])
  q_ihr_by_age <- apply(ihr_by_age, 2, quantile, qs)

  plot(0, 0, type = "n", las = 1,
       xlim = c(0, at_chr), ylim = c(0, 50), xlab = "Age", ylab = "(%)",
       xaxt = "n")
  segments(x0 = 0, x1 = at_chr, y0 = seq(0, 50, 10), col = grey(0.9))

  CI_bands(y = ages, quantiles = q_ifr_by_age, leg = FALSE, horiz = FALSE,
           cols = paper_seq_palette(4)[-4])
  lines(ages, q_ifr_by_age["50%", ], col = col_ifr, lend = 1)
  CI_bands(y = ages, quantiles = q_ihr_by_age, leg = FALSE, horiz = FALSE,
           cols = paper_seq_palette2(4)[-4])
  lines(ages, q_ihr_by_age["50%", ], col = col_ihr, lend = 1)

  plot_CI_bar(tail(ifr_t$england$IFR_CHR, 1), at = at_chr, horiz = FALSE,
              col = col_ifr, pt_col = col_ifr)
  plot_CI_bar(tail(ifr_t$england$IHR_CHR, 1), at = at_chr, horiz = FALSE,
              col = col_ihr, pt_col = col_ihr)
  axis(1, seq(0, 80, 10), labels = c(seq(0, 70, 10), "80+"))
  axis(1, at_chr, labels = "CHR")
  legend("topleft", inset = c(0.01, 0.02),
         legend = c("Infection hospitalisation ratio (IHR)",
                    "Infection fatality ratio (IFR)"),
         fill = c(col_ihr, col_ifr), y.intersp = 1.5, text.font = 3, bty = "n")
  title("C) England severity by age and in care home residents",
        font.main = 1, adj = 0, line = 1)

  ifr_end <- lapply(ifr_t[nms], function(x) tail(x$IFR_AR_all, 1))
  ihr <- lapply(ifr_t[nms], function(x) tail(x$IHR_AR_all, 1))
  ifr_start <- lapply(ifr_t[nms], function(x) head(x$IFR_AR_all, 1))
  par(mar = c(2.5, 3, 2, 0), mgp = c(1.7, 0.7, 0), bty = "n")

  compare_ci_bars(list(ihr), c(col_ihr), regions = regions,
                  xlab = "(%)", yaxt = TRUE, grid = 1,
                  labels = regions[nms, "label"], mgp = c(1.7, 0.7, 0))
  title(main = "D) Regional IHR", line = 0,  font.main = 1, adj = 0)

  ## Plot IFR_t
  par(mar = c(2.5, 3, 2, 1), mgp = c(1.7, 0.7, 0), bty = "n")
  plot_ifr_t(ifr_t$england, ifr_name = "IFR_AR_all", ylim = c(0, 2.5),
             forecast_until = "2020-12-01", col = col_ifr,
             include_forecast = FALSE, ylab = "(%)")

  x_start <- as.Date(sircovid::sircovid_date_as_date(ifr_t$england$date[2]))
  x_end <- as.Date("2020-12-01")

  add_point <- function(x, jitter, ifr_date) {
    y <- sapply(ifr_date, mean)
    y0 <- sapply(ifr_date, quantile, 0.025)
    y1 <- sapply(ifr_date, quantile, 0.975)

    x0 <- x + jitter
    names(x0) <- names(y0) <- names(y0) <- names(y)

    segments(x0 = x0, y0 = y0, y1 = y1)
    points(x0, y, pch = 23,
           bg = region_cols[names(ifr_date)])

    text(x0["london"], y0["london"], "LON", font = 3, pos = 1, cex = 0.8)

  }

  jitter_start <- c(4, -6, 0, 0, 7, -3, 4) + 1
  jitter_end <- c(9, 6, 1, -6, 0, -3, 3)

  add_point(x_start, jitter_start, ifr_start)
  add_point(x_end, jitter_end, ifr_end)

  title("E) England age-aggregated infection fatality ratio (IFR)",
        font.main = 1, adj = 0,
        line = 0)

  mean_ifr <- sapply(ifr_t[nms], function(x)
    sapply(x[sprintf("IFR_%s", seq(0, 80, 5))], function(x) mean(tail(x, 1))))
  mean_ihr <- sapply(ifr_t[nms], function(x)
    sapply(x[sprintf("IHR_%s", seq(0, 80, 5))], function(x) mean(tail(x, 1))))
}


plot_ifr_t <- function(sample_ifr_t, ifr_name, ylim = c(0, 2.5), title = "",
                       forecast_until, col,
                       alpha = 0.3, include_forecast = TRUE, ylab = "",
                       add = FALSE) {
  ## Rrop the first value at the moment the first time period is much
  ## more than one day
  ifr_t <- sample_ifr_t[[ifr_name]][-1L, ]
  date <- sample_ifr_t$date[-1L]

  if (!include_forecast) {
    i_date <- date <= sircovid::sircovid_date(forecast_until)
    ifr_t <- ifr_t[i_date, ]
    date <- date[i_date]
  }

  x <- sircovid::sircovid_date_as_date(date)
  rownames(ifr_t) <- as.character(x)

  ps <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  qs <- apply(ifr_t,  MARGIN = 1, FUN = quantile, ps, na.rm = TRUE)

  if (!add) {
    par(mgp = c(1.7, 0.5, 0), bty = "n")
    xlim <- c(x[1], as.Date(forecast_until))
    plot(xlim[1], 0, type = "n",
         xlim = xlim,
         ylim = ylim,
         main = toupper(title),
         font.main = 1,
         xlab = "", ylab = ylab,
         xaxt = "n")
    segments(x0 = xlim[1], x1 = xlim[2], y0 = seq(0, 4, 0.5), col = grey(0.9))
    axis.Date(1, at = seq(as.Date("2020-04-01"), as.Date(forecast_until),
                          by = "2 month"), format = "%b")
  }
  cols <- add_alpha(rep(col, 2), alpha)

  CI_bands(qs[c("2.5%", "25%", "75%", "97.5%"), ], x, cols = cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50%", ], col = col, lty = 1, lwd = 1.5, lend = 1)
}


plot_cum_inf_ch <- function(samples, data, date = NULL,
                            col_comm, col_ch, regions) {
  samples <- rev(samples)

  cum_inf_comm <- mapply(FUN = calc_cum_inf_age,
                         sample = samples, data = data,
                         MoreArgs = list(a = 17, date = date),
                         SIMPLIFY = FALSE)
  cum_inf_chr <- mapply(FUN = calc_cum_inf_age,
                        sample = samples, data = data,
                        MoreArgs = list(a = 19, date = date),
                        SIMPLIFY = FALSE)

  labels <- regions[names(samples), "label"]

  n_regions <- length(samples)
  at <- seq_len(n_regions)

  xmax <- ceiling(max(sapply(cum_inf_chr, quantile, 0.975)) * 10) / 10

  ylim <- c(0.5, n_regions + 0.5)
  plot(0, 0, type = "n", ylab = "", xlab = "",
       xlim = c(0, xmax), ylim = ylim, axes = FALSE)
  segments(x0 = seq(0, 100, 5), y0 = ylim[1], y1 = ylim[2], col = grey(0.9))

  axis(side = 2, at =  at, labels = labels, las = 1, font = 3,
       mgp = c(1.5, 0.7, 0))
  xaxt <- pretty(seq(0, xmax, length.out = 5))
  axis(side = 1, at = xaxt, labels = xaxt * 100, mgp = c(1.5, 0.5, 0))
  mtext("Proportion ever infected (%)", side = 1, line = 1.5, xpd = TRUE,
        cex = 0.7)

  mapply(plot_CI_bar, res = cum_inf_comm, at = at, col = col_comm,
         width = 0.2, cex = 1.7, lwd = 3)
  mapply(plot_CI_bar, res = cum_inf_chr, at = at, col = col_ch,
         width = 0.2, cex = 1.7, lwd = 3)
  at <- max(xaxt) * c(0.15, 0.6)
  mtext(side = 3, text = c("80+", "CHR"), line = -0.5, font = 3, cex = 0.7,
        at = at)
}


compare_ci_bars <- function(res, col,
                            yaxt = TRUE, grid = 1e2,
                            labels = NULL, xlab = "", xsf = 1,
                            width = 0.15, cex = 1.5, lwd = 1, yadj = -0.25,
                            mgp = c(1.7, 0.5, 0),
                            regions) {
  n_regions <- length(res[[1]])
  at <- seq_len(n_regions) + yadj
  xmax <- max(unlist(res, recursive = TRUE), na.rm = TRUE)
  xmax <- ceiling(xmax / xsf) * xsf
  xlim <- c(0, xmax)

  if (is.null(labels)) {
    labels <- regions[names(res[[1]]), "label"]
  }

  ylim <- c(0.5, n_regions + 0.5)
  plot(0, 0, type = "n", ylab = "", xlab = "",
       xlim = xlim, ylim = ylim, axes = FALSE)
  segments(y0 = ylim[1], y1 = max(at) + 0.1, x0 = seq(0, 100, grid),
           col = grey(0.9))
  if (yaxt) {
    axis(side = 2, at = at, labels = labels, las = 1, font = 3, mgp = mgp)
  }
  xaxt <- pretty(xlim, min.n = 2)
  axis(side = 1, at = xaxt, labels = xaxt, mgp = mgp)
  mtext(xlab, side = 1, line = 1.5, xpd = TRUE, cex = 0.7)

  mapply(function(res, col) {
    mapply(plot_CI_bar, res = res, at = at, col = col,
           width = width, cex = cex, lwd = lwd)
  }, res = res, col = col)
}


plot_paper_fig5 <- function(samples, data, date, region_names, maps, regions) {
  layout(mat = matrix(c(1, 2, 3, 4,
                        5, 6, 9, 10,
                        7, 8, 9, 10), byrow = TRUE, nrow = 3),
         widths = c(1, 1, 1), heights = c(1, 1, 1))

  ## Plot A-G
  par(mar = c(2, 2.5, 1.5, 0.5))
  title <- sprintf("%s) %s (%s)", LETTERS[seq_along(region_names)],
                   regions[region_names, "name"],
                   regions[region_names, "label"])

  mapply(FUN = plot_serology,
         sample = samples[region_names],
         data = data[region_names],
         region = region_names,
         title = title,
         date = date,
         ymax = 20, pos_col = paper_cols["purple"],
         inf_col = paper_cols["nowcast"], alpha = 0.3,
         plot_legend = c(TRUE, rep(FALSE, length(region_names) - 1)))

  ## Plot H
  par(mar = c(2, 2.5, 2, 1), mgp = c(1.5, 0.5, 0))
  plot_cum_inf_ch(samples[region_names], data[region_names],
                  date = date,
                  col_comm = paper_cols["hosp"], col_ch = paper_cols["comm"],
                  regions = regions)
  title(main = "H)", font.main = 1, adj = 0)
  par(mar = c(0, 0, 3, 0))
  plot_map(samples, regions = maps, date = date, legend = FALSE)
  title(main = "I) Total population", font.main = 1, line = 1, adj = 0.05)
  plot_map(samples, regions = maps, date = date, what = "S_CHR", legend = TRUE,
           add_label = FALSE)
  title(main = "J) Care home residents (CHR)", font.main = 1, line = 1,
        adj = 0.05)
}


plot_serology <- function(sample, data, date, region, title = "",
                           inf_col, pos_col,
                           alpha = 0.35,
                           ymax = 40,
                           plot_legend = FALSE) {
  extract_serodates <- function(region, data, cap = 150) {
    sero <- data$full$ntot_15_64 > cap
    sero[is.na(sero)] <- FALSE
    n <- length(sero)

    lag <- seq(4, n)

    start <- sero[lag] & !sero[lag - 1] & !sero[lag - 2] & !sero[lag - 3]
    end <- sero[lag - 3] & !sero[lag - 2] & !sero[lag - 1] & !sero[lag]

    dates <- data$full$date_string
    res <- data.frame(region = region,
                      start = dates[lag][start])

    end_dates <- dates[lag - 3][end]
    if (length(end_dates) < length(res$start)) {
      end_dates <- c(end_dates, date)
    }
    res$end <- end_dates
    res
  }

  sero_dates <- extract_serodates(region, data)
  tol <- 2

  sero_dates$start <- sero_dates$start - tol
  sero_dates$end <- sero_dates$end + tol

  rownames(data$fitted) <- data$fitted$date_string

  sero_data <- data$fitted[, c("ntot_15_64", "npos_15_64")]
  colnames(sero_data) <- c("ntot", "npos")
  sero_data[is.na(sero_data)] <- 0

  summ_serodata <- sapply(X = seq_len(nrow(sero_dates)), function(i) {
    w <- seq(from = as.Date(sero_dates[i, "start"]),
             to = as.Date(sero_dates[i, "end"]), 1)
    w <- as.character(as.Date(w))
    colSums(sero_data[w, ], na.rm = TRUE)
  })

  summ_serodata <- data.frame(sero_dates,
                              t(summ_serodata),
                              stringsAsFactors = FALSE)

  summ_serodata <- cbind(summ_serodata,
                         with(summ_serodata,
                              Hmisc::binconf(x = npos, n = ntot) * 100))
  summ_serodata$mid <- (as.numeric(as.Date(summ_serodata$start)) +
                          as.numeric(as.Date(summ_serodata$end))) / 2

  params <- sample$predict$transform(sample$pars[1, ])

  N0 <- sum(params$N_tot)

  trajectories <- sample$trajectories$state
  res <- (params$sero_sensitivity * trajectories["sero_pos", , ] +
            (1 - params$sero_specificity) *
            (params$N_tot_15_64 - trajectories["sero_pos", , ])) /
    params$N_tot_15_64 * 100
  res_infs <- trajectories["infections", , ] / N0 * 100

  ps <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)
  qs_infs <- apply(res_infs,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  pos_cols <- add_alpha(rep(pos_col, 2), alpha)
  inf_cols <- add_alpha(rep(inf_col, 2), alpha)

  ylim <- c(0, ymax)
  par(mgp = c(1.7, 0.5, 0), bty = "n")
  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)
  x <- x[x <= date]
  xlim <- c(min(x[-1L]), max(x[-1L]))
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim, las = 1,
       main = "",
       font.main = 1,
       xlab = "", ylab = "Cumulative proportion (%)")
  title(main = title, adj = 0, font.main = 1, line = 0.5,
        cex.main = 1)

  lapply(X = seq_len(nrow(summ_serodata)), FUN = function(i) {

    xx <- rep(unlist(summ_serodata[i, c("start", "end")]), each = 2)
    yy <- c(ylim, rev(ylim))

    polygon(x = xx, y = yy, col = grey(0.9), border = NA)
  })

  p_labels <- c("2.5%", "25%", "75%", "97.5%")

  CI_bands(qs[p_labels, seq_along(x)], x, cols = pos_cols,
           horiz = FALSE, leg = FALSE)
  CI_bands(qs_infs[p_labels, seq_along(x)], x, cols = inf_cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs_infs["50%", seq_along(x)], col = inf_col, lty = 1, lwd = 1.5,
        lend = 1)
  lines(x, qs["50%", seq_along(x)], col = pos_col, lty = 1, lwd = 1.5,
        lend = 1)
  with(summ_serodata, {
    segments(x0 = mid, y0 = Lower, y1 = Upper, lty = 1, lend = 1)
    points(x = mid, y = PointEst, pch = 18)
    points(x = rep(mid, 2), y = c(Lower, Upper), pch = "-")
  })

  if (plot_legend) {
    leg_cols <- c(inf_col, pos_col)
    legend("top", legend = c("Infected", "Seropositive"),
           cex = 1, x.intersp = 2, ncol = 2,
           fill = add_alpha(leg_cols, alpha * 2),
           border = leg_cols,
           box.col = "white")
  }
}


plot_map <- function(samples, regions, date, legend = TRUE, what = "total",
                     max_incid = NULL, add_label = TRUE) {

  extract_cumincid <- function(sample) {
    index_S <- grep("^S_", names(sample$predict$index))
    S <- sample$trajectories$state[index_S, , ]
    n_groups <- sample$predict$transform(sample$pars[1, ])$n_groups
    if (n_groups != length(index_S)) {
      ## sum across vacc classes
      S <- array(S, c(n_groups, dim(S)[1] / n_groups, dim(S)[-1]))
      S <- apply(S, c(1, 3, 4), sum)
      rownames(S) <- grep("^S_", names(sample$predict$index),
                          value = TRUE)[seq.int(n_groups)]
    }
    N0 <- S[, , 1]
    sircovid_date <- sircovid::sircovid_date(date)
    I <- N0 - S[, , sample$trajectories$date == sircovid_date]
    cumincid <- c(rowMeans(I / N0), total = mean(colSums(I) / sum(N0[, 1])))

    cumincid
  }

  cumincid <- sapply(samples, extract_cumincid)

  max_incid <- ceiling(max(cumincid[c("S_CHR", "total"), ]) * 100 / 5) * 5

  r <- c("london", "east_of_england", "midlands", "north_east_and_yorkshire",
         "north_west", "south_east", "south_west")

  dp <- 10
  ncol <- max_incid * dp
  palette <- paper_seq_palette(max_incid * dp)
  cols <- palette[round(cumincid[what, r] * 100 * dp)]
  names(cols) <- r

  w <- seq(2.5, max_incid, by = 2.5)

  w_lab <- paste0(w, "%")
  w_lab[seq(1, length(w), by = 2)] <- ""
  w_lab

  plot_lims <- regions@bbox

  raster::plot(regions, col = cols, border = "white")

  calc_pos <- function(x, w) {
    x["min"] * (1 - w) + x["max"] * w
  }

  if (add_label) {
    x_lab <- calc_pos(plot_lims["x", ],
                      c(0.68, 0.51, 0.99, 0.85, 0.6, 0.63, 0.37))
    y_lab <- calc_pos(plot_lims["y", ],
                      c(0.24, 0.22, 0.288, 0.4, 0.45, 0.7, 0.66))

    text(c("SE", "SW", "LON", "EE", "MID", "NE", "NW"),
         x = x_lab, y = y_lab, font = 3, col = "black")
  }

  if (legend) {
    legend("topleft", title = "Cumulative\nincidence\n", title.adj = 0,
           legend = rev(w_lab), pch = 15, adj = 0.5, x.intersp = 2,
           cex = 1, pt.cex = 2, inset = c(0, 0.1),
           col = rev(palette[w * dp]), y.intersp = 0.8, bty = "n", xpd = NA)
  }
}


plot_supplement_testing <- function(sample, data, title, col, col_data,
                                    alpha = 0.3) {
  ylim <- c(0, 5)
  xlim <- as.Date(c("2020-06-01", "2020-12-01"))

  ## Row 1: Pillar 2 testing
  plot_pillar2(sample, data, col, col_data, title = title,
               xlim, ylim = c(0, 25), ylab = "Pillar-2 positive (%)")

  ## Row 2: REACT-1 prevalence
  days_to_agg <- 3
  min_agg_tot <- 200
  npos <- data$fitted[, "react_pos"]
  ntot <- data$fitted[, "react_tot"]

  npos[is.na(npos)] <- 0
  ntot[is.na(ntot)] <- 0

  dx <- as.Date(data$fitted$date_string)

  ## aggregate
  agg_dates <- data.frame(start = dx[seq(1, length(dx), days_to_agg)])
  agg_dates$end <- agg_dates$start + days_to_agg - 1
  agg_dates$end[length(agg_dates$end)] <- dx[length(dx)]
  agg_dates$mid <- floor(rowMeans(apply(agg_dates, 2, sircovid::sircovid_date)))
  agg_dates$mid <- sircovid::sircovid_date_as_date(agg_dates$mid)

  agg_dates$end <- as.Date(agg_dates$end)

  aggregate_react <- function(i) {
    agg_pos <- sum(npos[dx >= agg_dates$start[i] & dx <= agg_dates$end[i]])
    agg_tot <- sum(ntot[dx >= agg_dates$start[i] & dx <= agg_dates$end[i]])
    if (agg_tot > min_agg_tot) {
      agg_out <- c(agg_pos, agg_tot)
    } else {
      agg_out <- c(0, 0)
    }
    agg_out
  }

  agg_dates$npos <- agg_dates$ntot <- rep(0, length(agg_dates$mid))

  agg_dates[, c("npos", "ntot")] <- t(sapply(seq_len(length(agg_dates$mid)),
                                             aggregate_react))

  cis <- Hmisc::binconf(x = agg_dates$npos, n = agg_dates$ntot) * 100
  dy <- cis[, "PointEst"]
  lower <- cis[, "Lower"]
  upper <- cis[, "Upper"]
  dy[agg_dates$ntot == 0] <- NA
  dx <- agg_dates$mid

  trajectories <- sample$trajectories$state

  model_params <- sample$predict$transform(sample$pars[1, ])

  pos <- trajectories["react_pos", , ]
  neg <- (sum(model_params$N_tot[2:18]) - pos)

  res <- (pos * model_params$react_sensitivity +
          neg * (1 - model_params$react_specificity)) / (pos + neg) * 100

  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)

  ps <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  ## Remove 1st one currently as this is time period is much more
  ## than one day
  plot(xlim[1], 0, type = "n",
       xlim = as.Date(c("2020-04-25", date)),
       ylim = ylim, las = 1,
       xlab = "", ylab = "REACT-1 positive (%)")

  CI_bands(qs[c("2.5%", "25%", "75%", "97.5%"), ], x,
           cols = add_alpha(rep(col, 2), alpha),
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50%", ], col = col, lty = 1, lwd = 1.5, lend = 1)
  segments(x0 = dx, y0 = lower, y1 = upper, col = col_data)
  points(dx, dy, pch = 23, bg = col_data, cex = 0.8, lwd = 0.6)
}
