create_forest_plot <- function(samples, date) {
  ## Order by start_date
  mean_start_date <- vapply(samples, function(x)
    mean(x$pars[, "start_date"]), numeric(1))


  ordered_regions <- c(names(sort(mean_start_date[region_names],
                                    decreasing = TRUE)))
  samples <- samples[ordered_regions]

  ## Formats for region labels
  labels <- regions[ordered_regions, "label"]
  par_names <- colnames(samples[[1]]$pars)
  gamma_names <- par_names[grep(pattern = "gamma", x = par_names)]
  ch_names <- c("eps", "C_1", "C_2")
  mu_names <- c("mu_death_hosp", "mu_ICU_hosp")
  beta_names <- par_names[substr(par_names, 1, 4) == "beta"]
  p_names <- par_names[substr(par_names, 1, 2) == "p_"]

  n_regions <- length(samples)

  extract_sample <- function(par_name, factor = 1) {
    lapply(X = samples, FUN = function(x)
      as.numeric(x$pars[, par_name]) / factor)
  }

  numeric_start_date <- extract_sample(par_name = "start_date")

  if (length(gamma_names) > 0) {
    gamma_durations <- lapply(X = gamma_names, FUN = extract_duration)
    names(gamma_durations) <- gamma_names
  }

  ch <- lapply(ch_names, FUN = extract_sample)
  names(ch) <- ch_names

  ps <- lapply(p_names, FUN = extract_sample)
  names(ps) <- p_names

  betas <- lapply(beta_names, FUN = extract_sample)
  names(betas) <- beta_names

  mus <- lapply(mu_names, FUN = extract_sample)
  names(mus) <- mu_names

  prop_noncovid_sympt <- extract_sample("prop_noncovid_sympt")
  rho_pillar2_tests <- extract_sample("rho_pillar2_tests")

  ylim <- c(0.5, n_regions + 0.5)

  plot_axis <- function() {
    yax_pos <- 0.9
    plot(0, 0, type = "n",
         ylab = "",
         xlab = "",
         xlim = c(0, 0),
         ylim = ylim,
         axes = FALSE)
    axis(side = 2, at =  at, labels = labels, las = 1, pos = yax_pos)
  }

  col_line <- "red3"
  plot_p <- function(p_name, xlim  = c(0, 1)) {
    plot(0, 0, type = "n",
         ylab = "",
         xlab = p_name,
         xlim = xlim,
         ylim = ylim,
         yaxt = "n")
    jitter <- 0.5
    regions <- names(ps[[p_name]])
    hp <- subset(hps, name == p_name)
    if (nrow(hp) > 0) {
      if (nrow(hp) > 1) {
        rownames(hp) <- hp$region
        hp <- hp[regions, ] # sort in correct order
        shape1 <- hp$beta_shape1
        shape2 <- hp$beta_shape2
      } else {
        shape1 <- rep(hp$beta_shape1, length(ps))
        shape2 <- rep(hp$beta_shape2, length(ps))
      }

      prior <- mapply(qbeta,
                      shape1 = shape1,
                      shape2 = shape2,
                      MoreArgs = list(p = c(0.025, 0.975)),
                      SIMPLIFY = TRUE)

      segments(x0 = prior[1, ],
               y0 = at - jitter, y1 = at + jitter,
               col = col_line, lty = 2, lwd = 1, lend = 2)
      segments(x0 = prior[2, ],
               y0 = at - jitter, y1 = at + jitter,
               col = col_line, lty = 2, lwd = 1, lend = 2)
    }
    mapply(plot_CI_bar, res = ps[[p_name]], at = at, width = 0.1)
  }

  plot_ch <- function(ch_name, xmax = 0.1, prior_scale, prior_shape) {
    plot(0, 0, type = "n",
         ylab = "",
         xlab = ch_name,
         xlim = c(0, xmax),
         ylim = ylim,
         yaxt = "n")
    if (ch_name %in% c("C_1", "C_2")) {
      jitter <- 0.5
      regions <- names(ch[[ch_name]])
      hp <- subset(hps, name == ch_name)
      if (nrow(hp) > 0) {
        if (nrow(hp) > 1) {
          rownames(hp) <- hp$region
          hp <- hp[regions, ] # sort in correct order
          scale <- hp$gamma_scale
          shape <- hp$gamma_shape
        } else {
          scale <- rep(hp$gamma_scale, length(ps))
          shape <- rep(hp$gamma_shape, length(ps))
        }

        prior <- mapply(qgamma,
                        shape = shape,
                        scale = scale,
                        MoreArgs = list(p = c(0.025, 0.975)),
                        SIMPLIFY = TRUE)


        segments(x0 = prior[1, ],
                 y0 = at - jitter, y1 = at + jitter,
                 col = col_line, lty = 2, lwd = 1, lend = 2)
        segments(x0 = prior[2, ],
                 y0 = at - jitter, y1 = at + jitter,
                 col = col_line, lty = 2, lwd = 1, lend = 2)
      }
    }
    mapply(plot_CI_bar, res = ch[[ch_name]], at = at, width = 0.1)
  }

  plot_prop <- function(prop_name, xmax = 0.1) {
    plot(0, 0, type = "n",
         ylab = "",
         xlab = prop_name,
         xlim = c(0, xmax),
         ylim = ylim,
         yaxt = "n")
    mapply(plot_CI_bar, res = prop_noncovid_sympt, at = at, width = 0.1)
  }

  plot_mu <- function(mu_name, xmax = 0.1) {
    plot(0, 0, type = "n",
         ylab = "",
         xlab = mu_name,
         xlim = c(0, xmax),
         ylim = ylim,
         yaxt = "n")
    mapply(plot_CI_bar, res = mus[[mu_name]], at = at, width = 0.1)
  }

  plot_rho <- function(prop_name, xmax = 0.1) {
    plot(0, 0, type = "n",
         ylab = "",
         xlab = prop_name,
         xlim = c(0, xmax),
         ylim = ylim,
         yaxt = "n")
    mapply(plot_CI_bar, res = rho_pillar2_tests, at = at, width = 0.1)
  }

  plot_gamma <- function(gamma_name, xmax = 25, prior_scale, prior_shape) {
    xlab <- paste(substr(x = gamma_name, start = 7, stop = nchar(gamma_name)),
                  "(days)")
    plot(0, 0, type = "n",
         ylab = "",
         xlab = xlab,
         xlim = c(0, xmax),
         ylim = ylim,
         yaxt = "n")
    add_gamma_prior(scale = prior_scale, shape = prior_shape, ylim = ylim)

    mapply(plot_CI_bar, res = gamma_durations[[gamma_name]], at = at,
           width = 0.1)
  }

  plot_betas <- function(beta_name, xlim = c(0, 0.12)) {
    plot(0, 0, type = "n",
         ylab = "",
         xlab = beta_name,
         xlim = xlim,
         ylim = ylim,
         yaxt = "n"
    )
    jitter <- 0.5
    regions <- names(betas[[beta_name]])
    hp <- subset(hps, name == beta_name)
    if (nrow(hp) > 0) {
      if (nrow(hp) > 1) {
        rownames(hp) <- hp$region
        hp <- hp[regions, ] # sort in correct order
        scale <- hp$gamma_scale
        shape <- hp$gamma_shape
      } else {
        scale <- rep(hp$gamma_scale, length(ps))
        shape <- rep(hp$gamma_shape, length(ps))
      }

      prior <- mapply(qgamma,
                      shape = shape,
                      scale = scale,
                      MoreArgs = list(p = c(0.025, 0.975)),
                      SIMPLIFY = TRUE)

      segments(x0 = prior[1, ],
               y0 = at - jitter, y1 = at + jitter,
               col = col_line, lty = 2, lwd = 1, lend = 2)
      segments(x0 = prior[2, ],
               y0 = at - jitter, y1 = at + jitter,
               col = col_line, lty = 2, lwd = 1, lend = 2)
    }
    mapply(plot_CI_bar, res = betas[[beta_name]], at = at, width = 0.1)
  }

  hps <- read.csv("parameters_prior.csv", stringsAsFactors = FALSE)

  par(bty = "n", mar = c(3, 0, 1, 0), mgp = c(2, 0.75, 0), oma = c(2, 0, 3, 3))

  npar <- length(par_names)
  nrow <- 4
  plot_per_row <- ceil(npar / nrow)
  colwidth <- 64 / plot_per_row

  reps <- rep(c(3.5, rep(colwidth, plot_per_row)), nrow)

  layout(mat = matrix(rep(seq_along(reps), reps),
                      nrow = nrow, byrow = TRUE),
         heights = rep(3, nrow), widths = c(4, rep(1, plot_per_row * colwidth))
  )
  at <- seq_len(n_regions)
  plot_axis()

  ## start_date plot
  par(mar = c(3, 0, 1, 0.5))
  plot(x = sircovid::sircovid_date("2020-01-01"),
       y = 1,
       type = "n",
       ylab = "",
       xlab = "start_date",
       xlim = sircovid::sircovid_date(c("2020-01-15", "2020-03-01")),
       ylim = ylim,
       yaxt = "n")
  mapply(plot_CI_bar, res = numeric_start_date, at = at, width = 0.1)

  mapply(plot_betas, beta_name = beta_names[c(1, 5:9)])

  ## 2nd row
  plot_axis()
  par(mar = c(3, 0, 1, 0.5))

  mapply(plot_betas, beta_name = beta_names[c(10:12, 2:4)])

  ## gamma plot
  hps_global <- subset(hps, region == "england")
  rownames(hps_global) <- hps_global$name
  if (length(gamma_names) > 0) {
    mapply(plot_gamma, gamma_name = gamma_names, xmax = c(8, 15, 6, 6, 6, 6),
           prior_shape = hps_global[gamma_names, "gamma_shape"],
           prior_scale = hps_global[gamma_names, "gamma_scale"])
  }

  plot_prop("prop_noncovid_sympt", xmax = 5e-3)

  ## 3rd row
  plot_axis()

  plot_rho("rho_pillar2_tests", xmax = 0.01)
  ## p plot
  mapply(plot_p, p_name = p_names[1:6])

  ## 4th row
  plot_axis()
  plot_mu("mu_death_hosp", xmax = 1)
  plot_mu("mu_ICU_hosp", xmax = 1)
  mapply(plot_ch, ch_name = ch_names, xmax = c(1, 2e-5, 2e-5),
         prior_shape = hps_global[ch_names, "gamma_shape"],
         prior_scale = hps_global[ch_names, "gamma_scale"])

  mtext(text = paste("Inferred epidemic parameters for NHS regions at", date),
        side = 3, line = 0, outer = TRUE, cex = 1.1)
}


## Add prior polygons to betas
add_gamma_prior <- function(shape, scale, p = c(0.025, 0.975), ylim,
                            col_line = grey(0.2)) {
  col_shade <- "white"

  prior <- qgamma(p = p, scale = scale, shape = shape)
  med <- qgamma(p = 0.5, scale = scale, shape = shape)
  polygon(x = c(prior, rev(prior)),
          y = rep(ylim, each = 2),
          col = col_shade, border = col_shade)
  segments(x0 = c(med, prior), y0 = ylim[1], y1 = ylim[2],
           col = col_line, lty = c(1, 2, 2), lwd = 1, lend = 2)
}
