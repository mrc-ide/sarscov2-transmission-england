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

  par_labels <-
    expression(alpha_D = alpha[D], alpha_H = alpha[H],
               beta1 = beta[1], beta2 = beta[2], beta3 = beta[3],
               beta4 = beta[4], beta5 = beta[5], beta6 = beta[6],
               beta7 = beta[7], beta8 = beta[8], beta9 = beta[9],
               beta10 = beta[10], beta11 = beta[11], beta12 = beta[12],
               eps = epsilon, m_CHR = m[CHR], m_CHW = m[CHW],
               mu_D = mu[D], mu_ICU = mu[ICU],
               p_G_D = p[G[D]], p_G_D_CHR = p[G[D]] ^ {CHR},
               p_H = p[H], p_H_CHR = p[H] ^ {CHR}, p_H_D = p[H[D]],
               p_ICU = p[ICU], p_ICU_D = p[ICU[D]],
               p_NC = p[NC], p_W_D = p[W[D]],
               rho_pillar2_tests = rho[P2[test]],
               start_date = t[0])

  stopifnot(length(setdiff(par_names, names(par_labels))) == 0)

  ch_names <- c("eps", "m_CHW", "m_CHR")
  mu_names <- par_names[substr(par_names, 1, 3) == "mu_"]
  alpha_names <- par_names[substr(par_names, 1, 6) == "alpha_"]
  beta_names <- par_names[substr(par_names, 1, 4) == "beta"]
  beta_names <- beta_names[order(as.numeric(gsub("beta", "", beta_names)))]
  p_names <- par_names[substr(par_names, 1, 2) == "p_"]

  p_max <- rep(1, length(p_names))
  p_max[p_names == "p_NC"] <- 0.01
  p_max[p_names == "p_G_D"] <- 0.05

  n_regions <- length(samples)

  extract_sample <- function(par_name, factor = 1) {
    lapply(X = samples, FUN = function(x)
      as.numeric(x$pars[, par_name]) / factor)
  }

  numeric_start_date <- extract_sample(par_name = "start_date")

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

  plot_par <- function(par_name, xmin = 0, xmax = 0.1, nxaxis = 4) {
    par <- extract_sample(par_name)

    plot(0, 0, type = "n",
         ylab = "",
         xlab = par_labels[par_name],
         xaxt = "n",
         xlim = c(xmin, xmax),
         ylim = ylim,
         yaxt = "n")
    xaxis_labels <- pretty(c(xmin, xmax), n = nxaxis, min.n = 2)
    axis(1, xaxis_labels, prettyNum(xaxis_labels))

    jitter <- 0.5
    regions <- names(par)
    hp <- subset(hps, name == par_name)
    rownames(hp) <- hp$region
    hp <- hp[regions, ] # sort in correct order
    if (hp$type[1] == "beta") {
      shape1 <- hp$beta_shape1
      shape2 <- hp$beta_shape2
      if (!(all(shape1 == 1) && all(shape2 == 1))) {
        prior <- mapply(qbeta,
                        shape1 = shape1,
                        shape2 = shape2,
                        MoreArgs = list(p = c(0.025, 0.975)),
                        SIMPLIFY = TRUE)


        segments(x0 = prior[1, ],
                 y0 = at- jitter, y1 = at + jitter,
                 col = col_line, lty = 2, lwd = 1, lend = 2)
        segments(x0 = prior[2, ],
                 y0 = at- jitter, y1 = at + jitter,
                 col = col_line, lty = 2, lwd = 1, lend = 2)
      }
    }
    if (hp$type[1] == "gamma") {
      shape <- hp$gamma_shape
      scale <- hp$gamma_scale
      prior <- mapply(qgamma,
                      shape = shape,
                      scale = scale,
                      MoreArgs = list(p = c(0.025, 0.975)),
                      SIMPLIFY = TRUE)


      segments(x0 = prior[1, ],
               y0 = at- jitter, y1 = at + jitter,
               col = col_line, lty = 2, lwd = 1, lend = 2)
      segments(x0 = prior[2, ],
               y0 = at- jitter, y1 = at + jitter,
               col = col_line, lty = 2, lwd = 1, lend = 2)
    }

    mapply(FUN = plot_CI_bar, res = par, at = at, width = 0.1)

  }

  hps <- read.csv("parameters_prior.csv", stringsAsFactors = FALSE)


  npar <- length(par_names)
  nrow <- 6
  plot_per_row <- ceiling(npar / nrow)
  nplots <- (plot_per_row + 1) * nrow

  layout(mat = matrix(seq_len(nplots), nrow = nrow, byrow = TRUE),
         widths = c(1, rep(2, plot_per_row)))

  at <- seq_len(n_regions)
  plot_axis()

  ## start_date plot

  plot(x = sircovid::sircovid_date("2020-01-01"),
       y = 1,
       type = "n",
       ylab = "",
       xlab = par_labels["start_date"],
       xlim = sircovid::sircovid_date(c("2020-01-01", "2020-03-01")),
       ylim = ylim,
       axes = FALSE)
  mapply(plot_CI_bar, res = numeric_start_date, at = at, width = 0.1)
  axis_dates <- c("2020-01-01", "2020-02-01", "2020-03-01")
  axis(1, at = sircovid::sircovid_date(axis_dates),
        labels = format.Date(axis_dates, "%b"))

  beta_xmax <- 0.15
  mapply(plot_par, par_name = beta_names[1:4], xmax = beta_xmax)

  ## 2nd row
  plot_axis()
  mapply(plot_par, par_name = beta_names[5:9], xmax = beta_xmax)

  ## 3rd row
  plot_axis()
  mapply(plot_par, par_name = beta_names[10:12], xmax = beta_xmax)
  mapply(plot_par, par_name = alpha_names[1:2], xmax = c(1, 0.1))

  ## 4th row
  plot_axis()
  plot_par(par_name = "rho_pillar2_tests", xmax = 0.01)
  idx_p <- 1:4
  mapply(plot_par, par_name = p_names[idx_p], xmax = p_max[idx_p])

  ## 5th row
  plot_axis()
  idx_p <- 5:9
  mapply(plot_par, par_name = p_names[idx_p], xmax = p_max[idx_p])

  ## 5th row
  plot_axis()
  mapply(plot_par, par_name = mu_names[1:2], xmax = 1)
  mapply(plot_par, par_name = ch_names[1:3], xmax = c(1, 2.5e-5, 2.5e-5),
         nxaxis = c(5, 2, 2))

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
