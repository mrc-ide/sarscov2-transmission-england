plot_traces <- function(samples, samples_thin) {
  if (is.null(samples$chain)) {
    n_chains <- 1L
  } else {
    ## We assume this below
    n_chains <- length(unique(samples$chain))
    stopifnot(
      identical(samples$chain,
                rep(seq_len(n_chains),
                    each = length(samples$chain) / n_chains)))
  }
  cols <- rev(viridisLite::viridis(n_chains))

  eff_size <- function(samples_thin) {
    chains <- lapply(unname(split(data.frame(samples$pars), samples$chain)),
                     coda::as.mcmc)
  }

  plot_traces1 <- function(name) {
    if (name == "log_likelihood") {
      traces <- matrix(samples$probabilities[, "log_likelihood"],
                       ncol = n_chains)
    } else {
      traces <- matrix(samples$pars[, name], ncol = n_chains)
      thinned <- matrix(samples_thin$pars[, name], ncol = n_chains)
      ess <- round(sum(coda::effectiveSize(coda::as.mcmc(thinned))))

    }


    main <- ifelse(name == "log_likelihood", "",
                   paste("ess =", ess))
    matplot(traces, type = "l", lty = 1,
              xlab = "Iteration", bty = "n",
              ylab = name, col = cols,
              main = main,
              font.main = 1)
  }

  par(mfrow = c(5, 7),
      mar = c(3, 3, 2, 1),
      mgp = c(2, 0.5, 0),
      oma = c(1, 1, 1, 1))
  for (nm in colnames(samples$pars)) {
    plot_traces1(nm)
  }
  plot_traces1("log_likelihood")
  plot.new()
  legend("left", fill = cols, bty = "n",
         legend = paste("chain", seq_len(n_chains)))
  rhat <- gelman_diagnostic(samples)
  if (!is.null(rhat)) {
    mtext(side = 3, text = paste("Rhat =", round(rhat$mpsrf, 1)))
  }
}

gelman_diagnostic <- function(samples) {
  if (is.null(samples$chain)) {
    return(NULL)
  }
  chains <- lapply(unname(split(data.frame(samples$pars), samples$chain)),
                   coda::as.mcmc)
  tryCatch(
    coda::gelman.diag(chains),
    error = function(e) NULL)
}
