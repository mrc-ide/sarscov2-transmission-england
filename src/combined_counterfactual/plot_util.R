## Colour palettes
paper_cols <- c(red = "#E75B64FF", orange = "#E48C2AFF", yellow = "#DCCA2CFF",
                green = "#58A449FF", blue1 = "#67B9E9FF", blue2 = "#AE93BEFF",
                purple = "#5C5992FF", carehomes = "#CEC917FF",
                hosp = "#4D6D93FF", nowcast = "#278B9AFF",
                data_fitted = "#E48C2AFF", data_unfitted = "#228B22FF")
paper_seq_palette <- colorRampPalette(c("#C0DDE1FF", "#14454CFF"))
paper_seq_palette2 <- colorRampPalette(c("#D7CADEFF", "#1F1935FF"))
paper_seq_palette3 <- colorRampPalette(c("#213148FF", "#A1B1C8FF"))

## Make a colour semi-transparent
add_alpha <- function(col, alpha = 1) {
  apply(sapply(col, col2rgb) / 255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha = alpha))
}

## Create shaded polygons denoting density
CI_bands <- function(quantiles, y, palette = NULL, cols = NULL, leg = TRUE,
                     leg_y = 0, leg_x = 1, horiz = TRUE, ...) {
  yy <- c(y, rev(y))
  yy <- c(yy, yy[1])
  n_bands <- (nrow(quantiles) - 1) / 2 + 1
  if (!is.null(palette)) {
    cols <- do.call(what = palette,
                    args = list(n = n_bands))
  }

  for (band in seq_len(n_bands)) {
    x1 <- quantiles[band, ]
    x2 <- quantiles[nrow(quantiles) + 1 - band, ]

    x2 <- rev(x2)
    x2 <- c(x2, x1[1])
    if (horiz) {
      polygon(y = yy,
              x = c(x1, x2),
              col = cols[band],
              border = NA)
    } else {
      polygon(x = yy,
              y = c(x1, x2),
              col = cols[band],
              border = NA)
    }

  }
  if (leg) {
    leg_cols <- which(row.names(quantiles) %in% leg)
    leg <- c(row.names(quantiles)[1], seq(5, 50, 5), "%")
    leg[seq(2, 10, 2)] <- ""
    legend(y = leg_y,
           x = leg_x,
           pch = 15,
           col = cols[leg_cols],
           legend = leg,
           border = NA,
           btx = "n", ...)
  }
}


## Add shaded bars to forest plot
plot_CI_bar <- function(res,  at, width = 1,
                        min = 0.025, max = 0.975, col = "grey20",
                        segments = FALSE, pt_col = NULL, horiz = TRUE, ...) {
  cols  <- c("grey80", col)
  qs <- quantile(res,
                 probs = seq(min, max, by = 0.005),
                 na.rm = TRUE)

  palette <- colorRampPalette(cols)
  if (segments) {
    segments(y0 = at, x0 = min(res), x1 = max(res), col = cols[2])
    points(y = rep(at, 2), x = range(res), col = cols[2], pch = "-")
  }
  CI_bands(quantiles = cbind(qs, qs),
           y = at + c(-1, 1) * width,
           palette = palette, leg = FALSE, horiz = horiz)

  if (is.null(pt_col)) pt_col <- col
  if (horiz) {
    points(y = at, x = mean(res), col = "white", pch = 23, bg = pt_col, ...)
  } else {
    points(x = at, y = mean(res), col = "white", pch = 23, bg = pt_col, ...)
  }

}


gradient_fill <- function(x, y, x0 = NULL, y0 = NULL, palette = gray.colors,
                          border = NA) {
  n <- length(y)
  if (is.null(x0)) {
    x0 <- rev(x)
  }
  if (is.null(y0)) {
    y0 <- rep(0, n)
  }
  cols <- palette(n)
  xx0 <- c(x, x0)
  yy0 <- c(y, y0)

  for (i in seq_along(x)) {
    j <- n - i + 1
    xx <- c(head(x, j), tail(x0, j))
    yy <- c(head(y, j), tail(y0, j))
    polygon(x = xx, y = yy, border = NA, col = cols[i])
  }
  polygon(x = xx0, y = yy0, border = border, col = rgb(0, 0, 0, 0))
}
