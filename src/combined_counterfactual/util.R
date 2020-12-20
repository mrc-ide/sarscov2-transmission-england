switch_levels <- function(x) {
  nms <- names(x[[1]]) %||% seq_along(x[[1]])
  y <- lapply(nms, function(z) lapply(x, "[[", z))
  names(y) <- nms
  y
}


`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}


aggregate_data <- function(data, region_names, type = "full") {
  d <- lapply(data[region_names], "[[", type)
  x <- d[[1]]

  date_names <- grep("date", names(x), value = TRUE)
  data_names <- setdiff(names(x), date_names)

  dim_t <- min(vapply(d, nrow, integer(1)))

  d <- lapply(d, function(x) x[seq_len(dim_t), data_names])

  ret <- Reduce("+", d)
  ret$deaths <- pmax(ret$deaths, ret$deaths_hosp + ret$deaths_comm,
                     na.rm = TRUE)
  data.frame(ret, x[, date_names])
}
