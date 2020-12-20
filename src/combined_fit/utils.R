colMedians <- function(x) {
  apply(x, 2, median)
}


rowMedians <- function(x) {
  apply(x, 1, median)
}


vnapply <- function(x, fun, ...) {
  vapply(x, fun, numeric(1), ...)
}
