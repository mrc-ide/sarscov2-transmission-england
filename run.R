main <- function() {
  regions <- c("north_west",
               "north_east_and_yorkshire",
               "midlands",
               "east_of_england",
               "london",
               "south_west",
               "south_east")
  for (r in regions) {
    orderly::orderly_run("regional_fit",
                         parameters = list(region = r, short_run = TRUE))
  }
  for (r in regions) {
    orderly::orderly_run("regional_counterfactual",
                         parameters = list(region = r, short_run = TRUE),
                         use_draft = TRUE)
  }

  orderly::orderly_run("combined_fit",
                       parameters = list(short_run = TRUE),
                       use_draft = TRUE)
  orderly::orderly_run("combined_counterfactual",
                       parameters = list(short_run = TRUE),
                       use_draft = TRUE)
}

main()
