
#### Import the pipe operator from magrittr ####
#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#### Import functions from packages ####
#' @importFrom dplyr bind_rows case_when desc mutate n n_distinct rename
#' @importFrom stats quantile rmultinom runif sd setNames
#' @importFrom tibble tibble
#' @importFrom tidyr separate
#' @importFrom utils read.table write.table




#### Declare names of columns  to keep CRAN checks from barking  ####
# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if(getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "AbsoluteIndex",
      "AlleIdx",
      "AlleLine",
      "Allele",
      "Chrom",
      "Freq",
      "Kid",
      "LocIdx",
      "LocLine",
      "LocScrub",
      "Locus",
      "Ma",
      "Manew",
      "Pa",
      "Panew",
      "Pos",
      "Sex",
      "Sexnew",
      "TwinStatus",
      "X1",
      "X2",
      "X3",
      "X4",
      "X5",
      "X6",
      "alleidx",
      "chrom",
      "ender",
      "list_name",
      "loc_name",
      "locidx",
      "newfreq",
      "num_bases",
      "num_markers",
      "pedname",
      "pos",
      "scaled_length"
    )
  )
}

