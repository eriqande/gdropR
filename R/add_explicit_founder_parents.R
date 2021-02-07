#' explicitly add parents of 0 to the founders of a pedigree
#'
#' This seems to be required by Mendel.
#' @param P pedigree with columns Kid Pa Ma
#' @export
add_explicit_founder_parents <- function(P) {
  founders <- setdiff(c(P$Ma, P$Pa), P$Kid)

  ret <- bind_rows(
    tibble(
      Kid = founders,
      Pa = "0",
      Ma = "0"
    ),
    P
  )
  ret
}
