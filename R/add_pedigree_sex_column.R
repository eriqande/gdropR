
#' given a pedigree as Kid Pa Ma, this adds a column called Sex
#'
#' This just uses the occurrence of each each individual in the
#' pedigree as Ma or Pa to figure out its sex.  Individuals with
#' no offspring are randomly assigned to be females (2) with probability
#' fem_prob, and are males (1) otherwise.
#' This function does not check for sex inconsistencies.
#' @param P the pedigree with columns Kid Pa Ma
#' @param fem_prob prob that unknown sex indiv is assigned as female.
#' @export
add_pedigree_sex_column <- function(
  P,
  fem_prob = 0.5
) {

  P2 <- P %>%
    mutate(
      Sex = case_when(
        Kid %in% P$Pa ~ 1L,
        Kid %in% P$Ma ~ 2L,
        TRUE ~ NA_integer_
      )
    )

  na <- is.na(P2$Sex)
  nn <- sum(na)

  P2$Sex[na] <- sample(
    x = 1:2,
    size = nn,
    replace = TRUE,
    prob = c(1 - fem_prob, fem_prob)
  )

  P2
}
