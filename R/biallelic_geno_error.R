#' Simulate genotyping error on biallelic markers
#'
#' Simple function to simulate genotyping errors.
#'
#' @param genos a vector of genotypes. The genotypes _must_ be given as `1/1`, `1/2`, `2/1` or `2/2`.
#' @param het_miscall the rate at which heterozygotes are called as
#' homozygotes.
#' @param hom_miscall the rate at which a homozygote is called as a het.
#' @export
biallelic_geno_error <- function(
  genos,
  het_miscall = 0.02,
  hom_miscall = 0.005
) {
  L <- length(genos)
  hm <- sample(
    x = c(T, F),
    size = L,
    replace = TRUE,
    prob = c(het_miscall, 1 - het_miscall)
  )
  bm <- sample(
    x = c(T, F),
    size = L,
    replace = TRUE,
    prob = c(hom_miscall, 1 - hom_miscall)
  )
  switchy <- sample(
    x = c(T, F),
    size = L,
    replace = TRUE,
    prob = c(0.5, 0.5)
  )

  ret <- case_when(
    (genos %in% c("1/2", "2/1")) & hm == FALSE  ~ genos,
    (genos %in% c("1/1", "2/2")) & bm == FALSE  ~ genos,
    hm == TRUE & (genos == "1/2" | genos == "2/1") & switchy == TRUE ~ "1/1",
    hm == TRUE & (genos == "1/2" | genos == "2/1") & switchy == FALSE ~ "2/2",
    bm == TRUE & (genos %in% c("1/1", "2/2")) ~ "1/2",
    TRUE ~ genos
  )
}
