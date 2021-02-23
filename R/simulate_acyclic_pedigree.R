# This file has a few function used to simulate
# acyclic, multigenerational pedigrees.  The approach
# is simple:
#    1. start with Nf male founders and Nm female founders.  Nf = Nm.
#    1.5 We will have N = (Nf + Nm) newborns each generation, and we will assign
#        these to the Nf possible monogamous pairs by a
#        simple multinomial + 1, i.e., each pair gets at least one offspring,
#        and the remaining offspring are divided up by a multinomial.  Unless,
#        there end up being more pairs than N because there are some founders....
#    2. randomly pair the individuals and keep them monogamous.
#    3. Let them have between 1 and X offspring, where X is 1 - Poisson(mu)
#    4. Keep track of who is in each individuals connected component, so
#       that we don't mate anyone with someone from within their connected
#       component.
#    5. If there are no unconnected mates for someone, then create a new founder
#       to be their mate.
#
#    * The first set of founders start at time 1.  Every generation after that
#      gets incremented by 1.
#    * We will name individuals as follows: S_G_I, where S is M or F for male or
#      female.  G is the generation number, and I is the index of that sex within
#      the generation.  If it is a founder, we follow it with _f.
#    * Each generation we keep a named list of the males and the females (named by their
#      IDs).  The contents of the list will be the individuals connected to each
#      of these individuals.



#' quick function to create Nm males and Nf female founders.
#' This also initializes their connected components, etc.
#' And, it returns them in a list, like what we shall use for
#' going forward...
create_founders <- function(Nm = 50, Nf = Nm) {
  male_ids <- paste0("M_1_", 1:Nm, "_f")
  female_ids <- paste0("F_1_", 1:Nm, "_f")

  male_ccs <- as.list(male_ids)
  female_ccs <- as.list(female_ids)

  names(male_ccs) <- male_ids
  names(female_ccs) <- female_ids

  ret <- list(
    males = male_ccs,
    females = female_ccs
  )

  list(ret)
}




#' Determine pairs.  Given the individuals in a given year
#' determine the pairs that will be formed, add in founders
#' as necessary to make pairings.
#' @param L a list with elements "males" and "females" which are named
#' lists of connected components.
#' @param t the year that the pairs that are being formed were born.
make_pairs <- function(L, t) {
  # get the names of the males and females. Then randomly order the male names to start
  # going through them and assigning pairs.
  m <- sample(names(L$males))
  f <- names(L$females)
  Fem <- L$females
  Mal <- L$males

  # now, cycle over males and pair each one up with an available female.
  # If no females are available, make a new founder.
  nextFem <- length(f) + 1  # the index of the next female founder we might add.
  nextMale <- length(m) + 1
  usedFem <- rep(NA, 2 * nextFem) # for dropping females that have already been paired up

  # get some variables to store the pairs.  Two parallel character vectors.
  # longer than then need to be, but the same length
  mp <- rep(NA_character_, nextMale + nextFem)
  fp <- rep(NA_character_, nextMale + nextFem)

  # also some vectors to store the IDs of the additional founders
  am <- rep(NA_character_, nextFem)
  af <- rep(NA_character_, nextMale)

  idx <- 0
  num_ff <- 0 # for counting the number of new founder females
  for(i in m) { # cycle over the males
    idx <- idx + 1
    # find eligible females (those in different connected components)
    eFem <- f[sapply(Fem, function(x) length(intersect(Mal[[i]], x)) == 0)]

    # and also remove those that have already been used:
    eFem <- setdiff(eFem, usedFem)

    # now, if, eFem is empty, make a new female founder
    if(length(eFem) == 0) {
      newf <- paste0("F_", t, "_", nextFem, "_f")
      nextFem <- nextFem + 1
      num_ff <- num_ff + 1
      af[num_ff] <- newf
    } else { # otherwise, just sample it
      newf <- sample(eFem, 1)
    }

    # and add this pair to the output vectors of pairs
    mp[idx] <- i
    fp[idx] <- newf

    # and mark that female as having been used
    usedFem[idx] <- newf


  }

  # now, if there are any females that didn't get paired up, we
  # assign them to new founder males.
  leftover_females <- setdiff(f, usedFem)
  if(length(leftover_females) > 0) {
    num_fm <- 0 # for counting the number of new founder males
    for(i in leftover_females) {

      newm <- paste0("M_", t, "_", nextMale, "_f")
      nextMale <- nextMale + 1
      num_fm <- num_fm + 1
      am[num_fm] <- newm


      # and add this pair to the output vectors of pairs
      idx <- idx + 1
      mp[idx] <- newm
      fp[idx] <- i
    }
  }

  # Finally make a tibble of the pairs that we managed to make:
  pTib <- tibble::tibble(
    pa = mp,
    ma = fp
  ) %>%
    dplyr::filter(!is.na(ma) & !is.na(pa))

  pTib
}




#' this little function takes the tibble that comes out of
#' make_pairs, and tells us how many offspring each will have.
#' For now, I am setting it to that multinom + 1, but that can
#' easily change in the future.
#' @param P the tibble that comes out of make_pairs
num_offspring <- function(P) {
  P %>%
    dplyr::mutate(num_offs = 1 + rmultinom(1, size = dplyr::n(), prob = rep(1/dplyr::n(), length.out = dplyr::n()))[,1])
    #dplyr::mutate(num_offs = 1 + rpois(n(), 1))

}


#' simple function to expand pairs into named individuals
#' with sex given randomly at 50/50
#' @param P the tibble that comes out of num_offspring
#' @param t the year that the newborns are being born
expand_offspring <- function(P, t) {
  ret <- P %>%
    mutate(
      sex_list = map(
        .x = num_offs,
        .f = function(x) sample(
          c("M", "F"),
          size = x,
          replace = TRUE
        )
      ),
      prefix = map(
        .x = sex_list,
        .f = function(x) paste0(x, "_", t, "_"))
    ) %>%
    unnest(cols = c(prefix)) %>%
    arrange(prefix, pa, ma) %>%
    mutate(sex = str_sub(prefix, 1, 1)) %>%
    group_by(sex) %>%
    mutate(offspring = paste0(prefix, 1:n())) %>%
    select(-sex_list, -prefix) %>%
    arrange(pa, ma) %>%
    ungroup()
}



#' set the initial connected component members of all newborns.
#' This will be the connected components of the parents plus
#' all of their siblings.
#' @param P a tibble like what comes out of expand_offspring
#' @param I the Indivs list component for the parent (t-1) generation
set_init_cc_members <- function(P, I) {
  P$pa_cc <- I$males[P$pa]
  P$ma_cc <- I$females[P$ma]

  # for the new founders, set their connected components to just themselves
  Pa_new_founders <- sapply(P$pa_cc, is.null)
  P$pa_cc[Pa_new_founders] <- as.list(P$pa[Pa_new_founders])

  Ma_new_founders <- sapply(P$ma_cc, is.null)
  P$ma_cc[Ma_new_founders] <- as.list(P$ma[Ma_new_founders])

  P %>%
    group_by(pa, ma) %>%
    mutate(
      sibs = list(offspring),
      init_cc = pmap(
        .l = list(
          x = pa_cc,
          y = ma_cc,
          z = sibs
          ),
        .f = function(x, y, z) unique(c(x, y, z))
      )
    ) %>%
    ungroup()
}

#' function that iteratively expands the connected component members
#' of each individual, until they reflect all the individuals they
#' are connected to.
#' @param P a tibble like that which comes out of expand_offspring
iteratively_expand_conn_comps <- function(P) {
  cc <- P$init_cc

  do_it <- 1
  while (do_it == 1) {
    old_cc <- cc
    for(i in seq_along(cc)) {
      for(j in seq_along(cc)) {
        if(length(intersect(cc[[i]], cc[[j]])) > 0) {
          the_u <- union(cc[[i]], cc[[j]])
          cc[[i]] <- the_u
          cc[[j]] <- the_u
        }
      }
    }
    # check to see if old and new are different here
    #if(all(unlist(lapply(1:length(old_cc), function(x) all(old_cc[[x]] == cc[[x]]))))) {
    if(identical(cc, old_cc)) {
      do_it <- 0
    }
  }

  P$full_cc <- cc

  # return it with a column showing how many more members were added to the
  # full_cc over the init_cc
  P %>%
    mutate(
      full_minus_init = map2(
        .x = full_cc,
        .y = init_cc,
        .f = function(x, y) setdiff(x, y)
      )
    )

}


#' this function takes the tibble coming out of
#' iteratively _expand_conn_comps, and it turn it
#' into a list item that can be returned into Indiv[[t]].
translate_tib_to_list_element <- function(Pairs) {
  # now, store the males and females with their connected component members
  # back into Indivs

  ret <- list()

  pm_tmp <- Pairs %>%
    filter(sex == "M")
  ret$males <- pm_tmp$full_cc
  names(ret$males) <- pm_tmp$offspring

  pf_tmp <- Pairs %>%
    filter(sex == "F")
  ret$females <- pf_tmp$full_cc
  names(ret$females) <- pf_tmp$offspring

  ret$pTib <- Pairs

  ret
}


#' this is how it all comes together
simulate_acyclic_pedigree <- function(Nm = 10, Nf = 10, end_time = 4) {

  stopifnot(end_time >= 2)
  Indivs <- create_founders(Nm = Nm, Nf = Nf)

  for(t in 2:end_time) {

    Pairs <- make_pairs(Indivs[[t-1]], t-1) %>%
      num_offspring() %>%
      expand_offspring(t = t) %>%
      set_init_cc_members(I = Indivs[[t - 1]]) %>%
      iteratively_expand_conn_comps()

    Indivs[[t]] <- translate_tib_to_list_element(Pairs)

  }

  # now, return all the stuff.  With the clean pedigree being the first thing.
  pTibs_all <- map(Indivs[-1], "pTib") %>%
    bind_rows()

  pedigree <- pTibs_all %>%
    select(pa, ma, offspring) %>%
    rename(
      Pa = pa,
      Ma = ma,
      Kid = offspring
    )

  return(
    list(
      pedigree = pedigree,
      pTibs_all = pTibs_all,
      Indivs_list = Indivs
    )
  )

}
