

#' convert a long-format data frame of loci into a vector of lines for a Mendel def file
#'
#' This makes a vector of text that has carriage returns and things in it as necessary.
#' So, to include it in a def file, if the return from this function is called
#' D you can just do cat(D, sep = "", file = myDef.txt), for example.
#'
#' Because it appears that Mendel requires unique locus names we will name them
#' Chrom-underscore-Loc.
#'
#' @param df  A data frame in the format of \code{\link{long_markers}}.
#' @keywords internal
#' @export
#' @examples
#' data(long_markers)
#' D <- markers2mendel_def_lines(long_markers)
#' cat(D[1:100], sep = "")
markers2mendel_def_lines <- function(df) {

  t1 <- df %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Chrom, LocIdx, AlleIdx) %>%
    dplyr::group_by(Chrom, LocIdx)

  # form all the locus names and make sure none of them are more than 16 characters
  illegal_names <- t1 %>%
    dplyr::summarise(loc_name = paste0(Chrom, "_", LocIdx)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(nchar(loc_name) > 16)

  if(nrow(illegal_names) > 0) {
    top10 <- paste(unique(illegal_names$loc_name)[1:10], collapse = ", ")
    outstr <- paste("Names formed as Chrom_LocusIndex to pass to Mendel are > 16 characters.\n\nExamples: ", top10, "\n\nPlease give your chromosomes shorter names.\n\n", collapse = "")
    stop(outstr)
  }

  tmp <- t1 %>%
    dplyr::mutate(LocLine = paste(Chrom, "_", LocIdx, ",,", n_distinct(AlleIdx), sep = ""),
                  LocScrub = ifelse(AlleIdx == 1, paste(LocLine,"\n"), ""),
                  AlleLine = paste(AlleIdx, ",", Freq, "\n", sep = ""),
                  TextVec = paste(LocScrub, AlleLine, sep = "   ")
    )

  # Note, in the LocLine above we leave the chromosome name blank (the two commas next to one another)
  # so that Mendel knows that it is autosomal.

  tmp$TextVec
}


#' create a vector of Mendel formatted map file strings
#'
#' This makes a vector of text that has carriage returns and things in it as necessary.
#' So, to include it in a map file, if the return from this function is called
#' D you can just do cat(D, sep = "", file = myMap.txt), for example.
#'
#' @param df  A data frame in the format of \code{\link{long_markers}}.
#' @keywords internal
#' @export
#' @examples
#' data(long_markers)
#' D <- markers2mendel_map_lines(long_markers)
#' cat(D[1:100], sep = "")
markers2mendel_map_lines <- function(df) {
  tmp <- df %>%
    dplyr::arrange(Chrom, LocIdx, AlleIdx) %>%
    dplyr::group_by(Chrom) %>%
    dplyr::mutate(ender = ifelse(LocIdx == max(LocIdx), "\n\n", "\n")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Chrom, LocIdx) %>%
    dplyr::summarise(map_lines = paste(Chrom, "_", LocIdx, ", ", Pos, ender , sep = "")[1])

  tmp$map_lines
}


#' for making a pedigree file with no observed genotypes on it
#' @keywords internal
pedigree2mendel_ped_file <- function(df, name, filename = NA) {


  df2 <- df %>%
    tibble::as_tibble() %>%
    dplyr::mutate(pedname = name,
                  MZtwin_stat = "",
                  Panew = ifelse(Pa == 0, "", Pa),
                  Manew = ifelse(Ma == 0, "", Ma),
                  Sexnew = ifelse(Sex == 0, "", c("M", "F")[Sex]),
                  TwinStatus = ""
    ) %>%
    dplyr::select(pedname, Kid, Panew, Manew, Sexnew, TwinStatus)

  if(is.na(filename)) {
    return(df2)
  }

  write.table(df2, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

  invisible(NULL)
}


#' prepend an ID to a mendel input or output file
#'
#' This is just a little helper fucntion to ensure consistency in file naming.
#' @param ID The ID to be prepended.
#' @param suff The final part of the file name (like Def.in, etc.)
id_prepend <- function(ID, suff) {
  paste(ID, "-", suff, sep = "")
}


#' return a list of values for a Mendel definitions file for gene dropping
#'
#' @param ID the name/prefix to be given to all the files involved in this run.
#' @param Reps the number of replicate simulations to perform
#' @param Seed an integer number to be used as the seed for the simulations.
mendel_control_list <- function(ID, Reps, Seed) {
  ret <- list(
    # input files
    DEFINITION_FILE = id_prepend(ID, "Def.in"),
    MAP_FILE = id_prepend(ID, "Map.in"),
    PEDIGREE_FILE = id_prepend(ID, "Ped.in"),

    # output files
    NEW_PEDIGREE_FILE = id_prepend(ID, "Ped.out"),
    OUTPUT_FILE = id_prepend(ID, "Mendel.out"),
    SUMMARY_FILE = id_prepend(ID, "Summary.out"),

    # analysis options
    ANALYSIS_OPTION = "Gene_dropping",
    REPETITIONS = sprintf("%d", unname(Reps)),
    SEED = Seed,
    MODEL = 2,
    KEEP_FOUNDER_GENOTYPES = "False",
    MISSING_DATA_PATTERN = "None",
    GENE_DROP_OUTPUT = "Unordered",
    MAP_DISTANCE_UNITS = "bp"
  )

  ret
}




#' Write all necessary files for Mendel to do a gene-dropping run
#'
#' Writes all files in the directory Dir but does not launch mendel.
#' @param ID  The ID or identifier of the run.  This string will be prepended to all input and output files.
#' @param Reps The number of replicate data sets to simulate.
#' @param Seed The random seed to use.  This must be a single integer between 1 and 30,000.
#' @param Markers A data frame in the format of \code{\link{long_markers}} that gives information
#' about marker positions and allele frequencies.
#' @param Pedigree A data frame formatted like any component of \code{\link{pedigrees}} which holds the
#' pedigree to drop genes upon.
#' @param Dir The directory in which to write the files.  It must be created already.  Defaults to the
#' current working directory.
#' @keywords internal
#' @export
#' @examples
#' # make a temp directory to place things into:
#' tmpDir <- tempdir()
#' write_all_mendel_files("mendel-example", 10, 1234, long_markers, pedigrees$HS, Dir = tmpDir)
write_all_mendel_files <- function(ID, Reps, Seed, Markers, Pedigree, Dir = ".") {
  # get the lines
  DL <- markers2mendel_def_lines(Markers)
  ML <- markers2mendel_map_lines(Markers)
  MCL <- mendel_control_list(ID, Reps, Seed)

  # write the pedigree file
  pedigree2mendel_ped_file(df = Pedigree, name = "xyz", filename = file.path(Dir, id_prepend(ID, "Ped.in")))

  # write the control file
  cat(paste(names(MCL), MCL, sep = " = "), sep = "\n", file = file.path(Dir, id_prepend(ID, "Control.in")))

  # write the definitions file
  cat(DL, sep = "", file = file.path(Dir, id_prepend(ID, "Def.in")))

  # write the map file
  cat(ML, sep = "", file = file.path(Dir, id_prepend(ID, "Map.in")))
}


#' Given the operating system, check for mendel program in the typical location
#'
mendelBin <- function() {

  sys_name <- Sys.info()['sysname']

  if(sys_name == "Darwin") {  # we are on a Mac
    path <- "/Applications/Mendel/mendel"
  } else if(sys_name == "Windows") { # we are on a PC
    path <- "C:/Program Files (x86)/Mendel-16/Mendel.exe"
  } else {
    stop("Sorry, we don't know where to expect the mendel binary on an operating system of type ", sys_name, ". Send email to eric.anderson@noaa.gov and he can work with you to get things configured.")
  }
  if(!file.exists(path)) stop("Didn't find mendel binary where we expected it at ", path, " Please download, then install, MENDEL v16.0 from http://software.genetics.ucla.edu/")
  path
}

#' Run the mendel binary in directory Dir using control file Control
#'
#' @param Dir the directory in which to run mendel
#' @param Control the "control file" with which to run mendel.  If not an absolute path
#' it should be given relative to the Dir directory.
#' @examples
#' \dontrun{
#' # first prepare the input files and define tmpDir
#' example(write_all_mendel_files)
#'
#' # then run mendel
#' run_mendel(tmpDir, "mendel-example-Control.in")
#' }
run_mendel <- function(Dir, Control) {

  nowdir <- getwd();
  setwd(Dir)

  system2(command = mendelBin(),
          args = paste("-c", Control)
  )

  setwd(nowdir)

  # COMM <- paste("cd",
  #               Dir,
  #               ";",
  #               mendelBin(),
  #               "-c",
  #               Control
  #               )
  # system(COMM)
}




#' simulate genotype-pairs from linked markers using Mendel
#'
#' This is the main function to use.  Pass it a data frame of markers (indexed in order)
#' and it will pass them off to Mendel, drop the genes, do genotyping error and then return
#' a list of vectors which hold the genotype index that you can use to subscript the
#' joint probability vectors with.  This function assumes that the loci in df are ordered
#' appropriately (i.e have been run through reindex_markers()) and that the components in
#' C are named Chrom.Locus.Pos, as is typical.  Obviously the C list should corresponds
#' exactly to the markers/alleles in df.
#' @param df  A data frame of markers with the columns Chrom, Locus, Pos, Allele,
#' LocIdx, AlleIdx, and Freq.  This should have been run through `reindex_markers()`.
#' @param ped  The pedigree to be simulating from
#' @return  A wide-format tibble of genotypes.  Each genotype occupies one
#' column with the index of each allele separated by a slash, like: "1/2" or
#' "1/1".
#' @export
simulate_linked_genotypes <- function(df, ped) {

  # this is just hardwired to do a single rep.  We could probably change
  # that later, if desired, but for our present applications, I think that
  # is OK>
  num <- 1

  # first get the number of alleles at each locus
  alle_nums <- df %>%
    dplyr::mutate(list_name = paste(Chrom, Locus, Pos, sep = ".")) %>%
    dplyr::mutate(list_name = factor(list_name, levels = unique(list_name))) %>%
    dplyr::group_by(list_name) %>%
    dplyr::tally()

  # this runs everything through Mendel and delivers outgenos, a list of two matrices. (One matrix for indiv1 and the other for indiv2)
  # In each, rows are reps and columns are loci and the entries are the indexes of genotypes.
  tmpDir = tempfile()  # make this run in a separate directory each time in case parallelizing...
  dir.create(tmpDir, recursive = TRUE)
  message("Mendel running in temp directory ", tmpDir)
  write_all_mendel_files("mendel-example", num, floor(runif(1, min = 100, max = 100000)), df, ped, Dir = tmpDir)
  run_mendel(tmpDir, "mendel-example-Control.in")
  sed_call <- paste("cd ", tmpDir, "; sed 's/ *//g;'  mendel-example-Ped.out > ped_out.csv;", sep = "")
  system(sed_call)

  # now, try to read that in.  It is quick with vroom.
  big_frame <- vroom::vroom(
    file = file.path(tmpDir, "ped_out.csv"),
    delim = ",",
    col_names = FALSE,
    col_types = c(.default = "c")
  ) %>%
    dplyr::select(-X1, -X3, -X4, -X6) %>%
    dplyr::rename(
      indiv = X2,
      sex = X5
    )

  # now, put the locus names on there
  loc_names <- unique(df$Locus)
  big_frame <- big_frame[, 1:(length(loc_names) + 2)]
  names(big_frame)[-(1:2)] <- loc_names

  big_frame
}
