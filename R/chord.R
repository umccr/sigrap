#' Run CHORD
#'
#' Runs CHORD given SNV and SV VCF files. **NOTE**: make sure you have
#' the BSgenome.Hsapiens.UCSC.hgXX installed.
#'
#' @param vcf.snv Path to VCF containing SNVs and INDELs.
#' @param vcf.sv Path to VCF containing SVs.
#' @param df.sv A data.frame object containing the columns
#'        'SVTYPE' and 'SVLEN' from a Manta SV VCF.
#' @param sample.name Name of sample to use.
#' @param ref.genome Human genome assembly. One of 'hg38' (default), 'hg19' or 'GRCh37'.
#' @param sv.caller manta (default) or gridss.
#' @param outpath File to write CHORD predictions to on disk
#' (should end in '.gz'). If not specified, results won't be written to disk.
#' @param ... Other arguments to be passed to
#'            <https://github.com/UMCUGenetics/CHORD/blob/d7c963/R/extractSigsChord.R>.
#'
#' @examples
#'
#' snv <- system.file("extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
#' sv <- system.file("extdata/umccrise/sv/manta.vcf.gz", package = "gpgr")
#' chord_res <- chord_run(
#'   vcf.snv = snv, df.sv = chord_mantavcf2df(sv),
#'   sample.name = "foo"
#' )
#' # chord_res2 <- chord_run(vcf.snv = snv, vcf.sv = sv, sample.name = "foo",
#' #                         outpath = "nogit/chord_results.json.gz")
#' @testexamples
#'
#' expect_equal(length(chord_res), 2)
#' expect_equal(names(chord_res), c("contexts", "prediction"))
#' expect_equal(chord_res[["prediction"]][, 1, drop = TRUE], "foo")
#'
#'
#' @return List with extracted signatures and HRD prediction.
#'
#' @export
chord_run <- function(vcf.snv = NULL, vcf.sv = NULL, df.sv = NULL,
                      sample.name = NULL, ref.genome = "hg38",
                      sv.caller = "manta", outpath = NULL, ...) {
  g <- get_genome_obj(ref.genome)

  contexts <- CHORD::extractSigsChord(
    vcf.snv = vcf.snv,
    vcf.sv = vcf.sv,
    df.sv = df.sv,
    sv.caller = sv.caller,
    sample.name = sample.name,
    ref.genome = g,
    ...
  )

  prediction <-
    CHORD::chordPredict(
      features = contexts,
      rf.model = CHORD::CHORD,
      do.bootstrap = TRUE, verbose = FALSE
    ) |>
    dplyr::mutate(
      dplyr::across(tidyselect::vars_select_helpers$where(is.numeric), round, 3)
    )

  # custom order of prediction cols
  col_order <- c(
    "sample", "p_hrd", "hr_status",
    "hrd_type", "p_BRCA1", "p_BRCA2",
    "remarks_hr_status", "remarks_hrd_type",
    "p_hrd.5%", "p_hrd.50%", "p_hrd.95%",
    "p_BRCA1.5%", "p_BRCA1.50%", "p_BRCA1.95%",
    "p_BRCA2.5%", "p_BRCA2.50%", "p_BRCA2.95%"
  )

  assertthat::assert_that(all(colnames(prediction) %in% col_order))
  prediction <- prediction[col_order]

  # write json.gz to file
  if (!is.null(outpath)) {
    gpgr::write_jsongz(x = prediction, path = outpath)
  }

  list(
    contexts = contexts,
    prediction = prediction
  )
}

#' Convert Manta VCF to data.frame
#'
#' Converts a Manta VCF to a data.frame for processing with CHORD.
#'
#' @param in_vcf Manta VCF.
#'
#' @return Tibble with two columns: sv_type and sv_len (INFO/SVTYPE and INFO/SVLEN from VCF).
#'
#' @examples
#' in_vcf <- system.file("extdata/umccrise/sv/manta.vcf.gz", package = "gpgr")
#' d <- chord_mantavcf2df(in_vcf)
#' @testexamples
#' expect_equal(d$sv_len[1], "-108")
#' expect_equal(d$sv_type[1], "DEL")
#'
#' @export
chord_mantavcf2df <- function(in_vcf) {
  assertthat::assert_that(file.exists(in_vcf))
  assertthat::assert_that(grepl("vcf$", in_vcf) | grepl("vcf\\.gz$", in_vcf))
  d <- bedr::read.vcf(in_vcf, split.info = TRUE, verbose = FALSE)
  tibble::tibble(
    sv_type = d$vcf$SVTYPE,
    sv_len = d$vcf$SVLEN
  )
}

#' Get BSgenome Object
#'
#' Gets BSgenome object given a human genome assembly string.
#'
#' @param genome Human genome assembly string: hg38 (default), hg19 or GRCh37.
#'
#' @return BSgenome object.
#'
#' @examples
#' \dontrun{
#' get_genome_obj("hg38")
#' }
#' @testexamples
#' expect_error(get_genome_obj("FOO"))
#'
#' @export
get_genome_obj <- function(genome = "hg38") {
  bsgenome <- c(
    hg19 = "BSgenome.Hsapiens.UCSC.hg19",
    hg38 = "BSgenome.Hsapiens.UCSC.hg38",
    GRCh37 = "BSgenome.Hsapiens.1000genomes.hs37d5"
  )
  pkg <- bsgenome[genome]
  assertthat::assert_that(
    genome %in% names(bsgenome),
    msg = glue::glue(
      "Instead of '{genome}', pick one of: ",
      "{paste(names(bsgenome), collapse = ', ')}"
    )
  )
  if (!gpgr::pkg_exists(pkg)) {
    stop(glue::glue(
      "{pkg} is not installed on your system.\n",
      "Please install with:\n'BiocManager::install(\"{pkg}\")'"
    ))
  }
  return(eval(parse(text = glue::glue("{pkg}::{pkg}"))))
}
