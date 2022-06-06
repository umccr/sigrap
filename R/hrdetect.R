#' Read VCF with SNVs/INDELs for use with HRDetect
#'
#' Reads a VCF with SNVs/INDELs for use with HRDetect.
#'
#' @param x Path to VCF.
#'
#' @return List containing CHROM, POS, REF and ALT columns
#'         for SNVs and INDELs in separate tibbles.
#'
#' @examples
#' x <- system.file("extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
#' (l <- hrdetect_read_snvindel_vcf(x))
#' @testexamples
#' expect_equal(length(l), 2)
#' expect_equal(names(l), c("snv", "indel"))
#' expect_equal(colnames(l$snv), c("chr", "position", "REF", "ALT"))
#' expect_equal(colnames(l$indel), c("chr", "position", "REF", "ALT"))
#'
#' @export
hrdetect_read_snvindel_vcf <- function(x) {
  assertthat::assert_that(file.exists(x))
  ALLOWED_BASES <- c("A", "C", "G", "T")
  readr::local_edition(1)
  d <- x |>
    readr::read_tsv(
      comment = "##",
      col_types = readr::cols_only("#CHROM" = "c", "POS" = "i", "REF" = "c", "ALT" = "c")
    ) |>
    dplyr::mutate(vartype = dplyr::case_when(
      .data$REF %in% ALLOWED_BASES & .data$ALT %in% ALLOWED_BASES ~ "SNV",
      TRUE ~ "INDEL"
    )) |>
    dplyr::rename(chr = "#CHROM", position = "POS")

  snv <- d |>
    dplyr::filter(.data$vartype == "SNV") |>
    dplyr::select(-.data$vartype)
  indel <- d |>
    dplyr::filter(.data$vartype == "INDEL") |>
    dplyr::select(-.data$vartype)

  list(
    snv = snv,
    indel = indel
  )
}

#' Read VCF with SVs for use with HRDetect
#'
#' Reads a VCF with SVs for use with HRDetect.
#'
#' @param x Path to VCF.
#' @param nm Sample name.
#' @param genome Human genome version (default: hg38. hg19 means GRCh37).
#'
#' @return Tibble with following BEDPE-like columns:
#' - chrom1, start1, end1
#' - chrom2, start2, end2
#' - sample
#' - strand1, strand2
#'
#'
#' @examples
#' x <- system.file("extdata/umccrise/sv/manta.vcf.gz", package = "gpgr")
#' sv_bedpe <- hrdetect_read_sv_vcf(x, nm = "SAMPLE")
#' head(sv_bedpe)
#' @testexamples
#' expect_equal(nrow(sv_bedpe), 190)
#' expect_equal(colnames(sv_bedpe), c("chrom1", "start1", "end1", "chrom2",
#'              "start2", "end2", "sample", "strand1", "strand2"))
#' @export
hrdetect_read_sv_vcf <- function(x, nm = NULL, genome = "hg38") {
  assertthat::assert_that(file.exists(x))
  assertthat::assert_that(!is.null(nm))
  assertthat::assert_that(genome %in% c("hg19", "hg38", "GRCh37"))
  if (genome == "GRCh37") {
    genome <- "hg19"
  }

  vcf <- VariantAnnotation::readVcf(x, genome)
  gr <- StructuralVariantAnnotation::breakpointRanges(vcf)
  bedpe <- StructuralVariantAnnotation::breakpointgr2bedpe(gr) |>
    dplyr::mutate_if(is.factor, as.character) |>
    dplyr::mutate(sample = nm) |>
    dplyr::select(
      .data$chrom1, .data$start1, .data$end1,
      .data$chrom2, .data$start2, .data$end2,
      .data$sample, .data$strand1, .data$strand2
    )

  bedpe
}

#' Read PURPLE Somatic CNVs for HRDetect
#'
#' Reads PURPLE somatic CNVs for use with HRDetect.
#'
#' @param x Path to `purple.cnv.somatic.tsv` file.
#'
#' @return Tibble containing following columns:
#' - chromosome, start, end
#' - copyNumber (total)
#' - minorAllelePloidy
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
#' (cnv <- hrdetect_read_purple_cnv(x))
#' @testexamples
#' expect_equal(colnames(cnv), c("Chromosome", "chromStart", "chromEnd",
#'                               "total.copy.number.inTumour",
#'                               "minor.copy.number.inTumour"))
#'
#' @export
hrdetect_read_purple_cnv <- function(x) {
  assertthat::assert_that(file.exists(x))
  cnv <- readr::read_tsv(x, col_types = readr::cols_only(
    "chromosome" = "c", "start" = "i", "end" = "i",
    "copyNumber" = "d", "minorAlleleCopyNumber" = "d"
  ))

  cnv |>
    dplyr::rename(
      Chromosome = .data$chromosome,
      chromStart = .data$start,
      chromEnd = .data$end,
      total.copy.number.inTumour = .data$copyNumber,
      minor.copy.number.inTumour = .data$minorAlleleCopyNumber,
    ) |>
    dplyr::mutate(Chromosome = sub("chr", "", .data$Chromosome))
}

#' Prepare VCF with SNVs/INDELs for use with HRDetect
#'
#' Prepares VCF with SNVs/INDELs for use with HRDetect.
#'
#' @param x Path to VCF with SNVs and INDELs.
#' @param nm Sample name.
#' @param genome Human genome version (default: hg38. hg19 means GRCh37).
#' @param sigsToUse COSMIC signatures to use.
#'
#' @return List with two elements:
#' - snv_results: tibble with exposure score and p-value for chosen signatures.
#' - indel_results: tibble with a summary of the count of indels and their proportion.
#'
#' @examples
#' x <- system.file("extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
#' (l <- hrdetect_prep_snvindel(x, nm = "sampleA"))
#' @testexamples
#' expect_equal(c("snv_results", "indel_results"), names(l))
#' expect_equal(c("sig", "exposure", "pvalue"), colnames(l[["snv_results"]]))
#' expect_equal(colnames(l[["indel_results"]])[c(1, 7)], c("sample", "del.mh.prop"))
#'
#' @export
hrdetect_prep_snvindel <- function(x, nm = NULL, genome = "hg38",
                                   sigsToUse = c(1, 2, 3, 5, 6, 8, 13, 17, 18, 20, 26, 30)) {
  assertthat::assert_that(
    file.exists(x), !is.null(nm),
    all(sigsToUse %in% 1:30), all(c(3, 8) %in% sigsToUse)
  )
  assertthat::assert_that(genome %in% c("hg19", "hg38", "GRCh37"))
  if (genome == "GRCh37") {
    genome <- "hg19"
  }

  snvindel_tabs <- hrdetect_read_snvindel_vcf(x)

  ## --- SNVs ---##
  snv_catalogue <- signature.tools.lib::tabToSNVcatalogue(
    subs = snvindel_tabs[["snv"]],
    genome.v = genome
  )[["catalogue"]]
  assertthat::assert_that(inherits(snv_catalogue, "data.frame"),
                          ncol(snv_catalogue) == 1,
                          colnames(snv_catalogue) == "catalogue",
                          nrow(snv_catalogue) == 96)

  subs_fit_res <- signature.tools.lib::Fit(
    catalogues = snv_catalogue,
    signatures = signature.tools.lib::COSMIC30_subs_signatures[, sigsToUse],
    useBootstrap = TRUE,
    nboot = 200,
    nparallel = 2
  )

  assertthat::assert_that(
    length(subs_fit_res) == 16,
    "exposures" %in% names(subs_fit_res))

  snv_exp <- subs_fit_res[["exposures"]] |>
    t() |>
    tibble::as_tibble(rownames = "sig") |>
    dplyr::rename(exposure = .data$catalogue)

  ## --- INDELs ---##
  indel_count_proportion <- signature.tools.lib::tabToIndelsClassification(
    indel.data = snvindel_tabs[["indel"]],
    sampleID = nm,
    genome.v = genome
  )[["count_proportion"]] |>
    tibble::as_tibble()

  list(
    snv_results = snv_exp,
    indel_results = indel_count_proportion
  )
}

#' Prepare VCF with SVs for use with HRDetect
#'
#' Prepares VCF with SVs for use with HRDetect.
#'
#' @param x Path to VCF with SVs.
#' @param nm Sample name.
#' @param genome Human genome version (default: hg38. hg19 means GRCh37).
#'
#' @return Single-column data.frame (with rownames) with counts for each SV category.
#'
#' @examples
#' x <- system.file("extdata/umccrise/sv/manta.vcf.gz", package = "gpgr")
#' nm <- "SampleA"
#' (d <- hrdetect_prep_sv(x, nm))
#' @testexamples
#' expect_equal(colnames(d), nm)
#' expect_true(inherits(d, "data.frame"))
#'
#' @export
hrdetect_prep_sv <- function(x, nm = NULL, genome = "hg38") {
  assertthat::assert_that(file.exists(x))
  if (genome == "GRCh37") {
    genome <- "hg19"
  }
  sv_bedpe <- hrdetect_read_sv_vcf(x, nm = nm, genome = genome)

  res <- signature.tools.lib::bedpeToRearrCatalogue(sv_bedpe = sv_bedpe)[["rearr_catalogue"]]
  res
}

#' Prepare PURPLE Somatic CNVs for HRDetect
#'
#' Prepares PURPLE somatic CNVs for use with HRDetect.
#'
#' @param x Path to `purple.cnv.somatic.tsv` file.
#' @param nm Sample name.
#'
#' @return Tibble with sample name and HRD-LOH index.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
#' (l <- hrdetect_prep_cnv(x, nm = "SampleA"))
#' @testexamples
#' expect_equal(colnames(l), c("name", "hrdloh_index"))
#' expect_equal(nrow(l), 1)
#'
#' @export
hrdetect_prep_cnv <- function(x, nm = NULL) {
  assertthat::assert_that(file.exists(x))
  assertthat::assert_that(!is.null(nm))
  cnv <- hrdetect_read_purple_cnv(x)
  cnv_hrd <- signature.tools.lib::ascatToHRDLOH(ascat.data = cnv, SAMPLE.ID = nm)

  tibble::tibble(
    name = names(cnv_hrd),
    hrdloh_index = unname(cnv_hrd)
  )
}

#' Run HRDetect via signature.tools.lib
#'
#' Runs HRDetect as described in the
#' [signature.tools.lib repository](https://github.com/Nik-Zainal-Group/signature.tools.lib).
#'
#' @param snvindel_vcf Path to VCF with SNVs and INDELs.
#' @param sv_vcf Path to VCF with SVs.
#' @param cnv_tsv Path to `purple.cnv.somatic.tsv` file.
#' @param nm Sample name.
#' @param genome Human genome version (default: hg38. hg19 means GRCh37).
#' @param sigsToUse COSMIC SNV signatures to use.
#' @param outpath File to write HRDetect predictions to on disk
#' (should end in '.gz'). If not specified, results won't be written to disk.
#' @return Tibble with sample name and HRD probability in first two columns.
#'
#' @examples
#' snvindel_vcf <- system.file(
#'   "extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz",
#'   package = "gpgr"
#' )
#' sv_vcf <- system.file("extdata/umccrise/sv/manta.vcf.gz", package = "gpgr")
#' cnv_tsv <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
#' nm <- "SampleA"
#' genome <- "hg38"
#' (res <- hrdetect_run(nm, snvindel_vcf, sv_vcf, cnv_tsv, genome))
#' # hrdetect_run(nm, snvindel_vcf, sv_vcf, cnv_tsv, genome,
#' #              outpath = "nogit/hrdetect_results.json.gz")
#' @testexamples
#' expect_equal(colnames(res), c("sample", "Probability", "intercept", "del.mh.prop", "SNV3",
#'                               "SV3", "SV5", "hrdloh_index", "SNV8"))
#' expect_true(inherits(res, "data.frame"))
#'
#' @export
hrdetect_run <- function(nm, snvindel_vcf, sv_vcf, cnv_tsv, genome = "hg38",
                         sigsToUse = c(1, 2, 3, 5, 6, 8, 13, 17, 18, 20, 26, 30),
                         outpath = NULL) {
  assertthat::assert_that(all(file.exists(snvindel_vcf, sv_vcf, cnv_tsv)))
  if (genome == "GRCh37") {
    genome <- "hg19"
  }

  snvindel <- hrdetect_prep_snvindel(x = snvindel_vcf, nm = nm,
                                     genome = genome, sigsToUse = sigsToUse)
  snv <- snvindel[["snv_results"]] |>
    dplyr::filter(.data$sig %in% c("Signature3", "Signature8")) |>
    tidyr::pivot_wider(
      names_from = "sig", values_from = "exposure"
    )
  # make sure this is a single-row tibble, else all hell breaks loose
  assertthat::assert_that(nrow(snv) == 1, all(colnames(snv) %in% c("Signature3", "Signature8")))
  indel <- snvindel[["indel_results"]][["del.mh.prop"]]
  sv <- hrdetect_prep_sv(x = sv_vcf, nm = nm, genome = genome)
  cnv <- hrdetect_prep_cnv(x = cnv_tsv, nm = nm)

  tib <- tibble::tibble(
    "del.mh.prop" = indel,
    "SNV3" = snv$Signature3,
    "SV3" = NA,
    "SV5" = NA,
    "hrd" = cnv$hrdloh_index,
    "SNV8" = snv$Signature8
  )
  mat <- as.matrix(tib)
  rownames(mat) <- nm
  res <- signature.tools.lib::HRDetect_pipeline(
    data_matrix = mat,
    genome.v = genome,
    SV_catalogues = sv,
    nparallel = 2
  )

  if ("hrdetect_output" %in% names(res)) {
    res <- res[["hrdetect_output"]]
  } else { # no result
    intercept <- matrix(c(NA), dimnames = list(nm, "intercept"))
    Probability <- matrix(c(NA), dimnames = list(nm, "Probability"))
    res <- cbind(intercept, mat, Probability)
  }

  res <- res |>
    tibble::as_tibble(rownames = "sample", .name_repair = "check_unique") |>
    dplyr::relocate(.data$Probability, .after = .data$sample) |>
    dplyr::rename(hrdloh_index = .data$hrd) |>
    dplyr::mutate(
      dplyr::across(tidyselect::vars_select_helpers$where(is.numeric), round, 3)
    )

  # write json.gz to file
  if (!is.null(outpath)) {
    gpgr::write_jsongz(x = res, path = outpath)
  }

  res
}
