#' Get COSMIC 2015 Signatures
#'
#' @usage data(cosmic_signatures_2015)
#' @docType data
#'
#' @format A matrix with 96 rows (one for each somatic mutation type within a specific context)
#' and 30 columns (one for each signature).
#'
"cosmic_signatures_2015"


#' Calculate COSMIC Signature Contribution
#'
#' Finds the linear combination of COSMIC (2015 and 2020) mutation signatures that
#' most closely reconstructs the SNV mutation matrix by solving the
#' nonnegative least-squares constraints problem.
#'
#' @param mut_mat Mutation count matrix (dimensions: m rows (mutation types) X 1 column (sample)).
#' @param signatures Signature matrix (dimensions: m rows (mutation types) X n columns (signatures))
#'
#' @return A list with the COSMIC 2015 and 2020 signature contributions to the
#' sample's signature.
#'
#' @export
sig_contribution <- function(mut_mat, signatures) {
  # Fit mutation matrix to cancer signatures
  fit_res <-
    MutationalPatterns::fit_to_signatures(mut_mat, signatures)$contribution |>
    tibble::as_tibble(rownames = "sig") |>
    dplyr::rename(contr = 2) |>
    dplyr::filter(.data$contr > 0)

  if (nrow(fit_res) == 0) {
    fit_res <- tibble::tribble(
      ~sig, ~contr,
      "No Signatures found!", 0
    )
  }

  fit_res_contr <- fit_res |>
    dplyr::mutate(
      contr = round(.data$contr, 0),
      RelFreq = round(.data$contr / sum(.data$contr), 2),
      Rank = as.integer(base::rank(-.data$contr))
    ) |>
    dplyr::select("Rank",
      Signature = "sig",
      Contribution = "contr", "RelFreq"
    ) |>
    dplyr::arrange(.data$Rank)

  fit_res_contr
}

#' Create COSMIC Signature Table
#'
#' Creates a COSMIC signature table for display in an RMarkdown report.
#'
#' @param contr A tibble with Rank, Signature, Contribution and RelFreq columns.
#' @param type One of Sig (old COSMIC), SBS, DBS or ID.
#' @param outdir Relative directory to write signature plots to for incorporating
#' in an RMarkdown report.
#'
#' @return The `contr` tibble with an additional `Plot` column pointing to
#' the local path to the corresponding signature plot (in markdown syntax).
#' If outdir = NULL (default), simply return the `contr` tibble with the
#' description of the signature.
#'
#' @export
sig_contribution_table <- function(contr, type, outdir = NULL) {
  available_types <- c(
    "Sig" = "v2_2015/Sig",
    "SBS" = "v3.2_2021-march/SBS",
    "DBS" = "v3.2_2021-march/DBS",
    "ID" = "v3.2_2021-march/ID"
  )
  assertthat::assert_that(length(type) == 1, type %in% names(available_types))
  assertthat::assert_that(all(colnames(contr) == c("Rank", "Signature", "Contribution", "RelFreq")))

  sig_dir <- system.file(file.path("extdata/sigs", available_types[type]), package = "sigrap")

  # if an outdir is specified, copy out the COSMIC plots to be linked to
  if (!is.null(outdir)) {
    img_cp_dir <- file.path(outdir, "sig_plots", type)
    fs::dir_create(img_cp_dir)
    sig_table <-
      readr::read_tsv(file = file.path(sig_dir, "description.tsv"), col_types = "cc") |>
      dplyr::mutate(
        Plot_original = file.path(sig_dir, paste0(.data$signature, ".png")),
        Plot_copy = file.path(img_cp_dir, paste0(.data$signature, ".png")),
        signature = paste0(type, .data$signature)
      ) |>
      dplyr::rename(Signature = .data$signature)

    d <- contr |>
      dplyr::left_join(sig_table, by = "Signature") |>
      dplyr::rowwise() |>
      dplyr::mutate(cp = file.copy(from = .data$Plot_original, to = .data$Plot_copy)) |>
      dplyr::mutate(plot_in_md = paste0("![](", .data$Plot_copy, ")")) |>
      dplyr::select(
        "Rank", "Signature", "Contribution", "RelFreq",
        Description = "description", Plot = "plot_in_md"
      )
  } else { # don't copy out any COSMIC plots
    sig_table <-
      readr::read_tsv(file = file.path(sig_dir, "description.tsv"), col_types = "cc") |>
      dplyr::mutate(
        signature = paste0(type, .data$signature)
      ) |>
      dplyr::rename(Signature = .data$signature)
    d <- contr |>
      dplyr::left_join(sig_table, by = "Signature") |>
      dplyr::rename(Description = .data$description)
  }
  return(d)
}

#' Count SNV Contexts
#'
#' Counts SNV Contexts.
#'
#' @param vcf_gr GRanges containing all mutation types from a single sample.
#' @param ref_genome The BSGenome reference genome to use.
#'
#' @return A list with two elements:
#' - snv_counts: matrix containing the number of SNVs per COSMIC context per gr.
#' - gr_snv: GRanges object containing the SNVs.
#'
#' @export
sig_count_snv <- function(vcf_gr, ref_genome) {
  cli::cli_h2(glue::glue("{date_log()} Counting SNV contexts"))
  gr_snv <- MutationalPatterns::get_mut_type(vcf_list = vcf_gr, type = "snv")
  snv_counts <- MutationalPatterns::mut_matrix(vcf_list = gr_snv, ref_genome = ref_genome)
  list(
    snv_counts = snv_counts,
    gr_snv = gr_snv
  )
}

#' Plot SNV Mutation Characteristics
#'
#' Plots SNV mutation characteristics.
#'
#' @param gr_snv GRanges containing SNVs from a single sample.
#' @param snv_counts A matrix with counts of SNV contexts.
#' @param ref_genome The BSGenome reference genome to use.
#' @param rainfall Logical. Whether to generate rainfall plot. Default is FALSE.
#'
#' @return A list with up to five ggplot2 objects:
#' - p_heatmap: a SNV mutation matrix as a heatmap.
#'   This is especially usefull when looking at a wide mutational context.
#' - p_river: a SNV mutation matrix as a riverplot.
#'   This is especially usefull when looking at a wide mutational context.
#' - p_96_profile: relative contribution of 96 trinucleotides.
#' - p_spectrum: point mutation spectrum.
#' - p_rainfall: a rainfall plot showing intermutational distances (if rainfall = TRUE).
#'
#' @export
sig_plot_snv <- function(gr_snv, snv_counts, ref_genome, rainfall = FALSE) {
  
  mut_to <- MutationalPatterns::mut_type_occurrences(
    vcf_list = gr_snv, ref_genome = ref_genome
  )

  mut_mat_ext_context <- MutationalPatterns::mut_matrix(
    vcf_list = gr_snv, ref_genome = ref_genome, extension = 2
  )

  p_spectrum <- MutationalPatterns::plot_spectrum(
    type_occurrences = mut_to, CT = TRUE,
    condensed = TRUE, error_bars = "none"
  ) +
    ggplot2::theme(legend.position = "top")
  p_96_profile <- MutationalPatterns::plot_96_profile(mut_matrix = snv_counts, condensed = TRUE)
  p_heatmap <- MutationalPatterns::plot_profile_heatmap(mut_matrix = mut_mat_ext_context) +
    ggplot2::theme(legend.position = "none")
  p_river <- MutationalPatterns::plot_river(mut_matrix = mut_mat_ext_context) +
    ggplot2::theme(legend.position = "none")
  
  if (rainfall) {
    p_rainfall <- sigrap::sig_plot_rainfall(vcf_gr = gr_snv, ref_genome = ref_genome)
  }

  result <- list(
    p_heatmap = p_heatmap,
    p_river = p_river,
    p_96_profile = p_96_profile,
    p_spectrum = p_spectrum
  )
  
  if (rainfall) {
    result$p_rainfall <- p_rainfall
  }
  
  return(result)
}

#' Count INDEL Contexts
#'
#' Counts INDEL Contexts.
#'
#' @param vcf_gr GRanges containing all mutation types from a single sample.
#' @param ref_genome The BSGenome reference genome to use.
#'
#' @return A tibble containing the number of INDELs per COSMIC context per gr.
#'
#' @export
sig_count_indel <- function(vcf_gr, ref_genome) {
  cli::cli_h2(glue::glue("{date_log()} Counting INDEL contexts"))
  gr_indel <- MutationalPatterns::get_mut_type(vcf_list = vcf_gr, type = "indel")
  gr_indel <- MutationalPatterns::get_indel_context(vcf_list = gr_indel, ref_genome = ref_genome)
  indel_counts <- MutationalPatterns::count_indel_contexts(vcf_list = gr_indel)
  indel_counts
}


#' Plot INDEL Mutation Characteristics
#'
#' Plots INDEL mutation characteristics.
#'
#' @param indel_counts INDEL context counts.
#'
#' @return A list with two ggplot2 objects:
#' - p_indel_main: the main INDEL contexts.
#' - p_indel_cont: the INDEL contexts.
#'
#' @export
sig_plot_indel <- function(indel_counts) {
  p_indel_main <- MutationalPatterns::plot_main_indel_contexts(counts = indel_counts) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.9, hjust = 1)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  p_indel_cont <- MutationalPatterns::plot_indel_contexts(counts = indel_counts, condensed = TRUE) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.9, hjust = 1),
      legend.position = "top"
    )

  list(
    p_indel_main = p_indel_main,
    p_indel_cont = p_indel_cont
  )
}

#' Count DBS Contexts
#'
#' Counts DBS Contexts.
#'
#' @param vcf_gr GRanges containing all mutation types from a single sample.
#'
#' @return A tibble containing the number of DBS per COSMIC context per gr.
#'
#' @export
sig_count_dbs <- function(vcf_gr) {
  cli::cli_h2(glue::glue("{date_log()} Counting DBS contexts"))
  gr_dbs <- MutationalPatterns::get_mut_type(vcf_list = vcf_gr, type = "dbs", predefined_dbs_mbs = TRUE)
  gr_dbs <- MutationalPatterns::get_dbs_context(vcf_list = gr_dbs)
  dbs_counts <- MutationalPatterns::count_dbs_contexts(vcf_list = gr_dbs)
  dbs_counts
}


#' Count mbs Contexts
#'
#' Counts mbs (Multi-Base Substitution) Contexts.
#'
#' @param vcf_gr GRanges containing all mutation types from a single sample.
#'
#' @return A matrix containing the number of MBSs per COSMIC context per gr.
#'
#' @export
sig_count_mbs <- function(vcf_gr) {
  cli::cli_h2(glue::glue("{date_log()} Counting MBS contexts"))
  gr_mbs <- MutationalPatterns::get_mut_type(vcf_list = vcf_gr, type = "mbs", predefined_dbs_mbs = TRUE)
  mbs_counts <- MutationalPatterns::count_mbs_contexts(gr_mbs)  # counts by length (3, 4, 5+)
  mbs_counts
}


#' Plot MBS Mutation Characteristics
#'
#' Plots MBS (Multi-Base Substitution) mutation characteristics.
#'
#' @param mbs_counts MBS context counts (matrix with counts by length).
#' @param same_y Logical. If TRUE, all facets have the same y-axis scale. Default is TRUE.
#'
#' @return A ggplot2 object showing MBS contexts.
#'
#' @export
sig_plot_mbs <- function(mbs_counts, same_y = TRUE) {
  p_mbs <- MutationalPatterns::plot_mbs_contexts(counts = mbs_counts, same_y = same_y)
  
  list(
    p_mbs = p_mbs
  )
}


#' Create MBS Table for JSON Output
#'
#' Creates a table of MBS counts by variant length (3bp, 4bp, 5+bp) for JSON export.
#' Unlike other signature types, MBS counts are simple counts categorized by length 
#' rather than fitted signature patterns.
#'
#' @param mbs_counts MBS context counts (matrix with counts by length).
#'
#' @return A tibble with Type (MBS length category) and Count columns.
#'
#' @export
sig_mbs_table <- function(mbs_counts) {
  # Sum across samples if multiple columns (though typically should be single sample)
  if (is.matrix(mbs_counts)) {
    total_counts <- rowSums(mbs_counts)
  } else {
    total_counts <- mbs_counts
  }
  
  # Create simple tibble with counts by MBS type
  mbs_table <- tibble::tibble(
    Type = names(total_counts),
    Count = as.integer(total_counts)
  ) |>
    dplyr::filter(.data$Count > 0) |>
    dplyr::arrange(dplyr::desc(.data$Count))
  
  mbs_table
}


#' Plot DBS Mutation Characteristics
#'
#' Plots DBS mutation characteristics.
#'
#' @param dbs_counts DBS context counts.
#'
#' @return A list with two ggplot2 objects:
#' - p_dbs_main: the main DBS contexts.
#' - p_dbs_cont: the DBS contexts.
#'
#' @export
sig_plot_dbs <- function(dbs_counts) {
  p_dbs_main <- MutationalPatterns::plot_main_dbs_contexts(counts = dbs_counts)
  p_dbs_cont <- MutationalPatterns::plot_dbs_contexts(counts = dbs_counts, condensed = TRUE)

  list(
    p_dbs_main = p_dbs_main,
    p_dbs_cont = p_dbs_cont
  )
}

#' Run MutationalPatterns Workflow
#'
#' Runs a MutationalPatterns Workflow.
#'
#' Writes plots and signature contribution tables to `outdir`.
#'
#' @param vcf VCF file with SNVs.
#' @param sample_nm Sample name.
#' @param ref_genome Human genome assembly string: hg38 (default), hg19 or GRCh37.
#' @param outdir Directory path to output all the results to.
#' @param rainfall Logical. Whether to generate rainfall plot. Default is FALSE.
#' @param strand_bias Logical. Whether to generate strand bias analysis. Default is FALSE.
#' @param predefined_dbs_mbs Logical. Whether DBS/MBS variants are predefined in the VCF. Default is TRUE.
#'
#' @export
sig_workflow_run <- function(vcf, sample_nm, ref_genome = "hg38", outdir, rainfall = FALSE, strand_bias = FALSE) {
  fs::dir_create(outdir)
  outdir <- normalizePath(outdir)
  ref_genome <- get_genome_obj(ref_genome)

  save_plot_list <- function(pl, outdir) {
    assertthat::assert_that(is.list(pl), length(pl) > 0)
    assertthat::assert_that(all(sapply(pl, inherits, "ggplot")))
    fs::dir_create(outdir)
    for (i in seq_len(length(pl))) {
      nm <- names(pl)[i]
      fn <- file.path(outdir, paste0(nm, ".png"))
      plot_obj <- pl[[i]]
      
      # Use larger dimensions for rainfall plots
      if (nm == "p_rainfall") {
        ggplot2::ggsave(filename = fn, plot = plot_obj, width = 25, height = 8, units = "in", dpi = 300)
      } else {
        ggplot2::ggsave(filename = fn, plot = plot_obj, width = 10, height = 6, units = "in", dpi = 300)
      }
    }
  }

  cli::cli_h2(glue::glue("{date_log()} Start of MutationalPatterns workflow"))
  # import VCF
  gr <- MutationalPatterns::read_vcfs_as_granges(
    vcf_files = vcf,
    sample_names = sample_nm,
    genome = ref_genome,
    group = "auto+sex",
    type = "all"
  )

  #---- SBS ----#
  # plots
  snv_counts <- sigrap::sig_count_snv(vcf_gr = gr, ref_genome = ref_genome)
  p_snv <- sigrap::sig_plot_snv(
    gr_snv = snv_counts$gr_snv, snv_counts = snv_counts$snv_counts,
    ref_genome = ref_genome, rainfall = rainfall
  )

  if (strand_bias) {
    p_strand <- sigrap::sig_plot_strand_bias(vcf_gr = gr, ref_genome = ref_genome)
  }

  # signature contributions (2015)
  sigs_snv_2015 <-
    sigrap::cosmic_signatures_2015 |>
    {
      \(sigs) sigrap::sig_contribution(mut_mat = snv_counts$snv_counts, signatures = sigs)
    }() |>
    sigrap::sig_contribution_table(type = "Sig")

  # signature contributions (2020)
  sigs_snv_2020 <-
    MutationalPatterns::get_known_signatures(
      muttype = "snv",
      incl_poss_artifacts = TRUE
    ) |>
    {
      \(sigs) sigrap::sig_contribution(mut_mat = snv_counts$snv_counts, signatures = sigs)
    }() |>
    sigrap::sig_contribution_table(type = "SBS")

  #---- DBS ----#
  # plots
  dbs_counts <- sigrap::sig_count_dbs(vcf_gr = gr)
  p_dbs <- sigrap::sig_plot_dbs(dbs_counts = dbs_counts)

  # signature contributions
  sigs_dbs <-
    MutationalPatterns::get_known_signatures(muttype = "dbs") |>
    {
      \(sigs) sigrap::sig_contribution(mut_mat = dbs_counts, signatures = sigs)
    }() |>
    sigrap::sig_contribution_table(type = "DBS")

  #---- MBS ----#
  # counts and plots
  mbs_counts <- sigrap::sig_count_mbs(vcf_gr = gr)
  p_mbs <- sigrap::sig_plot_mbs(mbs_counts = mbs_counts, same_y = TRUE)
  
  # Create MBS table for JSON export
  mbs_table <- sigrap::sig_mbs_table(mbs_counts = mbs_counts)

  #---- Indels ----#
  # plots
  indel_counts <- sigrap::sig_count_indel(vcf_gr = gr, ref_genome = ref_genome)
  p_indel <- sigrap::sig_plot_indel(indel_counts = indel_counts)

  # signature contributions
  sigs_indel <-
    MutationalPatterns::get_known_signatures(muttype = "indel") |>
    {
      \(sigs) sigrap::sig_contribution(mut_mat = indel_counts, signatures = sigs)
    }() |>
    sigrap::sig_contribution_table(type = "ID")

  cli::cli_h2(glue::glue("{date_log()} Saving MutationalPatterns results to\n'{outdir}'"))
  save_plot_list(p_snv, file.path(outdir, "plot/snv"))
  if (strand_bias) {
    save_plot_list(p_strand, file.path(outdir, "plot/strand"))
  }
  save_plot_list(p_dbs, file.path(outdir, "plot/dbs"))
  save_plot_list(p_mbs, file.path(outdir, "plot/mbs"))
  save_plot_list(p_indel, file.path(outdir, "plot/indel"))
  write_jsongz(x = sigs_snv_2015, path = file.path(outdir, "sigs/snv2015.json.gz"))
  write_jsongz(x = sigs_snv_2020, path = file.path(outdir, "sigs/snv2020.json.gz"))
  write_jsongz(x = sigs_dbs, path = file.path(outdir, "sigs/dbs.json.gz"))
  write_jsongz(x = mbs_table, path = file.path(outdir, "sigs/mbs.json.gz"))
  write_jsongz(x = sigs_indel, path = file.path(outdir, "sigs/indel.json.gz"))
  cli::cli_h2(glue::glue("{date_log()} End of MutationalPatterns workflow"))
}

#' Plot Strand Bias Analysis
#'
#' Plots transcriptional and replicative strand bias analysis.
#'
#' @param vcf_gr GRanges containing all mutation types from a single sample.
#' @param ref_genome The BSGenome reference genome to use.
#'
#' @return A list with four ggplot2 objects:
#' - p_transcriptional_bias: transcriptional strand bias plot
#' - p_transcriptional_effect: transcriptional strand bias effect size
#' - p_replicative_bias: replicative strand bias plot  
#' - p_replicative_effect: replicative strand bias effect size
#'
#' @export
sig_plot_strand_bias <- function(vcf_gr, ref_genome) {
  mutpat_gr_snv <- MutationalPatterns::get_mut_type(vcf_gr, type = "snv")
  
  ## ---- Transcriptional ---- ##
  # Only support hg38 for strand bias analysis
  txdb_pkg <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
  if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
    stop(txdb_pkg, " package not available. ",
         "Install with: BiocManager::install('", txdb_pkg, "')")
  }
  
  library(txdb_pkg, character.only = TRUE)
  genes_list <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  # Transcriptional strand bias
  mut_mat_s <- MutationalPatterns::mut_matrix_stranded(
    vcf_list = mutpat_gr_snv, ref_genome = ref_genome, 
    ranges = genes_list, mode = "transcription"
  )
  strand_counts <- MutationalPatterns::strand_occurrences(mut_mat_s, by = "all")
  strand_bias <- MutationalPatterns::strand_bias_test(strand_counts)
  
  p_transcriptional_bias <- MutationalPatterns::plot_strand(strand_counts, mode = "relative") +
    ggplot2::ggtitle("Transcriptional Strand Bias")
  p_transcriptional_effect <- MutationalPatterns::plot_strand_bias(strand_bias) +
    ggplot2::ggtitle("Transcriptional Bias Effect Size")
  
  ## ---- Replicative ---- ##
  repli_file <- system.file("extdata/ReplicationDirectionRegions.bed", package = "MutationalPatterns")
  repli_strand <- readr::read_tsv(repli_file, col_names = TRUE, col_types = "cddcc") |>
    dplyr::mutate_if(is.character, as.factor)
  
  repli_strand_granges <- GenomicRanges::GRanges(
    seqnames = repli_strand$Chr,
    ranges = IRanges::IRanges(start = repli_strand$Start + 1, end = repli_strand$Stop),
    strand_info = repli_strand$Class
  )
  GenomeInfoDb::seqlevelsStyle(repli_strand_granges) <- GenomeInfoDb::seqlevelsStyle(ref_genome)
  
  mut_mat_s_rep <- MutationalPatterns::mut_matrix_stranded(
    vcf_list = mutpat_gr_snv, ref_genome = ref_genome, 
    ranges = repli_strand_granges, mode = "replication"
  )
  strand_counts_rep <- MutationalPatterns::strand_occurrences(mut_mat_s_rep, by = "all")
  strand_bias_rep <- MutationalPatterns::strand_bias_test(strand_counts_rep)
  
  p_replicative_bias <- MutationalPatterns::plot_strand(strand_counts_rep, mode = "relative") +
    ggplot2::ggtitle("Replicative Strand Bias")
  p_replicative_effect <- MutationalPatterns::plot_strand_bias(strand_bias_rep) +
    ggplot2::ggtitle("Replicative Bias Effect Size")
  
  list(
    p_transcriptional_bias = p_transcriptional_bias,
    p_transcriptional_effect = p_transcriptional_effect,
    p_replicative_bias = p_replicative_bias,
    p_replicative_effect = p_replicative_effect
  )
}

#' Plot Rainfall Plot
#'
#' Generates a rainfall plot showing intermutational distances.
#'
#' @param vcf_gr GRanges containing all mutation types from a single sample.
#' @param ref_genome The BSGenome reference genome to use.
#'
#' @return A ggplot2 object representing the rainfall plot, or an informative 
#' placeholder plot if the rainfall plot cannot be generated.
#'
#' @export
sig_plot_rainfall <- function(vcf_gr, ref_genome) {
  gr_snv_rainfall <- MutationalPatterns::get_mut_type(vcf_gr, type = "snv")
  
  # When there is only 1 or lower number of variants on a chromosome,
  # MutationalPatterns::plot_rainfall will crash with an error. So need to check if it will work beforehand.
  chromosomes <- GenomeInfoDb::seqnames(ref_genome)[1:22]
  vcf <- gr_snv_rainfall[[1]]
  will_work <- FALSE
  for (i in 1:length(chromosomes)) {
    chr_subset <- vcf[GenomeInfoDb::seqnames(vcf) == chromosomes[i]]
    n <- length(chr_subset)
    if (n >= 2) {
      will_work <- TRUE
      break
    }
  }
  if (will_work) {
    p_rainfall <- MutationalPatterns::plot_rainfall(vcf,
      chromosomes = GenomeInfoDb::seqnames(ref_genome)[1:22],
      cex = 1.2, ylim = 1e+09
    )
  } else {
    # Create a placeholder when insufficient data
    p_rainfall <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, 
                       label = "Insufficient variants for rainfall plot\n(need ≥2 variants per chromosome)",
                       hjust = 0.5, vjust = 0.5, size = 4) +
      ggplot2::theme_void() +
      ggplot2::labs(title = paste("Rainfall Plot -", names(gr_snv_rainfall)[1]))
  }
  
  return(p_rainfall)
}
