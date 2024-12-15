chord_add_args <- function(subp) {
  chord <- subp$add_parser("chord", help = "CHORD help")
  chord$add_argument("--sample", help = "Sample name.", required = TRUE)
  chord$add_argument("--snv", help = "Input SNV (VCF format).", required = TRUE)
  chord$add_argument("--sv", help = "Input SV (VCF format).", required = TRUE)
  chord$add_argument("--out", help = "Output file ['./chord.json.gz']", default = "./chord.json.gz")
}

chord_parse_args <- function(args) {
  cli::cli_h1("Started running CHORD!")
  # CHORD requires this to be loaded...
  suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
  res <- sigrap::chord_run(
    vcf.snv = args$snv, vcf.sv = args$sv,
    sample.name = args$sample, outpath = args$out
  )
  print(t(res[["prediction"]]))
  cli::cli_h1("Finished running CHORD!")
}
