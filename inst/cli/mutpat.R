mutpat_add_args <- function(subp) {
  mutpat <- subp$add_parser("mutpat", help = "MutationalPatterns help")
  mutpat$add_argument("--sample", help = "Sample name.", required = TRUE)
  mutpat$add_argument("--snv", help = "Input SNV file (VCF format).", required = TRUE)
  mutpat$add_argument("--outdir", help = "Output directory to write results to.", required = TRUE)
}

mutpat_parse_args <- function(args) {
  cli::cli_h1("Started running MutationalPatterns!")
  sigrap::sig_workflow_run(vcf = args$snv, sample_nm = args$sample, outdir = args$outdir)
  cli::cli_h1("Finished running MutationalPatterns!")
}
