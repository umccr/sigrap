mutpat_add_args <- function(subp) {
  mutpat <- subp$add_parser("mutpat", help = "MutationalPatterns help")
  mutpat$add_argument("--sample", help = "Sample name.", required = TRUE)
  mutpat$add_argument("--snv", help = "Input SNV file (VCF format).", required = TRUE)
  mutpat$add_argument("--outdir", help = "Output directory to write results to.", required = TRUE)
  mutpat$add_argument("--rainfall", help = "Include rainfall plot.", action = "store_true")
  mutpat$add_argument("--strand-bias", help = "Include strand bias analysis.", action = "store_true")
}

mutpat_parse_args <- function(args) {
  cli::cli_h1("Started running MutationalPatterns!")
  
  # Handle potential missing arguments with defaults
  rainfall <- if(is.null(args$rainfall)) FALSE else args$rainfall
  strand_bias <- if(is.null(args$strand_bias)) FALSE else args$strand_bias
  
  sigrap::sig_workflow_run(
    vcf = args$snv, 
    sample_nm = args$sample, 
    outdir = args$out,
    rainfall = rainfall,
    strand_bias = strand_bias
  )
  cli::cli_h1("Finished running MutationalPatterns!")
}
