hrdetect_add_args <- function(subp) {
  hrdetect <- subp$add_parser("hrdetect", help = "HRDetect help")
  hrdetect$add_argument("--sample", help = "Sample name.", required = TRUE)
  hrdetect$add_argument("--snv", help = "Input SNV (VCF format).", required = TRUE)
  hrdetect$add_argument("--sv", help = "Input SV (VCF format).", required = TRUE)
  hrdetect$add_argument("--cnv", help = "Input CNV (TSV format).", required = TRUE)
  hrdetect$add_argument("--out", help = "Output file ['hrdetect.json.gz'].", default = "hrdetect.json.gz")
}


hrdetect_parse_args <- function(args) {
  cli::cli_h1("Started running HRDetect!")
  res <- sigrap::hrdetect_run(
    nm = args$sample, snvindel_vcf = args$snv, sv_vcf = args$sv,
    cnv_tsv = args$cnv, outpath = args$out
  )
  print(t(res))
  cli::cli_h1("Finished running HRDetect!")
}
