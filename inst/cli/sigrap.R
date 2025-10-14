#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(cli))
suppressPackageStartupMessages(require(glue))

devtools::load_all("~/github/sigrap")
prog_nm <- "sigrap"
sigrap_version <- as.character(packageVersion("sigrap"))
p <- argparse::ArgumentParser(description = "Somatic signature wrappers", prog = prog_nm)
p$add_argument("-v", "--version", action = "version", version = glue::glue("{prog_nm} {sigrap_version}"))
subparser_name <- "subparser_name"
subp <- p$add_subparsers(help = "sub-command help", dest = subparser_name)

source(system.file("cli/hrdetect.R", package = "sigrap"))
source(system.file("cli/chord.R", package = "sigrap"))
source(system.file("cli/mutpat.R", package = "sigrap"))

hrdetect_add_args(subp)
chord_add_args(subp)
mutpat_add_args(subp)

args <- p$parse_args()
if (length(args$subparser_name) == 0) {
  p$print_help()
} else if (args$subparser_name == "hrdetect") {
  hrdetect_parse_args(args)
} else if (args$subparser_name == "chord") {
  chord_parse_args(args)
} else if (args$subparser_name == "mutpat") {
  mutpat_parse_args(args)
} else {
  stop("NO IDEA HOW IT GOT TO THIS...")
}
