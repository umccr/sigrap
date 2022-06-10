
-   <a href="#id_-sigrap" id="toc-id_-sigrap">🎶 sigrap</a>
    -   <a href="#installation" id="toc-installation">Installation</a>
    -   <a href="#main-modules" id="toc-main-modules">Main Modules</a>
        -   <a href="#id_-hrdetect" id="toc-id_-hrdetect">🔍 HRDetect</a>
        -   <a href="#id_-chord" id="toc-id_-chord">🎸 CHORD</a>
        -   <a href="#id_-mutationalpatterns" id="toc-id_-mutationalpatterns">🐾
            MutationalPatterns</a>
    -   <a href="#id_-cli" id="toc-id_-cli">💻 CLI</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 🎶 sigrap

Wrappers for somatic mutation signature analysis tools (HRDetect, CHORD,
MutationalPatterns).

-   Docs: <https://umccr.github.io/sigrap/>

<!-- badges: start -->

[![Conda
install](https://anaconda.org/umccr/r-sigrap/badges/installer/conda.svg)](https://anaconda.org/umccr/r-sigrap)
<!-- badges: end -->

## Installation

``` r
remotes::install_github("umccr/sigrap")
```

-   Or if used inside a conda environment:

``` bash
conda install r-sigrap -c umccr -c conda-forge -c bioconda
```

## Main Modules

### 🔍 HRDetect

Wraps functionality from the
[HRDetect](https://github.com/Nik-Zainal-Group/signature.tools.lib)
framework - see vignette at
<https://umccr.github.io/sigrap/articles/hrdetect.html>.

### 🎸 CHORD

Wraps functionality from
[CHORD](https://github.com/UMCUGenetics/CHORD) - see vignette at
<https://umccr.github.io/sigrap/articles/chord.html>.

### 🐾 MutationalPatterns

Wraps functionality from
[MutationalPatterns](https://github.com/UMCUGenetics/MutationalPatterns) -
see vignette at
<https://umccr.github.io/sigrap/articles/mutationalpatterns.html>.

## 💻 CLI

A `sigrap` command line interface is available for convenience.

-   If you’re using the conda package, the `sigrap.R` command will
    already be set up inside an activated conda environment.
-   If you’re *not* using the conda package, you need to export the
    `sigrap/inst/cli/` directory to your `PATH` in order to use
    `sigrap.R`.

``` bash
sigrap_cli=$(Rscript -e 'x = system.file("cli", package = "sigrap"); cat(x, "\n")' | xargs)
export PATH="${sigrap_cli}:${PATH}"
```

    $ sigrap.R --version
    sigrap.R 0.0.4

    $ sigrap.R --help
    usage: sigrap [-h] [-v] {hrdetect,chord,mutpat} ...

    Somatic signature wrappers

    positional arguments:
      {hrdetect,chord,mutpat}
                            sub-command help
        hrdetect            HRDetect help
        chord               CHORD help
        mutpat              MutationalPatterns help

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit



    #------- HRDetect -------#


    $ sigrap.R hrdetect --help
    usage: sigrap hrdetect [-h] --sample SAMPLE --snv SNV --sv SV --cnv CNV
                           [--out OUT]

    optional arguments:
      -h, --help       show this help message and exit
      --sample SAMPLE  Sample name.
      --snv SNV        Input SNV (VCF format).
      --sv SV          Input SV (VCF format).
      --cnv CNV        Input CNV (TSV format).
      --out OUT        Output file ['hrdetect.json.gz'].



    #------- CHORD -------#


    $ sigrap.R chord --help
    usage: sigrap chord [-h] --sample SAMPLE --snv SNV --sv SV [--out OUT]

    optional arguments:
      -h, --help       show this help message and exit
      --sample SAMPLE  Sample name.
      --snv SNV        Input SNV (VCF format).
      --sv SV          Input SV (VCF format).
      --out OUT        Output file ['./chord.json.gz']



    #------- MutationalPatterns -------#


    $ sigrap.R mutpat --help
    usage: sigrap mutpat [-h] --sample SAMPLE --snv SNV --outdir OUTDIR

    optional arguments:
      -h, --help       show this help message and exit
      --sample SAMPLE  Sample name.
      --snv SNV        Input SNV file (VCF format).
      --outdir OUTDIR  Output directory to write results to.
