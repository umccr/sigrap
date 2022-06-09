
-   <a href="#sigrap" id="toc-sigrap">sigrap</a>
    -   <a href="#installation" id="toc-installation">Installation</a>
    -   <a href="#main-modules" id="toc-main-modules">Main Modules</a>
        -   <a href="#hrdetect" id="toc-hrdetect">HRDetect</a>
        -   <a href="#chord" id="toc-chord">CHORD</a>
        -   <a href="#mutationalpatterns"
            id="toc-mutationalpatterns">MutationalPatterns</a>
    -   <a href="#cli" id="toc-cli">CLI</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigrap

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

### HRDetect

Wraps functionality from the
[HRDetect](https://github.com/Nik-Zainal-Group/signature.tools.lib)
framework - see vignette at
<https://umccr.github.io/sigrap/articles/hrdetect.html>.

### CHORD

Wraps functionality from
[CHORD](https://github.com/UMCUGenetics/CHORD) - see vignette at
<https://umccr.github.io/sigrap/articles/chord.html>.

### MutationalPatterns

Wraps functionality from
[MutationalPatterns](https://github.com/UMCUGenetics/MutationalPatterns) -
see vignette at
<https://umccr.github.io/sigrap/articles/mutationalpatterns.html>.

## CLI

    $ sigrap --help
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



    $ sigrap --version
    sigrap 0.0.4



    #------- HRDetect -------#


    $ sigrap hrdetect --help
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


    $ sigrap chord --help
    usage: sigrap chord [-h] --sample SAMPLE --snv SNV --sv SV [--out OUT]

    optional arguments:
      -h, --help       show this help message and exit
      --sample SAMPLE  Sample name.
      --snv SNV        Input SNV (VCF format).
      --sv SV          Input SV (VCF format).
      --out OUT        Output file ['./chord.json.gz']



    #------- MutationalPatterns -------#


    $ sigrap mutpat --help
    usage: sigrap mutpat [-h] --sample SAMPLE --snv SNV --outdir OUTDIR

    optional arguments:
      -h, --help       show this help message and exit
      --sample SAMPLE  Sample name.
      --snv SNV        Input SNV file (VCF format).
      --outdir OUTDIR  Output directory to write results to.
