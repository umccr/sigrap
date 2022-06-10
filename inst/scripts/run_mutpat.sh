#!/usr/bin/env bash

timestamp=$(date +%s)
d=$(Rscript -e 'x = system.file("extdata", package="gpgr"); cat(x, "\n")' | xargs)
o="../../nogit/test"

../cli/sigrap.R mutpat \
    --sample "sampleA" \
    --snv "${d}/umccrise/snv/somatic-ensemble-PASS.vcf.gz" \
    --outdir "${o}/mutpat_results_${timestamp}"
