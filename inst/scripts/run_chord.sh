#!/usr/bin/env bash

timestamp=$(date +%s)
d=$(Rscript -e 'x = system.file("extdata", package="gpgr"); cat(x, "\n")' | xargs)
o="../../nogit/test"

../cli/sigrap.R chord \
    --sample "sampleA" \
    --snv "${d}/umccrise/snv/somatic-ensemble-PASS.vcf.gz" \
    --sv "${d}/umccrise/sv/manta.vcf.gz" \
    --out "${o}/chord_results_${timestamp}.json.gz"
