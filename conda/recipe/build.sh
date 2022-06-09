#!/bin/bash

export DISABLE_AUTOBREW=1
${R} CMD INSTALL --build . ${R_ARGS}

# Copy CLI to conda bin
mkdir -p ${PREFIX}/bin
cp ${SRC_DIR}/inst/cli/sigrap.R ${PREFIX}/bin/sigrap.R
chmod +x ${PREFIX}/bin/sigrap.R
