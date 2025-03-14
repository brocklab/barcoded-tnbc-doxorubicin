#!/usr/bin/env bash

mkdir -p data/references
curl -o data/references/refdata-gex-GRCh38-2024-A.tar.gz "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xvzf data/references/refdata-gex-GRCh38-2024-A.tar.gz -C data/references
