# Barcoded MDA-MB-231 fastq2h5ad

## Setup

Previously GRCh38 reference genome was downloaded from 10X:

```sh
mkdir -p data/references
curl -o data/references/refdata-gex-GRCh38-2024-A.tar.gz "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xvzf data/references/refdata-gex-GRCh38-2024-A.tar.gz -C data/references/refdata-gex-GRCh38-2024-A
```

Downloaded `cellranger-8.0.1` from [10X genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads)
and extracted to `./cellranger-8.0.1/`.

Symlinked all raw sequencing files to `data/raw/` from in-house data store with `./scripts/link-raw.py`.

## Usage
