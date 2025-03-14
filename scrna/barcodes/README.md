# DM2103 Single-Cell Lineages (Internal Use Only)

## Setup

Initialize the environment using `mamba`:

```bash
mamba env create -p ./env -f environment.yml
mamba activate ./env
```

## Pipeline

Update the `./config.yml` with the requisite metadata, then run `snakemake all`

## Notebooks

There are a number of notebooks that have been generated to perform EDA.

These need to be streamlined and their output saved.

## DAG

![DAG](./dag.svg)