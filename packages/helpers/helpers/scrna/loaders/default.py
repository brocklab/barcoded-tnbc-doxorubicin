import helpers as hp

project = hp.Project("scrna/downstream", "helpers.scrna.loader.default")


def main():
    args = hp.scrna.dataset_parser().parse_args()
    dataset = hp.data.get_sc_dataset(args.dataset)

    project.log.info(f"reading raw data for dataset {dataset.name}")
    project.log.info(dataset)

    adata = dataset.load_raw()

    project.log.info(f"writing anndata to {dataset.raw}")
    adata.write(dataset.raw)


if __name__ == "__main__":
    main()
