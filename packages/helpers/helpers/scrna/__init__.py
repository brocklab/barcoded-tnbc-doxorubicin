import argparse
import logging
from pathlib import Path

logging.getLogger("scrna").setLevel(logging.DEBUG)


def parse_input_output(input_help, output_help) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help=input_help, type=Path, required=True)
    parser.add_argument(
        "--output",
        help=output_help,
        type=Path,
        required=True,
    )

    return parser.parse_args()


def dataset_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", help="name of dataset")
    return parser
