"""Command line functions.

"""
import argparse

import poliastro
from poliastro.neos.dastcom5 import download_dastcom5


def main():
    parser = argparse.ArgumentParser(
        prog="poliastro",
        description="Command line tools for the poliastro Python library.")
    parser.add_argument("--version", action='version',
                        version=poliastro.__version__)
    parser.add_argument("--download-dastcom5", action="store_true",
                        help="Downloads DASTCOM5 database")

    args = parser.parse_args()

    if args.download_dastcom5:
        download_dastcom5()
    else:
        parser.print_help()
