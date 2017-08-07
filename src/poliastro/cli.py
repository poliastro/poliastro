# coding: utf-8
"""Command line functions.

"""
import argparse

import poliastro
from poliastro.neos.neos_dastcom5 import download_DASTCOM5

def main():
    parser = argparse.ArgumentParser(
        prog="poliastro",
        description="Command line tools for the poliastro Python library.")
    parser.add_argument("--version", action='version',
                        version=poliastro.__version__)

    subparsers = parser.add_subparsers(help='Available commands')

    download_parser = subparsers.add_parser(
        "download-dastcom5",
        help="Downloads DASTCOM5 database")
    download_parser.add_argument(
        "-p", "--path", help="Download path")

    args = parser.parse_args()
    try:
        if args.path:
            download_DASTCOM5(args.path)
        else:
            download_DASTCOM5()
    except AttributeError:
        parser.print_help()

