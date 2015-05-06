# coding: utf-8
"""Command line functions.

"""
import argparse

import poliastro
from poliastro import ephem


def main():
    parser = argparse.ArgumentParser(
        prog="poliastro",
        description="Command line tools for the poliastro Python library.")
    parser.add_argument("--version", action='version',
                        version=poliastro.__version__)

    subparsers = parser.add_subparsers(help='Available commands')

    download_parser = subparsers.add_parser(
        "download-spk",
        help="Downloads .bsp SPICE kernel files")
    download_parser.add_argument(
        "-n", "--name",
        dest="name", default="de430",
        help="Name of a .bsp SPICE kernel file to download")
    download_parser.set_defaults(func=download_kernel)

    args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError:
        parser.print_help()


def download_kernel(args):
    """Downloads SPICE kernel file.

    """
    ephem.download_kernel(args.name)
