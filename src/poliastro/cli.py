"""Command line functions.

"""
import argparse

import poliastro


def main():
    parser = argparse.ArgumentParser(
        prog="poliastro",
        description="Command line tools for the poliastro Python library.",
    )
    parser.add_argument("--version", action="version", version=poliastro.__version__)

    args = parser.parse_args()

    parser.print_help()
