#!/usr/bin/env python3
"""Stub: BKT / logarithmic-FSS fits are Milestone 3 scope.

Not implemented until the ordinary-FSS core passes its synthetic
validation tests (RESEARCH_PLAN.md, milestone constraint).
"""

import _path  # noqa: F401

import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file", nargs="?")
    args = ap.parse_args()
    raise SystemExit("bkt_fit.py is a Milestone 3 stub; not yet implemented.")


if __name__ == "__main__":
    main()
