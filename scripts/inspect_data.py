#!/usr/bin/env python3
"""Inspect a canonical FSS data file: schema, sizes, control range, observables."""

import _path  # noqa: F401

import argparse
import json

import pandas as pd

import fss


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file", help="data file (csv/tsv)")
    ap.add_argument("--control", help="control column override")
    ap.add_argument("--size", help="size column override")
    ap.add_argument("--volume", help="volume column override")
    ap.add_argument("--observable", action="append", default=[], metavar="NAME=COL",
                    help="observable mapping (repeatable)")
    ap.add_argument("--dimension", type=float, help="spatial dimension d")
    ap.add_argument("--geometry", help="geometry (hypercubic, periodic, ...)")
    ap.add_argument("--json", action="store_true", help="emit info as JSON")
    args = ap.parse_args()

    observables = dict(o.split("=", 1) for o in args.observable) if args.observable else None
    data = fss.FSSData.load(
        args.file,
        control=args.control, size=args.size, volume=args.volume,
        observables=observables,
        dimension=args.dimension, geometry=args.geometry,
    )
    info = data.info()
    if args.json:
        print(json.dumps(info, indent=2))
        return
    print(f"rows:        {info['n_rows']}")
    print(f"size col:    {info['size_col']}   volume col: {info['volume_col']}")
    print(f"  sizes ({info['n_sizes']}): {info['sizes']}")
    print(f"  range:     {info['size_min']} .. {info['size_max']}")
    print(f"control col: {info['control_col']}   range {info['control_range']}"
          f"   ({info['n_control_values']} values)")
    print(f"observables: {info['observables']}")
    print(f"errors:      {info['errors']}")
    print(f"with errors: {info['with_errors']}   missing: {info['has_missing']}")
    print(f"dimension:   {info['dimension']}   geometry: {info['geometry']}")
    print(f"volume_from_size: {info['volume_from_size']}")
    print("\nfirst 5 rows:")
    print(data.df.head().to_string(index=False))


if __name__ == "__main__":
    main()
