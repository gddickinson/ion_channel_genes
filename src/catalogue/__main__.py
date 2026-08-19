"""`python -m src.catalogue` — validate the catalogue and print its headline counts.

Run at the top of every S0/S1 script and as the fastest possible check that
an edit to a division file did not break the catalogue. Exits non-zero on any
validation problem, so it works as a pre-commit gate.
"""

from __future__ import annotations

import sys

from .registry import (CATALOGUE, SUPERFAMILIES, census_families,
                       shared_signatures, stats, validate)


def main(argv: list[str]) -> int:
    problems = validate()
    st = stats()
    print("ion-channel catalogue")
    for k, v in st.items():
        print(f"  {k:26s} {v}")

    if "--families" in argv:
        print("\nfamilies by superfamily:")
        for sf_key, sf in SUPERFAMILIES.items():
            members = [f for f in CATALOGUE.values() if f.superfamily == sf_key]
            flag = "" if sf.alignable else "   [not alignable — fold network only]"
            print(f"  {sf_key:14s} {sf.name}{flag}")
            for f in members:
                n = len(f.human_genes)
                print(f"      {f.key:24s} {f.status.value:20s} "
                      f"{n:3d} human gene(s)  {f.name}")

    if "--shared" in argv:
        print("\nsignatures carried by more than one family:")
        for acc, fams in sorted(shared_signatures().items(),
                                key=lambda kv: -len(kv[1])):
            print(f"  {acc}  {len(fams):2d}  {', '.join(fams)}")

    if "--census" in argv:
        print("\ncensus families:")
        for f in census_families():
            print(f"  {f.key:24s} {len(f.human_genes):3d}  {f.name}")

    if problems:
        print(f"\n{len(problems)} validation problem(s):")
        for p in problems:
            print(f"  ! {p}")
        return 1
    print("\nvalidation: clean")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
