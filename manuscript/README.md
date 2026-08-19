# manuscript/ — the submission package

**Nothing here yet.** The package is built in S14a; this file records how,
so the build is not reinvented.

## How it builds

```
python3 scripts/s14_assemble.py            # figures → claims → stitch → deposit
python3 scripts/s14_assemble.py --list     # the stages
python3 scripts/s14_pdf.py                 # typeset PDF (pandoc + xelatex)
```

`s14_assemble.py` exits non-zero if any figure, section or load-bearing
number is missing — which it will, until the analyses have run. That is the
intended behaviour, not a bug to work around.

## Rules

- **Never edit `manuscript.md`.** It is stitched from the numbered section
  files listed in `scripts/s14_lib.py:SECTION_ORDER` and overwritten on
  every build.
- **Figures are copied, never re-plotted.** `s14_figures.py` copies each
  committed figure from its results directory under its publication number,
  so a manuscript figure cannot differ from the one the analysis produced.
- **Every load-bearing number needs a claim row** in `scripts/s14_claims.py`
  naming its source table and the operation that recovers it (D12).
- **Look at the figure** before writing its legend (D11).
- If the framing changes substantially, freeze the old draft in
  `manuscript_v1/` rather than overwriting it.

## What this paper has to get right that a single-family paper does not

1. **State the denominator in the abstract.** Published human channelome
   counts run 240–400 and the spread is almost all scope. The scope
   decisions (`docs/scope_and_boundaries.md`, D23) belong in the first
   figure, not in the supplement.
2. **Never draw one tree.** Every phylogenetic figure is a tier-1 or tier-2
   panel, labelled with its tier and, for tier 2, with the words *pore
   module*. The cross-superfamily panel is a network with no branch lengths
   (D27). A reviewer who reads a tier-2 panel as a protein tree has been
   misled by the figure, not by the reader.
3. **Report which tier made each call.** The per-tier attribution table from
   S1 is a main-text result: it is the difference between a classifier and a
   nearest-neighbour lookup.
4. **Report the untested hazards.** Sixteen are defined; the ones no panel
   member exercises are limits of the benchmark and are stated as such.

## Planned figure set

Mapped in `scripts/s14_lib.py` to the results paths that will produce them.
The working title is *A declared census, a benchmarked classification and a
forest phylogeny of the ion channels*; the real title is chosen in S14a from
what was actually found.
