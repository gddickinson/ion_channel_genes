# The Ion Channel Census — findings summary

*This is the illustrated summary page's source. `scripts/build_findings_page.py`
renders it to one self-contained HTML file with every linked figure inlined,
plus a live summary card read from the committed census tables (omitted until
those tables exist, so the page can never show a card of zeroes).*

*Rebuild:* `python3 scripts/build_findings_page.py`

---

## What this project is counting

Ion channels are not one family. They are roughly twenty-five independent
inventions of the same trick — a gated hole in a membrane — and the word
covers proteins that share no ancestor, no fold and no alignable position.
A nicotinic acetylcholine receptor and a potassium channel are both ion
channels and cannot be put in the same alignment.

So "how many ion channels are there?" has no answer until somebody says what
counts. Published totals for the human genome run from about 240 to about
400, and almost all of that spread is six decisions, not six disagreements
about biology: do the auxiliary subunits count, do the transporters built on
a channel fold count, do the scramblases, do the large-pore channels, do the
aquaporins, do the viral and bacterial channels. This project answers all
six in writing before counting anything.

Under those answers the catalogue holds **90 families in 25 superfamilies**,
of which 68 are ion channels covering **320 human pore-forming genes**. The
other 22 families are in the catalogue precisely so they can be excluded.

## Findings so far

*(The findings from each completed task are appended here as they land. The
current entries are in `FINDINGS.md`; this page is rebuilt from them at the
end of each session that produces a figure.)*

