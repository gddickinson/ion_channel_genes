"""Novel-paralog discovery pipeline.

Glues together the analysis module (sequence outlier detection) with
structure / domain / orthology cross-checks to score and rank candidates
for being an *unclassified member* of the scoped channel family — an
un-named gene model that carries the family's evidence but no symbol, or a
lineage-specific paralogue the databases have not recognised.

Contents:
    candidates.py — composite scoring + ranked report
"""

from .candidates import (
    DiscoveryConfig,
    DiscoveryReport,
    Candidate,
    build_signature_set,
    discover_novel_paralogs,
    write_discovery,
)
from .domain_scan import run_domain_scan
from .exhaustive import run_exhaustive_hunt

__all__ = [
    "DiscoveryConfig", "DiscoveryReport", "Candidate", "build_signature_set",
    "discover_novel_paralogs", "write_discovery", "run_domain_scan",
    "run_exhaustive_hunt",
]
