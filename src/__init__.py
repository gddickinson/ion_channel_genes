"""Ion Channel Census — identify, classify and place every ion channel.

Read `INTERFACE.md` before opening any module here.

Subpackages:
    catalogue     — the subject: 90 families, 25 superfamilies, 16 hazards
    classify      — three tiers → one ChannelCall with an audit trail
    phylo         — pore modules, tier-1/tier-2 forests, the tier-3 network
    core          — data models, search orchestrator, on-disk cache
    databases     — one client per source behind a common base
    analysis      — MSA, distance, tree, clusters, mutations, conservation
    discovery     — the candidate scorer and the domain-bait / exhaustive hunts
    investigation — seven evidence lines on a single accession
    gui           — tkinter/ttk widgets
    utils         — scope, MAFFT, species, data root, exporters, reports

The catalogue, classify and phylo subpackages need only the standard library
plus `requests`; the rest need Biopython and matplotlib (Decisions D18).
"""

__version__ = "2.4.0"
