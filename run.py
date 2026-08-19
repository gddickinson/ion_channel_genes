"""Entry point for the ion-channel census, classifier and phylogeny toolkit.

Classification and catalogue commands run in a bare interpreter (stdlib +
`requests`); the search / analysis / GUI paths need Biopython and matplotlib
(see `results/toolchain_manifest.txt` and decision **D18**).

Catalogue and classification:
    python run.py --catalogue                       # the whole catalogue
    python run.py --catalogue --scope ploop         # one superfamily
    python run.py --classify Q14524,O43497          # by UniProt accession
    python run.py --classify KCNQ1,KCTD1 --save-results
    python run.py --classify-preset channelome_human --save-results

Phylogeny (MAFFT → trimAl → IQ-TREE 2, D27-guarded):
    python run.py --phylo nav                       # tier 1: within a family
    python run.py --phylo cysloop --tier 2          # tier 2: pore module
    python run.py --phylo tmem16_like --tier 2      # refused, with the reason

Search (needs Biopython):
    python run.py                                    # GUI
    python run.py --headless --preset channelome_human --save-results
    python run.py --headless --preset ploop_domain_scan --save-results
    python run.py --investigate A0A0P1B5Q5 --scope ploop
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description="Ion-channel census & classifier")
    p.add_argument("--email", default="",
                   help="Email for NCBI Entrez identification.")
    p.add_argument("--preset", default="",
                   help="Preset name (see presets/).")
    p.add_argument("--scope", default="all",
                   help="Catalogue scope: 'all', a superfamily key (ploop, "
                        "cysloop, …) or a family key (nav, kir, itpr, …).")
    p.add_argument("--auto-search", action="store_true",
                   help="GUI mode: run the preset's search on startup.")
    p.add_argument("--headless", action="store_true",
                   help="Run a preset search without a GUI.")
    p.add_argument("--save-results", action="store_true",
                   help="Write a results bundle under results/.")
    p.add_argument("--species", default="",
                   help="Override preset species, or a panel name "
                        "(human, vertebrate, metazoan, eukaryote, "
                        "prokaryote, outgroup, tree_of_life).")
    p.add_argument("--max", type=int, default=0,
                   help="Override preset max_per_source.")
    p.add_argument("--label", default="", help="Results-folder label.")
    p.add_argument("--analyze", action="store_true",
                   help="Headless: run sequence analysis (needs Biopython).")
    p.add_argument("--discover", action="store_true",
                   help="Headless: rank unclassified candidates in scope.")
    p.add_argument("--interpro", action="store_true",
                   help="Discovery: fetch InterPro signatures per accession.")
    p.add_argument("--known-paralog", action="append", default=[],
                   help="Repeatable; genes that count as already known.")
    p.add_argument("--investigate", default="",
                   help="Deep-dive one UniProt accession.")
    p.add_argument("--vs", default="",
                   help="Comma-separated accessions for the comparison panel "
                        "(default: the scope's exemplars plus its sister "
                        "families, so the hazard tests have something to "
                        "measure against).")
    p.add_argument("--foldseek", action="store_true",
                   help="Investigation: also run Foldseek (slow, async).")
    # --- catalogue / classify / phylogeny -------------------------------
    p.add_argument("--catalogue", action="store_true",
                   help="Print the catalogue for --scope and exit.")
    p.add_argument("--detail", action="store_true",
                   help="--catalogue: also print architectures and exemplars.")
    p.add_argument("--classify", default="",
                   help="Comma-separated UniProt accessions or human gene "
                        "symbols to classify into the catalogue.")
    p.add_argument("--classify-preset", default="",
                   help="Classify every accession/symbol in a preset's "
                        "'classify' list.")
    p.add_argument("--no-reference", action="store_true",
                   help="Classification: skip the reference-identity tier.")
    p.add_argument("--phylo", default="",
                   help="Build a tree for a family (tier 1) or superfamily "
                        "(tier 2) key.")
    p.add_argument("--tier", type=int, default=0,
                   help="Force the phylogeny tier (1 or 2).")
    p.add_argument("--bootstrap", type=int, default=1000,
                   help="IQ-TREE ultrafast bootstrap replicates.")
    p.add_argument("--linsi", action="store_true",
                   help="MAFFT L-INS-i instead of --auto (slower, better).")
    args = p.parse_args()

    root = Path(__file__).resolve().parent

    if args.catalogue:
        from src.cli_channel import run_catalogue
        run_catalogue(args.scope, show_families=args.detail)
        return

    if args.classify or args.classify_preset:
        from src.cli_channel import run_classify
        targets: list[str] = []
        expected: dict = {}
        if args.classify:
            targets += [t for t in args.classify.split(",") if t.strip()]
        if args.classify_preset:
            data = json.loads(
                (root / "presets" / f"{args.classify_preset}.json").read_text())
            targets += list(data.get("classify", []))
            expected = dict(data.get("expected", {}))
        run_classify(root, targets, save_results=args.save_results,
                     label=args.label or args.classify_preset or "classify",
                     use_reference=not args.no_reference,
                     expected=expected or None)
        return

    if args.phylo:
        from src.cli_channel import run_phylo
        run_phylo(root, args.phylo, tier=args.tier,
                  bootstrap=args.bootstrap, linsi=args.linsi)
        return

    if args.investigate:
        from src.cli import run_investigate
        run_investigate(project_root=root, accession=args.investigate,
                        email=args.email, panel_csv=args.vs,
                        run_foldseek=args.foldseek,
                        save_results=args.save_results,
                        label_override=args.label, scope_key=args.scope)
        return

    if args.headless:
        if not args.preset:
            p.error("--headless requires --preset")
        from src.cli import run_headless
        run_headless(project_root=root, preset=args.preset, email=args.email,
                     species_override=args.species, max_override=args.max,
                     save_results=args.save_results,
                     label_override=args.label, analyze=args.analyze,
                     discover=args.discover,
                     known_paralogs=args.known_paralog or None,
                     interpro_lookup=args.interpro, scope_key=args.scope)
        return

    from src.gui import run as run_gui
    run_gui(project_root=root, email=args.email, preset=args.preset,
            auto_search=args.auto_search)


if __name__ == "__main__":
    main()
