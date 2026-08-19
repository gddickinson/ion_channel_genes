"""Analysis-menu controller.

Gives the GUI access to everything the headless CLI can do: the sequence
analysis pipeline, the discovery scorer (with optional InterPro lookup),
and analysis/discovery/report artefacts inside saved results bundles. The
one-off investigations — selection test, fold check, presence matrix,
domain scan, exhaustive hunt, deep dive, and **classification against the
ion-channel catalogue** — live in `tools_extra.py:ToolsExtraActions`, which
this controller inherits.

All heavy work runs through `TaskRunner` on a worker thread; the GUI never
blocks. Completed results are cached on the controller so "Save results
bundle" persists them alongside the raw hits, exactly like the CLI's
`--analyze --discover --save-results` path.
"""

from __future__ import annotations

import json
import tkinter as tk
import webbrowser
from pathlib import Path
from tkinter import filedialog, messagebox
from typing import Optional

from ..analysis import analyse
from ..analysis.pipeline import AnalysisResult, _label_for, write_analysis
from ..databases.interpro import batch_fetch_family_signatures
from ..discovery import (
    DiscoveryConfig,
    DiscoveryReport,
    build_signature_set,
    discover_novel_paralogs,
    write_discovery,
)
from ..utils.report import write_report
from ..utils.results_writer import make_bundle_dir
from .dialogs import AnalysisDialog, DiscoveryDialog
from .task_runner import TaskRunner
from .tools_extra import ToolsExtraActions
from .text_window import TextWindow


class ToolsController(ToolsExtraActions):
    """Owns the Analysis menu, its background tasks, and cached results."""

    def __init__(self, app) -> None:  # app: MainWindow (avoids circular import)
        self.app = app
        self.runner = TaskRunner(app.root)
        self.analysis: Optional[AnalysisResult] = None
        self.discovery: Optional[DiscoveryReport] = None
        self.presence = None                      # last PresenceReport
        self.novel_labels: set[str] = set()       # analysis labels to highlight

    def attach_menu(self, menubar: tk.Menu) -> None:
        m = tk.Menu(menubar, tearoff=False)
        m.add_command(label="Classify selected row (catalogue)…",
                      command=self.classify_selection)
        m.add_separator()
        m.add_command(label="Run sequence analysis…", command=self.run_analysis)
        m.add_command(label="Discover unclassified members…",
                      command=self.run_discovery)
        m.add_command(label="Domain-bait scan (hunt novel members)…", command=self.run_domain_scan)
        m.add_command(label="Exhaustive hunt (all sources + all analyses)…",
                      command=self.run_exhaustive)
        m.add_separator()
        m.add_command(label="Alignment & tree viewer", command=self.show_phylo_viewer)
        m.add_separator()
        m.add_command(label="Investigate accession (deep dive)…", command=self.investigate)
        m.add_command(label="Selection test (dN/dS)…", command=self.run_selection_test)
        m.add_command(label="Fold check (ESMFold)…", command=self.run_fold_check)
        m.add_command(label="Presence/absence matrix", command=self.show_presence_matrix)
        m.add_separator()
        m.add_command(label="View last analysis summary", command=self.show_analysis)
        m.add_command(label="View last discovery report", command=self.show_discovery)
        menubar.add_cascade(label="Analysis", menu=m)

    def reset(self) -> None:
        """A new search invalidates results computed for the old table."""
        self.analysis = None
        self.discovery = None
        self.presence = None
        self.novel_labels = set()

    # ---- shared helpers --------------------------------------------------

    def _guard_busy(self) -> bool:
        busy = self.runner.busy_with()
        if busy:
            messagebox.showinfo("Busy", f"A background task is already running: {busy}.")
            return True
        return False

    def _task_failed(self, tb: str) -> None:
        self.app.status("Background task failed — see details window.")
        TextWindow(self.app.root, "Task failed", tb, save_name="error.txt",
                   width=760, height=420)

    # ---- sequence analysis ----------------------------------------------

    def run_analysis(self) -> None:
        if self._guard_busy():
            return
        variants = self.app.results.all_variants()
        n_seq = sum(1 for v in variants if v.sequence)
        if not variants:
            messagebox.showinfo("Analysis", "No results to analyse — run a search first.")
            return
        if n_seq < 2:
            messagebox.showinfo(
                "Analysis",
                "Need at least 2 variants with sequences.\n"
                "Re-search with 'Fetch sequences' checked.")
            return
        dlg = AnalysisDialog(self.app.root, n_with_seq=n_seq)
        if not dlg.result:
            return
        opts = dlg.result

        def job(progress):
            progress(f"Aligning {min(n_seq, opts['max_variants'])} sequence(s) — this can take a while…")
            return analyse(
                variants,
                identity_threshold=opts["identity_threshold"],
                conservation_threshold=opts["conservation_threshold"],
                use_mafft=opts["use_mafft"],
                max_variants=opts["max_variants"],
            )

        if self.runner.run("sequence analysis", job, self._analysis_done,
                           self._task_failed, self.app.status):
            self.app.status(f"Running sequence analysis on {n_seq} sequence(s)…")

    def _analysis_done(self, result: AnalysisResult) -> None:
        self.analysis = result
        self.app.status(
            f"Analysis done: {result.n_analyzed}/{result.n_input} variants analysed. "
            "Included in the next saved bundle.")
        self.show_analysis()

    def show_analysis(self) -> None:
        if not self.analysis:
            messagebox.showinfo("Analysis", "No analysis yet — run 'Run sequence analysis…' first.")
            return
        TextWindow(self.app.root, "Sequence analysis", self.analysis.text_summary(),
                   save_name="analysis_summary.txt",
                   extra_actions=[("Save artefacts…", self._save_analysis_artefacts)])

    def _save_analysis_artefacts(self) -> None:
        target = filedialog.askdirectory(title="Folder for the analysis/ artefacts")
        if not target or not self.analysis:
            return
        written = write_analysis(Path(target) / "analysis", self.analysis)
        messagebox.showinfo("Analysis", f"Wrote {len(written)} file(s) to {Path(target) / 'analysis'}")

    # ---- unclassified-member discovery -----------------------------------

    def run_discovery(self) -> None:
        if self._guard_busy():
            return
        variants = self.app.results.all_variants()
        if not variants:
            messagebox.showinfo("Discovery", "No results to score — run a search first.")
            return
        # Honour the loaded preset's discovery tuning (same keys the
        # headless CLI reads), so a preset-driven method is repeatable
        # from the GUI.
        preset = getattr(self.app.search, "loaded_preset", {}) or {}
        default_known = list(preset.get("known_paralogs", [])) \
            or (list(self.app.last_query.gene_symbols) if self.app.last_query else [])
        dlg = DiscoveryDialog(self.app.root, default_known=default_known,
                              has_analysis=self.analysis is not None)
        if not dlg.result:
            return
        known = dlg.result["known_paralogs"]
        do_interpro = dlg.result["interpro"]
        config_kwargs = dict(preset.get("discovery", {}))
        analysis = self.analysis

        def job(progress):
            domain_hits: dict = {}
            if do_interpro:
                accs = sorted({v.accession for v in variants if v.source == "UniProt"})
                progress(f"InterPro: fetching Pfam signatures for {len(accs)} UniProt entrie(s)…")
                domain_hits = batch_fetch_family_signatures(accs)
            label_map = {f"{v.source}|{v.accession}|{v.gene_symbol}": _label_for(v)
                         for v in variants}
            sset = None
            if analysis and analysis.msa:
                sset = build_signature_set(variants, analysis.msa, label_map, known)
                if sset:
                    progress(f"Derived {len(sset.signatures)} family-signature blocks "
                             "from the known-paralog MSA…")
            progress("Scoring candidates…")
            return discover_novel_paralogs(
                variants,
                DiscoveryConfig(known_paralogs=known, **config_kwargs),
                analysis_label_for=label_map,
                distances=(analysis.distances if analysis else None),
                domain_hits=domain_hits,
                foldseek_hits=[v for v in variants if v.source == "Foldseek"],
                signature_set=sset,
            )

        if self.runner.run("discovery", job, self._discovery_done,
                           self._task_failed, self.app.status):
            self.app.status("Running discovery scorer…")

    def _discovery_done(self, report: DiscoveryReport) -> None:
        self.discovery = report
        self.novel_labels = self._novel_labels_from_discovery()
        bins = report.by_verdict()
        self.app.status(
            f"Discovery done: {len(bins['promising'])} promising, "
            f"{len(bins['worth manual review'])} worth review, {len(bins['weak'])} weak.")
        self.show_discovery()

    def show_discovery(self) -> None:
        if not self.discovery:
            messagebox.showinfo("Discovery", "No discovery report yet — run 'Discover unclassified members…' first.")
            return
        TextWindow(self.app.root, "Novel-paralog discovery", self.discovery.text_summary(),
                   save_name="discovery_report.txt",
                   extra_actions=[("Save TSV + report…", self._save_discovery_artefacts)])

    def _save_discovery_artefacts(self) -> None:
        target = filedialog.askdirectory(title="Folder for the discovery/ artefacts")
        if not target or not self.discovery:
            return
        paths = write_discovery(Path(target) / "discovery", self.discovery)
        messagebox.showinfo("Discovery", "Wrote:\n" + "\n".join(paths.values()))

    # ---- selection test (dN/dS) ------------------------------------------

    def augment_bundle(self, bundle_dir: Path, variants, query, notes: str = "") -> dict:
        """Write analysis/, discovery/, and report.md|.html into a bundle,
        mirroring the CLI's --analyze/--discover/--save-results path.
        Returns the report paths."""
        discovery_text = ""
        if self.analysis:
            write_analysis(bundle_dir / "analysis", self.analysis)
        if self.discovery:
            discovery_text = self.discovery.text_summary()
            write_discovery(bundle_dir / "discovery", self.discovery)
        return write_report(bundle_dir, variants, query,
                            analysis=self.analysis,
                            discovery_summary=discovery_text,
                            notes=notes)
