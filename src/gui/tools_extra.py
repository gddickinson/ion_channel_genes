"""Analysis-menu actions that are not the core search → analyse → discover path.

Split out of `tools.py` when that file reached the project's 500-line limit.
The seam is deliberate rather than arbitrary: `tools.py` owns the menu, the
task runner and the pipeline the headless CLI also runs, while this mixin
owns the *investigations* — the one-off questions asked of a selected row.

`ToolsController` inherits from `ToolsExtraActions`, so every method here is
a `ToolsController` method and `self` is the controller.
"""

from __future__ import annotations

import json
from pathlib import Path
from tkinter import messagebox

from ..analysis.esmfold import fold_segment
from ..analysis.pipeline import _label_for
from ..analysis.presence import build_presence_matrix
from ..analysis.selection import selection_test
from ..discovery import run_domain_scan, run_exhaustive_hunt
from ..investigation import InvestigationOptions, Investigator, write_investigation
from ..investigation.synthesis import (assess_signals, render_case_file,
                                       verdict_from_signals)
from .dialogs import (DomainScanDialog, FoldCheckDialog, InvestigateDialog,
                      SelectionTestDialog)
from .text_window import TextWindow


class ToolsExtraActions:
    """Selection test, fold check, presence matrix, scans, investigations."""

    def run_selection_test(self) -> None:
        if self._guard_busy():
            return
        sel = self.app.results.selected()
        default = ""
        if sel:
            # Compara rows carry an Ensembl protein id as accession.
            default = sel.accession if sel.accession.startswith("ENS") else \
                sel.raw.get("protein_id", "")
        dlg = SelectionTestDialog(self.app.root, default_accession=default)
        if not dlg.result:
            return
        acc, partner, k = (dlg.result["accession"], dlg.result["partner"],
                           dlg.result["max_partners"])

        def job(progress):
            return selection_test(acc, partner_id=partner, max_partners=k,
                                  log=progress)

        if self.runner.run(f"selection test on {acc}", job,
                           self._selection_done, self._task_failed,
                           lambda m: self.app.status(f"dN/dS {acc}: {m}")):
            self.app.status(f"Running selection test on {acc}…")

    def _selection_done(self, report) -> None:
        self.app.status(
            f"Selection test done: dN/dS "
            f"{'n/a' if report.result.ratio is None else f'{report.result.ratio:.2f}'}"
            f" — {report.result.verdict()}")
        TextWindow(self.app.root, f"Selection test — {report.query_id}",
                   report.summary(), save_name=f"dnds_{report.query_id}.txt",
                   width=780, height=360)

    # ---- fold check (ESMFold) --------------------------------------------

    def run_fold_check(self) -> None:
        if self._guard_busy():
            return
        sel = self.app.results.selected()
        if not sel:
            messagebox.showinfo("Fold check", "Select a row first.")
            return
        if not sel.sequence:
            messagebox.showinfo("Fold check",
                                "The selected row has no sequence — re-search with "
                                "'Fetch sequences' checked.")
            return
        dlg = FoldCheckDialog(self.app.root, seq_len=len(sel.sequence))
        if not dlg.result:
            return
        seq, where, length = sel.sequence, dlg.result["where"], dlg.result["length"]
        label = f"{sel.gene_symbol} ({sel.accession})"

        def job(progress):
            progress("Submitting segment to ESMFold (30–120 s)…")
            return fold_segment(seq, where=where, length=length)

        if self.runner.run(f"ESMFold check on {sel.accession}", job,
                           lambda r: self._fold_done(label, r),
                           self._task_failed, self.app.status):
            self.app.status(f"ESMFold: folding {length} aa of {label}…")

    def _fold_done(self, label: str, result) -> None:
        self.app.status(f"ESMFold done: mean pLDDT {result.mean_plddt:.1f} — "
                        f"{result.verdict()}")

        def save_pdb() -> None:
            path = filedialog.asksaveasfilename(
                defaultextension=".pdb", initialfile="esmfold_segment.pdb",
                filetypes=[("PDB", "*.pdb"), ("All files", "*.*")])
            if path:
                Path(path).write_text(result.pdb_text)

        TextWindow(self.app.root, f"ESMFold — {label}",
                   result.summary(), save_name="esmfold_check.txt",
                   extra_actions=[("Save PDB…", save_pdb)],
                   width=760, height=300)

    # ---- presence/absence matrix -----------------------------------------

    def show_presence_matrix(self) -> None:
        variants = self.app.results.all_variants()
        if not variants:
            messagebox.showinfo("Presence/absence", "No results — run a search first.")
            return
        label_map = {f"{v.source}|{v.accession}|{v.gene_symbol}": _label_for(v)
                     for v in variants}
        report = build_presence_matrix(
            variants,
            distances=(self.analysis.distances if self.analysis else None),
            analysis_label_for=label_map,
        )
        if not self.analysis:
            report.notes.append(
                "Tip: run 'Run sequence analysis…' first so unnamed cross-species "
                "orthologs merge into homolog clusters instead of singletons.")
        TextWindow(self.app.root, "Family presence/absence matrix",
                   report.summary(), save_name="presence_matrix.txt")

    # ---- exhaustive hunt -------------------------------------------------

    def run_special_preset(self, name: str, data: dict) -> None:
        """Dispatch presets carrying a 'mode' key (loaded via the search
        panel's preset picker) to the matching tool."""
        mode = data.get("mode", "")
        if mode == "domain_scan":
            self.run_domain_scan(preset=name)
        elif mode == "exhaustive":
            self.run_exhaustive(preset=name)
        else:
            messagebox.showinfo("Preset", f"Unknown preset mode: {mode!r}")

    def run_exhaustive(self, preset: str = "") -> None:
        if self._guard_busy():
            return
        presets_dir = self.app.project_root / "presets"
        candidates = []
        for p in sorted(presets_dir.glob("*.json")):
            try:
                if json.loads(p.read_text()).get("mode") == "exhaustive":
                    candidates.append(p.stem)
            except (OSError, json.JSONDecodeError):
                continue
        if not candidates:
            messagebox.showinfo("Exhaustive hunt",
                                'No presets with "mode": "exhaustive" found in presets/.')
            return
        name = preset if preset in candidates else candidates[0]
        preset_data = json.loads((presets_dir / f"{name}.json").read_text())
        if not messagebox.askyesno(
                "Exhaustive hunt",
                f"Run the full hunt from preset '{name}'?\n\n"
                "Harvests every sequence in the active scope (Compara mine + ortholog "
                "expansion + InterPro), saves the census, and runs every "
                "analysis. Takes 5–15 minutes against live APIs."):
            return
        project_root, email = self.app.project_root, self.app.email

        def job(progress):
            return run_exhaustive_hunt(
                project_root=project_root, preset_data=preset_data,
                preset_name=name, email=email, save_results=True,
                log=progress,
            )

        if self.runner.run("exhaustive hunt", job, self._exhaustive_done,
                           self._task_failed, self.app.status):
            self.app.status(f"Exhaustive hunt '{name}' running — Compara mine first…")

    def _exhaustive_done(self, out: dict) -> None:
        variants = out.get("variants") or []
        if not variants:
            self.app.status("Exhaustive hunt finished with no sequences.")
            return
        self.analysis = out.get("analysis")
        self.discovery = out.get("discovery_report")
        self.presence = out.get("presence_report")
        self.novel_labels = out.get("novel_labels") or set()
        self.app.load_variants(variants, query=out.get("query"), results=out.get("results"))
        msg = (f"Exhaustive hunt done: {len(variants)} candidate protein(s) in census. "
               f"Bundle: {out.get('dir', '?')}")
        self.app.status(msg)
        self.show_phylo_viewer()
        html = out.get("report_html")
        if html and messagebox.askyesno("Exhaustive hunt", msg + "\n\nOpen the HTML report?"):
            webbrowser.open(Path(html).resolve().as_uri())

    # ---- alignment & tree viewer -----------------------------------------

    def show_phylo_viewer(self) -> None:
        if not self.analysis or not self.analysis.msa:
            messagebox.showinfo(
                "Alignment & tree",
                "No alignment yet — run 'Run sequence analysis…' (or the "
                "exhaustive hunt) first.")
            return
        from .phylo_view import PhyloWindow
        PhyloWindow(self.app.root, analysis=self.analysis,
                    novel_labels=getattr(self, "novel_labels", set()) or
                    self._novel_labels_from_discovery())

    def _novel_labels_from_discovery(self) -> set[str]:
        if not self.discovery:
            return set()
        variants = self.app.results.all_variants()
        label_map = {f"{v.source}|{v.accession}|{v.gene_symbol}": _label_for(v)
                     for v in variants}
        return {label_map.get(c.label, "") for c in self.discovery.candidates
                if c.score >= 40} - {""}

    # ---- domain-bait scan ------------------------------------------------

    def run_domain_scan(self, preset: str = "") -> None:
        if self._guard_busy():
            return
        presets_dir = self.app.project_root / "presets"
        scan_presets = []
        for p in sorted(presets_dir.glob("*.json")):
            try:
                if json.loads(p.read_text()).get("mode") == "domain_scan":
                    scan_presets.append(p.stem)
            except (OSError, json.JSONDecodeError):
                continue
        if not scan_presets:
            messagebox.showinfo("Domain-bait scan",
                                'No presets with "mode": "domain_scan" found in presets/.')
            return
        dlg = DomainScanDialog(self.app.root, scan_presets, default=preset)
        if not dlg.result:
            return
        name = dlg.result["preset"]
        save = dlg.result["save_results"]
        preset_data = json.loads((presets_dir / f"{name}.json").read_text())
        project_root, email = self.app.project_root, self.app.email

        def job(progress):
            return run_domain_scan(
                project_root=project_root, preset_data=preset_data, preset_name=name,
                email=email, save_results=save, analyze=True, discover=True,
                log=progress,
            )

        if self.runner.run("domain-bait scan", job, self._domain_scan_done,
                           self._task_failed, self.app.status):
            self.app.status(f"Domain-bait scan '{name}' running — enumerating InterPro…")

    def _domain_scan_done(self, out: dict) -> None:
        variants = out.get("variants") or []
        if not variants:
            self.app.status("Domain scan finished with no candidates — see preset filters.")
            return
        self.analysis = out.get("analysis")
        self.discovery = out.get("discovery_report")
        self.app.load_variants(variants, query=out.get("query"), results=out.get("results"))
        msg = f"Domain scan done: {len(variants)} candidate(s) loaded into the table."
        if out.get("dir"):
            msg += f"  Bundle: {out['dir']}"
        self.app.status(msg)
        if self.discovery:
            self.show_discovery()
        html = out.get("report_html")
        if html and messagebox.askyesno("Domain-bait scan", msg + "\n\nOpen the HTML report?"):
            webbrowser.open(Path(html).resolve().as_uri())

    # ---- deep-dive investigation ----------------------------------------

    def investigate(self) -> None:
        if self._guard_busy():
            return
        sel = self.app.results.selected()
        dlg = InvestigateDialog(self.app.root,
                                default_accession=sel.accession if sel else "")
        if not dlg.result:
            return
        acc = dlg.result["accession"]
        panel = [(f"panel_{a}", a)
                 for a in (s.strip() for s in dlg.result["panel_csv"].split(",")) if a]
        options = InvestigationOptions(
            panel=panel or list(InvestigationOptions().panel),
            run_foldseek=dlg.result["foldseek"],
            email=self.app.email,
        )
        seq = sel.sequence if (sel and sel.accession == acc) else ""

        def job(progress):
            result = Investigator(acc, options, candidate_sequence=seq).run(on_progress=progress)
            signals = assess_signals(result)
            verdict, pct = verdict_from_signals(signals)
            return result, verdict, pct

        if self.runner.run(f"investigation of {acc}", job, self._investigation_done,
                           self._task_failed,
                           lambda msg: self.app.status(f"Investigating {acc}: {msg}")):
            self.app.status(f"Investigating {acc}…")

    def _investigation_done(self, payload) -> None:
        result, verdict, pct = payload
        self.app.status(f"Investigation of {result.accession} done — verdict: {verdict} ({pct}/100).")

        def save() -> None:
            results_root = self.app.project_root / "results"
            results_root.mkdir(exist_ok=True)
            bundle = make_bundle_dir(results_root, f"investigate_{result.accession}")
            paths = write_investigation(bundle / "investigations" / result.accession, result)
            messagebox.showinfo("Investigation", "Saved:\n" + "\n".join(paths.values()))

        TextWindow(self.app.root,
                   f"Case file — {result.accession}  ({verdict}, {pct}/100)",
                   render_case_file(result),
                   save_name=f"case_file_{result.accession}.md",
                   extra_actions=[("Save artefacts…", save)])

    # ---- bundle augmentation --------------------------------------------

    # ------------------------------------------------------------ classify
    def classify_selection(self) -> None:
        """Classify the selected row into the ion-channel catalogue.

        The one Analysis-menu action that is specific to this project rather
        than inherited: it runs the same three-tier classifier the CLI and
        the S1 benchmark use, on whatever row is selected, and shows the full
        audit trail rather than just the answer. The row's gene symbol is
        displayed and never used to decide (hazard H15).
        """
        if self._guard_busy():
            return
        variant = self.app.results.selected()
        if variant is None:
            messagebox.showinfo("Classify", "Select a result row first.")
            return

        def work():
            from ..classify import ChannelQuery, ReferenceSet, classify
            from ..classify.motifs import FOUR_REPEAT_ANCHOR
            from ..classify.reference import fetch_uniprot_sequence
            from ..cli_channel import build_query
            acc = variant.accession
            query = (build_query(acc) if variant.source == "UniProt"
                     else ChannelQuery(accession=acc,
                                       sequence=variant.sequence,
                                       gene_symbol=variant.gene_symbol,
                                       species=variant.species))
            panel = (self.app.project_root / "results" / "s0_baseline"
                     / "reference_panel.fasta")
            refs = ReferenceSet.from_fasta(panel) if panel.exists() else None
            nav = fetch_uniprot_sequence(FOUR_REPEAT_ANCHOR.reference_uniprot) \
                if refs else ""
            return classify(query, refs, nav)

        def done(call):
            TextWindow(self.app.root, f"Classification — {call.query}",
                       call.text()
                       + ("\n\n(no reference panel at "
                          "results/s0_baseline/reference_panel.fasta — run "
                          "scripts/s0_catalogue_verify.py to enable the "
                          "reference tier)" if not call.margin else ""))

        self.runner.run(work, done, self._task_failed,
                        status="classifying against the catalogue…")
