"""Git provenance for every file the dashboard serves.

Each panel on the site names the committed file it was read from, the last commit that changed that file and a GitHub
link to the file at that commit. Only files tracked by git are ever served; `require_clean()` refuses to export a file
with uncommitted local changes, so the public site can only show what is already in the public repository.
"""
from __future__ import annotations

import subprocess
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPO_URL = "https://github.com/gddickinson/ion_channel_genes"
#: Paths the dashboard reads from (git pathspecs).
WATCHED = ["results", "FINDINGS.md", "PUBLICATION_ROADMAP.md", "scripts/figure_notes.py",
           "scripts/figure_notes_late.py"]


def _git(*args: str) -> str | None:
    try:
        r = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout if r.returncode == 0 else None


@lru_cache(maxsize=1)
def tracked() -> frozenset[str]:
    """Every path tracked by git (repo-relative, '/' separators)."""
    out = _git("ls-files", "-z")
    return frozenset(p for p in (out or "").split("\0") if p)


@lru_cache(maxsize=1)
def last_commits() -> dict[str, dict]:
    """path -> {sha, date, subject} of the newest commit that touched it (one `git log` pass)."""
    out = _git("log", "--format=%x01%H%x09%cs%x09%s", "--name-only", "--", *WATCHED) or ""
    found: dict[str, dict] = {}
    for block in out.split("\x01")[1:]:
        lines = block.strip("\n").split("\n")
        sha, date, subject = (lines[0].split("\t") + ["", ""])[:3]
        for path in lines[1:]:
            if path and path not in found:
                found[path] = {"sha": sha, "date": date, "subject": subject}
    return found


@lru_cache(maxsize=1)
def modified() -> frozenset[str]:
    """Watched paths with uncommitted changes (staged, unstaged or untracked)."""
    out = _git("status", "--porcelain=v1", "-z", "-uall", "--", *WATCHED) or ""
    paths = set()
    for entry in out.split("\0"):
        if len(entry) > 3:
            paths.add(entry[3:])
    return frozenset(paths)


@lru_cache(maxsize=1)
def head() -> dict:
    out = (_git("log", "-1", "--format=%H%x09%cs") or "").strip()
    sha, _, date = out.partition("\t")
    return {"sha": sha, "commit": sha[:7], "date": date}


def of(path: str) -> dict:
    """Provenance record for one repo-relative path."""
    c = last_commits().get(path)
    rec = {"file": path, "committed": path in tracked(), "modified": path in modified()}
    if c:
        rec.update(commit=c["sha"][:7], date=c["date"], subject=c["subject"],
                   url=f"{REPO_URL}/blob/{c['sha']}/{path}")
    else:
        rec.update(commit=None, date=None, subject=None, url=f"{REPO_URL}/blob/main/{path}")
    return rec


def is_servable(path: str) -> bool:
    """A file may be served only if git tracks it."""
    return path in tracked() and (ROOT / path).is_file()


def require_clean(paths) -> None:
    """Raise if any path is untracked or has uncommitted changes (used by the static export)."""
    bad = sorted(p for p in paths if p not in tracked() or p in modified())
    if bad:
        raise RuntimeError("refusing to publish files that are not committed as-is: " + ", ".join(bad[:10])
                           + (f" (+{len(bad) - 10} more)" if len(bad) > 10 else ""))
