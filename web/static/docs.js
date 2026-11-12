// Stage reports (results/**/report*.md) and FINDINGS.md, rendered from the committed Markdown.
import { route, api, app, esc, list } from './core.js';
import { badges, pendingNote, prov, provBlock, tableHref, treeHref, reportHref } from './ui.js';
import { render } from './md.js';

function join(dir, rel) {
  const parts = dir ? dir.split('/') : [];
  for (const p of rel.split('/')) { if (p === '..') parts.pop(); else if (p && p !== '.') parts.push(p); }
  return parts.join('/');
}

// Link/image resolvers for a Markdown file at `path`: repo files open in the matching page, others on GitHub.
async function resolvers(path, p) {
  const [tb, tr, rp, fg] = await Promise.all([api('/tables'), api('/trees'), api('/reports'), api('/figures')]);
  const tables = new Set(tb.tables.map(t => t.path)), trees = new Set(tr.trees.map(t => t.path));
  const reports = new Set(rp.reports.map(r => r.path)), figs = new Set(fg.figures.map(f => f.path));
  const dir = path.includes('/') ? path.slice(0, path.lastIndexOf('/')) : '';
  const blob = p.url.replace(/\/blob\/([^/]+)\/.*$/, '/blob/$1/');
  const local = h => join(dir, h.split('#')[0]);
  return {
    link: h => {
      if (/^(https?:|mailto:)/.test(h)) return h;
      if (h.startsWith('#')) return '#/' + location.hash.replace(/^#\/?/, '');
      const f = local(h);
      if (tables.has(f)) return tableHref(f);
      if (trees.has(f)) return treeHref(f);
      if (reports.has(f)) return reportHref(f);
      if (f === 'FINDINGS.md') return '#/findings';
      return blob + f;
    },
    img: h => (figs.has(local(h)) ? 'files/' + local(h) : null),
  };
}

route('reports', async () => {
  const d = await api('/reports');
  app().innerHTML = `<p class="kicker">Data</p><h1>Stage reports</h1>
  <p class="lead">Each stage renders its report purely from its own committed tables (decision D13), so report and
  data cannot drift. ${d.reports.length} reports.</p>
  <div class="scroll"><table><thead><tr><th scope="col">Report</th><th scope="col">File</th><th scope="col">Task</th></tr></thead>
  <tbody>${list(d.reports, r => `<tr><td><a href="${reportHref(r.path)}">${esc(r.title)}</a></td><td class="mono small">${esc(r.path)}</td><td>${badges(r.stage)}</td></tr>`)}</tbody></table></div>`;
});

route('report', async args => {
  const path = args.join('/');
  const r = await api('/report/' + path);
  const opts = await resolvers(path, r.provenance);
  app().innerHTML = `<p class="kicker"><a href="#/reports">Stage reports</a> · ${esc(path)}</p>
  ${pendingNote(r.stage)}
  <article class="md">${render(r.markdown, opts)}</article>
  ${provBlock(prov(r.provenance, 'Report'), 'Written by ' + badges(r.stage))}`;
});

route('findings', async () => {
  const f = await api('/findings');
  const opts = await resolvers('FINDINGS.md', f.provenance);
  const html = render(f.markdown, opts);
  const toc = [...html.matchAll(/<h([23]) id="([^"]+)">(.*?)<\/h\1>/g)];
  app().innerHTML = `<p class="kicker">Narrative</p>
  <div class="callout small">${esc(f.note)} <b>${f.pending_markers.length}</b> pending markers in the file.</div>
  <details class="panel"><summary>Contents (${toc.length} sections)</summary><div class="toc">${list(toc, m => `<a href="#/findings?keep=1" data-to="${m[2]}">${m[3].replace(/<[^>]+>/g, '')}</a>`)}</div></details>
  <article class="md">${html}</article>
  ${provBlock(prov(f.provenance))}`;
  app().querySelector('.toc').addEventListener('click', e => {
    const a = e.target.closest('a[data-to]');
    if (!a) return;
    e.preventDefault();
    document.getElementById(a.dataset.to)?.scrollIntoView();
  });
});
