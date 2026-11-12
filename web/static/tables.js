// Tables: index of every committed results table, and one table in the sortable/filterable viewer.
import { route, api, app, esc, list } from './core.js';
import { fmt, badges, pendingNote, prov, provBlock, tableHref } from './ui.js';
import { tableView } from './table_view.js';

route('tables', async (_a, query) => {
  const idx = await api('/tables');
  const led = await api('/ledger');
  const state = Object.fromEntries(led.tasks.map(t => [t.id, t]));
  app().innerHTML = `<p class="kicker">Data</p><h1>Tables</h1>
  <p class="lead">${fmt(idx.tables.length)} committed results tables (TSV) and small JSON summaries under <code>results/</code>.
  Tables longer than ${fmt(idx.max_rows)} rows show their first ${fmt(idx.max_rows)} here, with a link to the full file.</p>
  <div class="row"><div><label for="tb-q">Find a table (name, folder or task)</label>
  <input id="tb-q" type="search" value="${esc(query.task || '')}" placeholder="e.g. cells, phylogeny, S5b"></div></div>
  <div id="tb-out"></div>`;
  const out = app().querySelector('#tb-out');
  const draw = q => {
    const words = q.toLowerCase().split(/\s+/).filter(Boolean);
    const hits = idx.tables.filter(t => { const s = `${t.path} ${t.tasks.join(' ')}`.toLowerCase(); return words.every(w => s.includes(w)); });
    const dirs = [...new Set(hits.map(t => t.dir))];
    out.innerHTML = hits.length ? list(dirs, d => `<h3>results/${esc(d)}</h3><div class="scroll"><table><thead><tr>
      <th scope="col">table</th><th scope="col">rows</th><th scope="col">task</th></tr></thead><tbody>
      ${list(hits.filter(t => t.dir === d), t => `<tr><td><a href="${tableHref(t.path)}">${esc(t.name)}</a></td><td>${fmt(t.rows)}</td>
      <td>${badges(t.tasks.map(id => state[id] || { id, state: 'unknown' }))}</td></tr>`)}</tbody></table></div>`)
      : '<p class="mute">No table matches.</p>';
  };
  app().querySelector('#tb-q').addEventListener('input', e => draw(e.target.value));
  draw(query.task || '');
});

route('table', async args => {
  const path = args.join('/');
  const t = await api('/table/' + path);
  app().innerHTML = `<p class="kicker"><a href="#/tables">Tables</a> · results/${esc(t.dir)}</p><h1>${esc(t.name)}</h1>
  ${pendingNote(t.stage)}
  <p class="lead">${fmt(t.total_rows)} ${t.kind === 'json' ? 'keys (JSON summary, flattened)' : 'rows'} · ${fmt(t.columns.length)} columns.
  ${t.truncated ? `<span style="color:var(--warn)">Showing the first ${fmt(t.rows.length)} rows in file order;
  <a href="${esc(t.provenance.url)}" target="_blank" rel="noopener">the full table is on GitHub</a>.</span>` : ''}</p>
  <div id="tv"></div>
  ${provBlock(prov(t.provenance), 'Written by ' + (t.stage.length ? badges(t.stage) : 'no ledger task mapped'))}`;
  tableView(app().querySelector('#tv'), t.columns, t.rows, { id: 'tv' });
});
