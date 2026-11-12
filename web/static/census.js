// Census: the S5 presence matrix (family x species, genome verdicts) and the census status tables of each version.
import { route, api, app, esc, list } from './core.js';
import { fmt, badges, pendingNote, prov, provBlock, tableHref } from './ui.js';

// Colours on the repo's best-evidence-first order (scripts/s5r4_figures.py:VERDICT), lightened for a dark page.
const COLOUR = { present: '#2f8cff', genome_found: '#5fb0ff', genome_present: '#a7d0ff', genome_weak: '#8a8f9e',
  partial: '#555c78', gap: '#555c78', trace: '#ef9a90', absent: '#ff4d5e', no_locus_unrescued: '#1a2038', unmatched: '#c79bff' };

const VERSIONS = [
  ['results/census_v2/census_status.tsv', 'Census v2 (S2): every UniProtKB record with a pore signature, by status and confidence'],
  ['results/census_v3/v3_status.tsv', 'Census v3a (S3a): the same records after best-profile assignment'],
  ['results/panel_sweep/census_status.tsv', 'Census v3 over the panel proteomes (S3b)'],
  ['results/genome_sweep/census_v4_status.tsv', 'Census v4 (S5b): census v3 plus genome-only loci'],
];

function matrixSvg(m) {
  const cs = 12, left = 190, top = 118, W = left + m.families.length * cs + 10, H = top + m.species.length * cs + 10;
  let s = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" role="img" aria-label="Presence matrix: ${m.species.length} species by ${m.families.length} census families">`;
  m.families.forEach((f, j) => {
    const x = left + j * cs + cs / 2 + 3;
    s += `<text transform="translate(${x},${top - 6}) rotate(-70)">${esc(f.id)}</text>`;
    if (j && f.superfamily !== m.families[j - 1].superfamily) s += `<line x1="${left + j * cs - 1}" x2="${left + j * cs - 1}" y1="${top - 4}" y2="${H - 10}" stroke="#55618a" stroke-width="1"/>`;
  });
  m.species.forEach((sp, i) => {
    const y = top + i * cs;
    s += `<text x="${left - 6}" y="${y + 9.5}" text-anchor="end" font-style="italic">${esc(sp.name.length > 28 ? sp.name.slice(0, 27) + '…' : sp.name)}</text>`;
    if (i && sp.group !== m.species[i - 1].group) s += `<line x1="${left - 180}" x2="${W - 10}" y1="${y - 1}" y2="${y - 1}" stroke="#55618a" stroke-width="1"/>`;
    m.cells[i].forEach((c, j) => {
      if (!c) return;
      const v = m.verdicts[c[0]];
      s += `<rect class="cell" data-i="${i}" data-j="${j}" x="${left + j * cs}" y="${y}" width="${cs - 1.5}" height="${cs - 1.5}" fill="${COLOUR[v.id] || '#fff'}"><title>${esc(sp.name)} · ${esc(m.families[j].id)}: ${esc(v.label)}</title></rect>`;
    });
  });
  return s + '</svg>';
}

route('census', async () => {
  const m = await api('/matrix');
  const counts = m.verdicts.map((_, k) => m.cells.flat().filter(c => c && c[0] === k).length);
  app().innerHTML = `<p class="kicker">Census</p><h1>Presence matrix</h1>
  <p class="lead">Every census family (columns, grouped by superfamily) in every panel species (rows, grouped by
  lineage), coloured by the genome sweep's verdict for that cell. Click a cell for its counts.</p>
  ${pendingNote(m.stage)}
  <div class="legend">${list(m.verdicts, (v, k) => `<span><i style="background:${COLOUR[v.id] || '#fff'}"></i>${esc(v.id)}: ${esc(v.label)}${v.in_repo_scale ? '' : ' (not on the repo figure scale)'} <span class="mute">(${fmt(counts[k])})</span></span>`)}</div>
  <div class="matrix-wrap" id="mx">${matrixSvg(m)}</div>
  <div class="panel" id="cell" aria-live="polite"><p class="mute small">Select a cell to see its counts.</p></div>
  ${provBlock(prov(m.provenance), `Written by ${badges(m.stage)} · full table: <a href="${tableHref(m.provenance.file)}">open in the table explorer</a>`)}
  <h2>Census versions</h2>
  <p class="note">Each census version's status table, as committed. Version numbers and statuses are the project's own.</p>
  <div id="versions"></div>`;
  const sel = app().querySelector('#mx');
  sel.addEventListener('click', e => {
    const r = e.target.closest('rect.cell');
    if (!r) return;
    sel.querySelectorAll('rect.sel').forEach(x => x.classList.remove('sel'));
    r.classList.add('sel');
    const i = +r.dataset.i, j = +r.dataset.j, c = m.cells[i][j], sp = m.species[i], f = m.families[j];
    app().querySelector('#cell').innerHTML = `<h3>${esc(sp.name)} · ${esc(f.id)}</h3>
      <table><tbody>${list(m.cell_fields, (k, n) => `<tr><th scope="row">${esc(k)}</th><td>${n === 0 ? esc(m.verdicts[c[0]].id + ' (' + m.verdicts[c[0]].label + ')') : fmt(c[n])}</td></tr>`)}
      <tr><th scope="row">superfamily</th><td>${esc(f.superfamily)}</td></tr><tr><th scope="row">species group</th><td>${esc(sp.group)}</td></tr></tbody></table>`;
  });
  const host = app().querySelector('#versions');
  for (const [path, title] of VERSIONS) {
    const t = await api('/table/' + path);
    const div = document.createElement('div');
    div.className = 'panel';
    div.innerHTML = `<h3>${esc(title)}</h3><div class="scroll"><table><thead><tr>${list(t.columns, c => `<th scope="col">${esc(c)}</th>`)}</tr></thead>
      <tbody>${list(t.rows, r => `<tr>${list(r, v => `<td>${/^\d+$/.test(v) ? fmt(+v) : esc(v)}</td>`)}</tr>`)}</tbody></table></div>
      ${provBlock(prov(t.provenance), badges(t.stage))}`;
    host.appendChild(div);
  }
});
