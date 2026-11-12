// Trees: index of the committed treefiles with their statistics, and a rectangular tree viewer (SVG).
import { route, api, app, esc, list } from './core.js';
import { fmt, badges, pendingNote, prov, provBlock, treeHref, tableHref } from './ui.js';
import { parse, layout, support, summary } from './newick.js';

const COLS = { tier1: ['superfamily', 'n_ingroup', 'model', 'root_rule', 'root_ufboot', 'ufboot_median', 'frac_ge95'],
  reroot: ['n_ingroup', 'model', 'outgroup_monophyletic', 'root_ufboot', 'ufboot_median', 'frac_ge95'],
  tier2: ['kind', 'n_families', 'model', 'rooted', 'root_ufboot', 'ufboot_median', 'frac_ge95'], other: [] };

route('trees', async () => {
  const idx = await api('/trees');
  const groups = Object.entries(idx.groups);
  app().innerHTML = `<p class="kicker">Phylogeny</p><h1>Trees</h1>
  <p class="lead">${fmt(idx.trees.length)} committed maximum-likelihood trees (IQ-TREE 2, UFBoot support on internal
  nodes). The forest has no single tree of ion channels: tier 1 is within a family, tier 2 within a superfamily.</p>
  ${list(groups, ([g, title]) => {
    const ts = idx.trees.filter(t => t.group === g);
    if (!ts.length) return '';
    const cols = COLS[g] || [];
    return `<h2>${esc(title)} <span class="mute small">(${ts.length})</span></h2>
    <p class="note">${ts[0].stats_table ? `Statistics columns from <a href="${tableHref(ts[0].stats_table)}">${esc(ts[0].stats_table)}</a>. ` : ''}Written by ${badges(ts[0].stage)}</p>
    <div class="scroll"><table><thead><tr><th scope="col">tree</th><th scope="col">tips</th>${list(cols, c => `<th scope="col">${esc(c)}</th>`)}</tr></thead>
    <tbody>${list(ts, t => `<tr><td><a href="${treeHref(t.path)}">${esc(t.name)}</a></td><td>${fmt(t.tips)}</td>${list(cols, c => `<td>${esc(t.stats ? t.stats[c] : '')}</td>`)}</tr>`)}</tbody></table></div>`;
  })}`;
});

function svgTree(root, { showSupport, hit }) {
  const L = layout(root), rowH = 13, byNode = new Map(L.nodes.map(d => [d.n, d])), left = 10, labelW = 260;
  const width = Math.max(720, Math.min(1400, 300 + L.rows * 2)), plotW = width - left - labelW;
  const sx = x => left + (L.maxX ? x / L.maxX : 0) * plotW, sy = y => 14 + y * rowH;
  const H = sy(L.rows) + 10;
  let s = `<svg width="${width}" height="${H}" viewBox="0 0 ${width} ${H}" role="img" aria-label="Tree with ${L.rows} tips">`;
  for (const d of L.nodes) {
    if (d.parent) {
      const v = d.n.children.length ? support(d.n.name) : null;
      const col = v == null ? '#8fa3c8' : v >= 95 ? '#3df5e0' : v >= 70 ? '#8fa3c8' : '#ff5fd2';
      s += `<path d="M${sx(d.parent.x)},${sy(d.y)}H${sx(d.x)}" stroke="${col}" stroke-width="1.2" fill="none"/>`;
    }
    if (d.n.children.length) {
      const ys = d.n.children.map(c => byNode.get(c).y);
      s += `<path d="M${sx(d.x)},${sy(Math.min(...ys))}V${sy(Math.max(...ys))}" stroke="#8fa3c8" stroke-width="1" fill="none"/>`;
      if (showSupport && d.parent && d.n.name) s += `<text class="sup" x="${sx(d.x) - 2}" y="${sy(d.y) - 2}" text-anchor="end">${esc(d.n.name)}</text>`;
    } else {
      const isHit = hit && d.n.name.toLowerCase().includes(hit);
      s += `<text x="${sx(d.x) + 4}" y="${sy(d.y) + 3.5}"${isHit ? ' class="hit"' : ''}>${esc(d.n.name)}</text>`;
    }
  }
  return s + '</svg>';
}

route('tree', async args => {
  const path = args.join('/');
  const t = await api('/tree/' + path);
  const root = parse(t.newick), sm = summary(root);
  app().innerHTML = `<p class="kicker"><a href="#/trees">Trees</a> · ${esc(t.group)}</p><h1>${esc(t.name)}</h1>
  ${pendingNote(t.stage)}
  <p class="lead">${fmt(sm.leaves)} tips · ${fmt(sm.internal)} internal edges, ${fmt(sm.ge95)} with UFBoot ≥ 95 (counted from the treefile).</p>
  <div class="row"><div><label for="tr-q">Highlight tips containing</label><input id="tr-q" type="search" placeholder="e.g. Homsap, Q9UHC3"></div>
  <div style="flex:0 0 auto"><label><input type="checkbox" id="tr-s"> show UFBoot support</label></div></div>
  <p class="small mute">Branch colour: <span style="color:#3df5e0">UFBoot ≥ 95</span> · <span style="color:#8fa3c8">70-94 or tip</span> ·
  <span style="color:#ff5fd2">below 70</span>. Drawn as stored in the file (rooting as written by the stage); x = substitutions per site.</p>
  <div class="tree-wrap" id="tr"></div>
  ${t.stats ? `<h2>Statistics row</h2><div class="scroll"><table><tbody>${list(Object.entries(t.stats), ([k, v]) => `<tr><th scope="row">${esc(k)}</th><td class="mono">${esc(v)}</td></tr>`)}</tbody></table></div>` : ''}
  ${provBlock(prov(t.provenance, 'Tree'), t.stats_provenance ? prov(t.stats_provenance, 'Statistics') : '', 'Written by ' + badges(t.stage))}`;
  const host = app().querySelector('#tr'), q = app().querySelector('#tr-q'), cb = app().querySelector('#tr-s');
  const draw = () => { host.innerHTML = svgTree(root, { showSupport: cb.checked, hit: q.value.trim().toLowerCase() }); };
  q.addEventListener('input', draw); cb.addEventListener('change', draw);
  draw();
});
