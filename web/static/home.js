// Home (hero, headline numbers read from committed summaries, ledger progress, entry cards) and About.
import { route, api, app, esc, list } from './core.js';
import { fmt, badge, meter, prov, provBlock, tableHref } from './ui.js';

function tile(h) {
  if (h.pending) {
    return `<a class="tile pending" href="#/ledger"><b>pending</b><span>${esc(h.label)}</span>
      <span class="src">${badge(h.stage)} no committed result yet</span></a>`;
  }
  const p = h.provenance;
  return `<a class="tile" href="${tableHref(p.file)}" title="${esc(p.file)}: ${esc(h.source)}">
    <b>${fmt(h.value)}${h.of != null ? ` <small>of ${fmt(h.of)}</small>` : ''}</b><span>${esc(h.label)}</span>
    <span class="src">${esc(p.file.replace(/^results\//, ''))} · ${esc(h.source)} · ${esc(p.commit || '?')}</span>
    <span class="src">${badge(h.stage)}</span></a>`;
}

const DOORS = [
  ['#/ledger', 'Progress', 'Every roadmap task with its status (completed, in progress, pending) and its headline, from the ledger.', 'tasks'],
  ['#/census', 'Presence matrix', 'Every census family in every panel species, with the genome sweep\'s verdict for each cell.', null],
  ['#/trees', 'Trees', 'The committed family and superfamily trees (IQ-TREE, UFBoot support), drawn in the browser.', 'trees'],
  ['#/figures', 'Figures', 'The headline figure of each task with its plain-English caption.', 'figures'],
  ['#/tables', 'Tables', 'Every committed results table, sortable and filterable, with its source commit.', 'tables'],
  ['#/reports', 'Stage reports', 'The report each stage renders from its own committed tables.', 'reports'],
  ['#/findings', 'Findings', 'The biological story task by task, with unconfirmed claims marked pending.', null],
];

route('home', async () => {
  const [info, hl, led] = await Promise.all([api('/info'), api('/headlines'), api('/ledger')]);
  const c = led.counts;
  app().innerHTML = `
  <section class="hero">
    <p class="kicker">Census · Classification · Phylogeny</p>
    <h1 class="glow">Ion channel census</h1>
    <p class="lead">Identify, classify and reconstruct the phylogeny of every ion channel. This is a read-only dashboard
    of the project's committed results: every number below is read from a committed file when the site is built, and
    each one names its file and commit.</p>
    <a class="btn" href="#/ledger">Progress</a><a class="btn alt" href="#/census">Presence matrix</a><a class="btn alt" href="#/trees">Trees</a>
  </section>
  <h2>Headline numbers</h2>
  <p class="note">Click a number to open the table or summary it was read from.</p>
  <div class="tiles">${list(hl.headlines, tile)}</div>
  <h2>Roadmap progress</h2>
  <div class="panel">${meter(c, led.total)}
    <p class="small"><span class="ok">${c.completed} completed</span> · <span style="color:var(--accent3)">${c.in_progress} in progress</span> ·
    <span style="color:var(--warn)">${c.pending} pending</span> of ${led.total} ledger tasks.</p>
    ${provBlock(prov(led.provenance, 'Ledger'))}</div>
  <h2>Explore</h2>
  <div class="grid wide">${list(DOORS, ([href, t, d, k]) => `<a class="card" href="${href}"><b>${t}${k ? ` <small style="display:inline">${fmt(info.counts[k])}</small>` : ''}</b><span>${d}</span></a>`)}</div>
  ${provBlock(`Site built from commit <code>${esc(info.head.commit || '?')}</code>${info.head.date ? ` (${esc(info.head.date)})` : ''} of <a href="${esc(info.repo)}">${esc(info.repo.replace('https://', ''))}</a>`)}`;
});

route('about', async () => {
  const info = await api('/info');
  app().innerHTML = `<p class="kicker">About</p><h1>About this dashboard</h1>
  <div class="panel"><p><b>The project.</b> A census, classification and phylogeny of the ion channels as a whole: not
  one family but the independent superfamilies and human pore-forming genes the term covers, plus their relatives
  across the tree of life. Three deliverables, in order: a census (what exists, in a declared search space), a
  classification (what each thing is, by positive tests with an audit trail) and a phylogeny (a forest of
  within-family and within-superfamily trees plus a structural network, because the superfamilies are not homologous
  to one another).</p></div>
  <h2>Where the numbers come from</h2>
  <div class="panel"><ul>
    <li>Every number on the site is read from a file under <code>results/</code> that git tracks, at build time. Nothing
    is typed in by hand and nothing is recomputed.</li>
    <li>Each panel names its source file and the last commit that changed it, linked to that file on GitHub.</li>
    <li>Task status comes from the roadmap ledger (<code>PUBLICATION_ROADMAP.md</code>): <span class="st st-completed">completed</span>
    <span class="st st-in_progress">in progress</span> <span class="st st-pending">pending</span>. Panels whose task is not finished say so.</li>
    <li>Claims in <code>FINDINGS.md</code> that wait on unfinished work carry the repo's own <mark class="pending">(pending: task)</mark> marker.</li>
    <li>The ledger headlines, the findings and the stage reports are the project's own prose and are shown as written.</li>
    <li>Large tables are shown up to their first ${fmt(1500)} rows on this site, with a link to the full file.
    Bulk data (proteomes, genomes, HMM databases, structures) is not in the repository and is not shown.</li>
    <li>The static export refuses any file with uncommitted changes, so the public site shows only what the public
    repository holds.</li></ul></div>
  <h2>Run it yourself</h2>
  <div class="panel"><pre>git clone https://github.com/gddickinson/ion_channel_genes
cd ion_channel_genes
pip install -r requirements-web.txt
python run_web.py                         # http://127.0.0.1:8180
python scripts/build_static.py --out site # the static site, as published</pre>
  <p class="small">The census tool itself (catalogue, classifier, phylogeny, GUI) is described in the README:
  <code>python3 run.py --catalogue</code>, <code>python3 run.py --classify KCNQ1</code>, <code>python3 scripts/selftest.py</code>.</p></div>
  <h2>Credits</h2>
  <div class="panel"><p>By George Dickinson (MIT licence). The app, figure style, dashboard, manuscript tooling and
  session protocol were ported from the IP3 receptor project, itself ported from the PIEZO project; no result was.</p>
  ${provBlock(`Built from commit <code>${esc(info.head.commit || '?')}</code> · ${fmt(info.counts.tables)} tables · ${fmt(info.counts.trees)} trees · ${fmt(info.counts.reports)} reports · ${fmt(info.counts.figures)} figures`)}</div>`;
});
