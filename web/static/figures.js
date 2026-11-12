// Figures: the headline figure of each task (committed PNGs) with its caption from scripts/figure_notes.py.
import { route, api, app, esc, list } from './core.js';
import { badges, pendingNote, prov, provBlock } from './ui.js';

const src = path => 'files/' + path.split('/').map(encodeURIComponent).join('/');

route('figures', async () => {
  const d = await api('/figures');
  app().innerHTML = `<p class="kicker">Results in figures</p><h1>Figures</h1>
  <p class="lead">${d.figures.length} headline figures, one or more per completed task, each drawn by the task's own
  figure script from its committed tables. Captions are the project's own (<code>scripts/figure_notes.py</code>).</p>
  <div class="grid wide">${list(d.figures, (f, i) => `<a class="card fig" href="#/figure/${i}">
    <img src="${src(f.path)}" alt="${esc(f.task)}: ${esc(f.title)}" loading="lazy" width="660" height="440">
    <b>${esc(f.task)} · ${esc(f.title)}</b><small>${badges(f.stage)}</small></a>`)}</div>
  ${provBlock(...d.captions.map(c => prov(c, 'Captions')))}`;
});

route('figure', async args => {
  const d = await api('/figures');
  const i = Math.max(0, Math.min(d.figures.length - 1, +args[0] || 0)), f = d.figures[i];
  app().innerHTML = `<p class="kicker"><a href="#/figures">Figures</a> · ${esc(f.task)}</p><h1>${esc(f.title)}</h1>
  ${pendingNote(f.stage)}
  <div class="figfull"><a href="${src(f.path)}" target="_blank" rel="noopener"><img src="${src(f.path)}" alt="${esc(f.task)}: ${esc(f.title)}"></a></div>
  <div class="panel"><h3>What it shows</h3><p>${esc(f.shows)}</p><h3>How it was made</h3><p>${esc(f.how)}</p>
  <h3>How to read it</h3><p>${esc(f.read)}</p>
  ${provBlock(prov(f.provenance, 'Figure'), `Drawn by <code>${esc(f.script)}</code> · ${badges(f.stage)}`, ...d.captions.map(c => prov(c, 'Caption')))}</div>
  <p>${i > 0 ? `<a class="btn alt" href="#/figure/${i - 1}">Previous</a>` : ''}${i < d.figures.length - 1 ? `<a class="btn alt" href="#/figure/${i + 1}">Next</a>` : ''}</p>`;
});
