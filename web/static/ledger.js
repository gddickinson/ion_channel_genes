// Progress: every roadmap ledger task with its status marker and its headline (quoted from the ledger).
import { route, api, app, esc, list } from './core.js';
import { badge, meter, prov, provBlock } from './ui.js';

route('ledger', async (_a, query) => {
  const led = await api('/ledger');
  const show = query.state || 'all';
  const tasks = led.tasks.filter(t => show === 'all' || t.state === show);
  const c = led.counts;
  app().innerHTML = `<p class="kicker">Roadmap ledger</p><h1>Progress</h1>
  <p class="lead">${c.completed} of ${led.total} tasks completed, ${c.in_progress} in progress, ${c.pending} pending.</p>
  ${meter(c, led.total)}
  <div class="seg" role="group" aria-label="Show tasks">${list(['all', 'completed', 'in_progress', 'pending'], s =>
    `<button class="${s === show ? 'on' : ''}" data-s="${s}" aria-pressed="${s === show}">${s.replace('_', ' ')}</button>`)}</div>
  <p class="note">${esc(led.note)}</p>
  <div class="scroll"><table class="ledger"><thead><tr><th scope="col">Task</th><th scope="col">What</th><th scope="col">Status</th></tr></thead>
  <tbody>${list(tasks, t => `<tr><td>${esc(t.id)}</td><td>${esc(t.title)}
    <div class="small mute">depends on ${esc(t.depends)} · <a href="#/tables?task=${esc(t.id)}">tables</a></div>
    ${t.results ? `<details class="hl"><summary>headline (ledger prose)</summary><div>${esc(t.results)}</div></details>` : ''}</td>
    <td>${badge(t)}</td></tr>`)}</tbody></table></div>
  ${provBlock(prov(led.provenance, 'Ledger'))}`;
  app().querySelector('.seg').addEventListener('click', e => {
    const b = e.target.closest('button[data-s]');
    if (b) location.hash = '#/ledger' + (b.dataset.s === 'all' ? '' : '?state=' + b.dataset.s);
  });
});
