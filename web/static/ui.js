// Shared UI pieces: number formatting, ledger status badges (the repo's own markers), provenance lines.
import { esc, list } from './core.js';

export const fmt = v => (typeof v === 'number' ? v.toLocaleString('en-US') : esc(v ?? ''));

const STATE_TEXT = { completed: 'completed', in_progress: 'in progress', pending: 'pending', unknown: 'not in ledger' };
// One ledger task as a badge: "S5b completed 2026-09-29" / "S13b pending". Hover shows the ledger's task title.
export function badge(b) {
  const when = b.state === 'completed' && b.date ? ' ' + b.date : '';
  return `<span class="st st-${esc(b.state)}" title="${esc(b.title || '')}">${esc(b.id)} ${STATE_TEXT[b.state] || esc(b.state)}${when}</span>`;
}
export const badges = items => list(items, badge);

// "Is any task behind this panel unfinished?" -> a sentence quoting the ledger, else ''.
export function pendingNote(stage) {
  const open = (stage || []).filter(b => b.state !== 'completed');
  if (!open.length) return '';
  return `<p class="callout small">${open.map(b => `<b>${esc(b.id)}</b> is <b>${STATE_TEXT[b.state] || esc(b.state)}</b> in the roadmap ledger${b.title ? ` (${esc(b.title.slice(0, 160))}${b.title.length > 160 ? '...' : ''})` : ''}`).join('; ')}. Treat what is shown here as work in progress.</p>`;
}

// Provenance line for one committed file: path linked to GitHub at the last commit that changed it.
export function prov(p, label = 'Source') {
  if (!p) return '';
  const warn = !p.committed ? ' <span class="warn">(not committed: local only)</span>'
    : p.modified ? ' <span class="warn">(uncommitted local changes: not what the public repo holds)</span>' : '';
  const commit = p.commit ? ` · commit <code>${esc(p.commit)}</code>${p.date ? ` (${esc(p.date)})` : ''}` : '';
  return `${esc(label)}: <a href="${esc(p.url)}" target="_blank" rel="noopener">${esc(p.file)}</a>${commit}${warn}`;
}
export const provBlock = (...items) => `<p class="prov">${items.filter(Boolean).join('<br>')}</p>`;

// Percentage bar of ledger states.
export function meter(c, total) {
  const pct = n => (100 * n / Math.max(total, 1)).toFixed(1);
  return `<div class="meter" role="img" aria-label="${c.completed} completed, ${c.in_progress} in progress, ${c.pending} pending of ${total} tasks">
    <i class="m-c" style="width:${pct(c.completed)}%"></i><i class="m-i" style="width:${pct(c.in_progress)}%"></i><i class="m-p" style="width:${pct(c.pending)}%"></i></div>`;
}

export const tableHref = path => '#/table/' + path.split('/').map(encodeURIComponent).join('/');
export const treeHref = path => '#/tree/' + path.split('/').map(encodeURIComponent).join('/');
export const reportHref = path => '#/report/' + path.split('/').map(encodeURIComponent).join('/');
