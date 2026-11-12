// Sortable, filterable, paged table for the committed TSVs (all rows already in the browser; capped server-side).
import { esc } from './core.js';

const NUM = /^-?(\d+\.?\d*|\.\d+)(e[-+]?\d+)?$/i;
export const isNum = s => NUM.test(String(s).trim());

// Sort key comparison: numbers numerically, blanks last, text alphabetically.
export function compare(a, b) {
  const ea = a === '' || a == null, eb = b === '' || b == null;
  if (ea || eb) return ea === eb ? 0 : ea ? 1 : -1;
  if (isNum(a) && isNum(b)) return parseFloat(a) - parseFloat(b);
  return String(a).localeCompare(String(b));
}

export function filterRows(rows, q) {
  const words = q.toLowerCase().split(/\s+/).filter(Boolean);
  if (!words.length) return rows;
  return rows.filter(r => { const t = r.join('\t').toLowerCase(); return words.every(w => t.includes(w)); });
}

// Render into host. opts: {pageSize, id}
export function tableView(host, columns, rows, { pageSize = 50, id = 'tv' } = {}) {
  const numCol = columns.map((_, j) => rows.length > 0 && rows.slice(0, 200).every(r => r[j] === '' || isNum(r[j])));
  let q = '', sortCol = -1, dir = 1, page = 0;
  host.innerHTML = `<div class="tv-bar"><div><label for="${id}-q">Filter rows (all words must match)</label>
    <input id="${id}-q" type="search" placeholder="type to filter"></div>
    <div class="pg"><button class="alt" data-p="-1" aria-label="Previous page">&lt;</button><span id="${id}-pi" aria-live="polite"></span>
    <button class="alt" data-p="1" aria-label="Next page">&gt;</button></div></div><div class="tv" id="${id}-t"></div>`;
  const out = host.querySelector(`#${id}-t`), pi = host.querySelector(`#${id}-pi`);
  function draw() {
    let view = filterRows(rows, q);
    if (sortCol >= 0) view = [...view].sort((a, b) => dir * compare(a[sortCol], b[sortCol]));
    const pages = Math.max(1, Math.ceil(view.length / pageSize));
    page = Math.min(Math.max(page, 0), pages - 1);
    const slice = view.slice(page * pageSize, (page + 1) * pageSize);
    pi.textContent = `${view.length ? page * pageSize + 1 : 0}-${page * pageSize + slice.length} of ${view.length}`;
    out.innerHTML = `<table><thead><tr>${columns.map((c, j) => `<th scope="col"><button data-c="${j}" aria-sort="${j === sortCol ? (dir > 0 ? 'ascending' : 'descending') : 'none'}">${esc(c)}${j === sortCol ? (dir > 0 ? ' ▲' : ' ▼') : ''}</button></th>`).join('')}</tr></thead>
      <tbody>${slice.map(r => `<tr>${columns.map((_, j) => `<td${numCol[j] ? ' class="num"' : ''} title="${esc(r[j])}">${esc(r[j])}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
  }
  host.querySelector(`#${id}-q`).addEventListener('input', e => { q = e.target.value; page = 0; draw(); });
  host.addEventListener('click', e => {
    const p = e.target.closest('button[data-p]'), c = e.target.closest('button[data-c]');
    if (p) { page += +p.dataset.p; draw(); }
    if (c) { const j = +c.dataset.c; dir = j === sortCol ? -dir : 1; sortCol = j; draw(); host.querySelector(`button[data-c="${j}"]`).focus(); }
  });
  draw();
}
