// Small Markdown renderer for the stage reports and FINDINGS.md: headings, paragraphs, lists, tables, code blocks,
// quotes, rules, inline code/bold/italic/links/images. All text is escaped first; links and images go through the
// caller's resolvers. The repo's "(pending: task)" markers are highlighted.
import { esc } from './core.js';

// Split a table row on '|' that is neither escaped (\|) nor inside `code`.
export function splitRow(line) {
  const s = line.trim().replace(/^\|/, '').replace(/\|$/, '');
  const cells = []; let cur = '', code = false;
  for (let i = 0; i < s.length; i++) {
    const ch = s[i];
    if (ch === '\\' && s[i + 1] === '|') { cur += '|'; i++; continue; }
    if (ch === '`') code = !code;
    if (ch === '|' && !code) { cells.push(cur.trim()); cur = ''; continue; }
    cur += ch;
  }
  cells.push(cur.trim());
  return cells;
}

export function inline(text, { link = h => h, img = () => null } = {}) {
  const codes = [];
  let s = esc(text).replace(/`([^`]+)`/g, (_, c) => { codes.push(c); return `\u0000${codes.length - 1}\u0000`; });
  s = s.replace(/!\[([^\]]*)\]\(([^)\s]+)\)/g, (_, alt, src) => {
    const u = img(src.replace(/&amp;/g, '&'));
    return u ? `<img src="${esc(u)}" alt="${alt}" loading="lazy">` : `<a href="${esc(link(src.replace(/&amp;/g, '&')))}">${alt || 'image'}</a>`;
  });
  s = s.replace(/\[([^\]]+)\]\(([^)\s]+)\)/g, (_, t, h) => {
    const u = link(h.replace(/&amp;/g, '&'));
    const ext = /^https?:/.test(u);
    return `<a href="${esc(u)}"${ext ? ' target="_blank" rel="noopener"' : ''}>${t}</a>`;
  });
  s = s.replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>').replace(/(^|[^*\w])\*([^*\s][^*]*?)\*(?!\w)/g, '$1<em>$2</em>');
  s = s.replace(/\((pending:[^)]*)\)/g, '<mark class="pending">($1)</mark>');
  return s.replace(/\u0000(\d+)\u0000/g, (_, i) => `<code>${codes[+i]}</code>`);
}

export const slug = t => t.toLowerCase().replace(/<[^>]+>/g, '').replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');

export function render(md, opts = {}) {
  const lines = md.replace(/\r/g, '').split('\n');
  const out = []; let i = 0;
  const isTableSep = l => /^\s*\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)*\|?\s*$/.test(l);
  const isBlockStart = l => /^(#{1,6} |```|> |\s*[-*] |\s*\d+\. |\s*\|)/.test(l) || /^-{3,}\s*$/.test(l);
  while (i < lines.length) {
    const l = lines[i];
    if (!l.trim()) { i++; continue; }
    if (l.startsWith('```')) {
      const buf = []; i++;
      while (i < lines.length && !lines[i].startsWith('```')) buf.push(lines[i++]);
      i++; out.push(`<pre>${esc(buf.join('\n'))}</pre>`); continue;
    }
    const h = l.match(/^(#{1,6}) (.*)$/);
    if (h) { const n = h[1].length, body = inline(h[2], opts); out.push(`<h${n} id="${slug(body)}">${body}</h${n}>`); i++; continue; }
    if (/^(-{3,}|\*{3,})\s*$/.test(l)) { out.push('<hr>'); i++; continue; }
    if (l.trim().startsWith('|') && i + 1 < lines.length && isTableSep(lines[i + 1])) {
      const head = splitRow(l); i += 2; const rows = [];
      while (i < lines.length && lines[i].trim().startsWith('|')) rows.push(splitRow(lines[i++]));
      out.push(`<div class="tblwrap"><table><thead><tr>${head.map(c => `<th scope="col">${inline(c, opts)}</th>`).join('')}</tr></thead><tbody>${rows.map(r => `<tr>${r.map(c => `<td>${inline(c, opts)}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`);
      continue;
    }
    if (l.startsWith('>')) {
      const buf = [];
      while (i < lines.length && lines[i].startsWith('>')) buf.push(lines[i++].replace(/^>\s?/, ''));
      out.push(`<blockquote>${render(buf.join('\n'), opts)}</blockquote>`); continue;
    }
    const li = l.match(/^(\s*)([-*]|\d+\.) (.*)$/);
    if (li) {
      const ordered = /\d/.test(li[2]); const items = [];
      while (i < lines.length) {
        const m = lines[i].match(/^(\s*)([-*]|\d+\.) (.*)$/);
        if (m) { items.push({ depth: Math.floor(m[1].length / 2), text: m[3] }); i++; continue; }
        if (lines[i].trim() && /^\s+/.test(lines[i]) && items.length) { items[items.length - 1].text += ' ' + lines[i].trim(); i++; continue; }
        break;
      }
      const tag = ordered ? 'ol' : 'ul';
      out.push(`<${tag}>${items.map(it => `<li${it.depth ? ` style="margin-left:${it.depth * 1.4}em"` : ''}>${inline(it.text, opts)}</li>`).join('')}</${tag}>`);
      continue;
    }
    const buf = [l]; i++;
    while (i < lines.length && lines[i].trim() && !isBlockStart(lines[i])) buf.push(lines[i++]);
    out.push(`<p>${inline(buf.join(' '), opts)}</p>`);
  }
  return out.join('\n');
}
