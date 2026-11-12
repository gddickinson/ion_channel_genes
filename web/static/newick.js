// Newick reader + rectangular layout. parse() follows scripts/s7_newick.py:parse (same grammar: labels without quotes
// or whitespace, ':' branch lengths, internal labels = IQ-TREE UFBoot); newick.test.mjs checks it against that code.

export function parse(text) {
  const s = text.trim().replace(/;$/, '');
  let pos = 0;
  function label() {
    let start = pos;
    while (pos < s.length && !',():'.includes(s[pos])) pos++;
    const name = s.slice(start, pos);
    let length = null;
    if (s[pos] === ':') {
      start = ++pos;
      while (pos < s.length && !',()'.includes(s[pos])) pos++;
      length = parseFloat(s.slice(start, pos));
    }
    return [name, length];
  }
  function node() {
    const n = { name: '', length: null, children: [] };
    if (s[pos] === '(') {
      pos++;
      n.children.push(node());
      while (s[pos] === ',') { pos++; n.children.push(node()); }
      if (s[pos] !== ')') throw new Error(`expected ')' at ${pos}`);
      pos++;
    }
    [n.name, n.length] = label();
    return n;
  }
  const root = node();
  if (pos !== s.length) throw new Error(`trailing text at ${pos}`);
  return root;
}

// UFBoot from an IQ-TREE internal label ('98', or 'SH/UFBoot' '85.1/98'), as s7_newick.support.
export function support(name) {
  if (!name) return null;
  const v = Number(name.split('/').pop());
  return Number.isFinite(v) ? v : null;
}

export function leaves(n, out = []) {
  if (!n.children.length) out.push(n.name);
  n.children.forEach(c => leaves(c, out));
  return out;
}

export function summary(root) {
  let internal = 0, ge95 = 0, total = 0;
  (function walk(n, isRoot) {
    total += n.length || 0;
    if (n.children.length && !isRoot) { internal++; const v = support(n.name); if (v != null && v >= 95) ge95++; }
    n.children.forEach(c => walk(c, false));
  })(root, true);
  return { leaves: leaves(root).length, internal, ge95, totalLength: total };
}

// x = distance from root (branch lengths, missing = 0), y = leaf order; internal y = mean of children's.
export function layout(root) {
  let row = 0, maxX = 0;
  const nodes = [];
  (function walk(n, x, parent) {
    const me = { n, x, parent, y: 0 };
    nodes.push(me);
    maxX = Math.max(maxX, x);
    if (!n.children.length) me.y = row++;
    else {
      const kids = n.children.map(c => walk(c, x + (c.length || 0), me));
      me.y = (kids[0].y + kids[kids.length - 1].y) / 2;
    }
    return me;
  })(root, 0, null);
  return { nodes, rows: row, maxX };
}
