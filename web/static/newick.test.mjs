// node web/static/newick.test.mjs: parse the committed trees listed in tests/fixtures/newick_fixture.json and compare
// with the values scripts/s7_newick.py produced (scripts/make_newick_fixture.py).
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import assert from 'node:assert/strict';
import { parse, summary, leaves, layout, support } from './newick.js';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
const fx = JSON.parse(readFileSync(join(ROOT, 'tests', 'fixtures', 'newick_fixture.json'), 'utf8'));
const close = (a, b) => Math.abs(a - b) <= 1e-9 + 1e-6 * Math.abs(b);

for (const t of fx.trees) {
  const root = parse(readFileSync(join(ROOT, t.path), 'utf8'));
  const s = summary(root);
  assert.equal(s.leaves, t.leaves, t.path + ' leaves');
  assert.equal(s.internal, t.internal, t.path + ' internal');
  assert.equal(s.ge95, t.ge95, t.path + ' ge95');
  assert.ok(close(s.totalLength, t.total_length), `${t.path} length ${s.totalLength} vs ${t.total_length}`);
  assert.deepEqual(leaves(root).slice(0, 5), t.first_leaves, t.path + ' first leaves');
  const L = layout(root);
  assert.equal(L.rows, t.leaves);
  assert.ok(L.maxX > 0);
}
assert.equal(support('85.1/98'), 98);
assert.equal(support(''), null);
assert.equal(support('abc'), null);
assert.throws(() => parse('(a,b'));
console.log(`newick: ${fx.trees.length} trees match scripts/s7_newick.py`);
