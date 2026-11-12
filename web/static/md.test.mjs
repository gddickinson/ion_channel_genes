// node web/static/md.test.mjs: the Markdown renderer on the constructs the reports and FINDINGS.md use.
import assert from 'node:assert/strict';
import { render, splitRow, inline } from './md.js';

assert.deepEqual(splitRow('| a | `x|y` | c\\|d |'), ['a', '`x|y`', 'c|d']);
const t = render('| h1 | h2 |\n|---|---:|\n| 1 | **b** |\n');
assert.ok(t.includes('<th scope="col">h1</th>') && t.includes('<td><strong>b</strong></td>'));
assert.ok(render('# Title').startsWith('<h1 id="title">Title</h1>'));
assert.ok(render('- one\n  more\n- two').includes('<li>one more</li><li>two</li>'));
assert.ok(render('1. a\n2. b').startsWith('<ol>'));
assert.ok(render('```\n<x>\n```').includes('<pre>&lt;x&gt;</pre>'));
assert.ok(render('> quoted').includes('<blockquote><p>quoted</p></blockquote>'));
assert.ok(inline('a <script>').includes('&lt;script&gt;'));
assert.ok(inline('see *Daphnia* and `a*b*c`').includes('<em>Daphnia</em>'));
assert.ok(inline('`a*b*c`').includes('<code>a*b*c</code>'));
assert.ok(inline('claim (pending: S13b)').includes('<mark class="pending">(pending: S13b)</mark>'));
assert.ok(inline('[r](report.md)', { link: h => '#/x/' + h }).includes('href="#/x/report.md"'));
assert.ok(inline('![f](figures/a.png)', { img: s => 'files/' + s }).includes('<img src="files/figures/a.png"'));
assert.ok(inline('![f](figures/a.png)').includes('<a href="figures/a.png">f</a>'));
assert.ok(render('p1 line\nline2\n\n---\n').includes('<p>p1 line line2</p>\n<hr>'));
console.log('md: ok');
