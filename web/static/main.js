// Entry point: registers the page modules (loaded on first visit), sets the mode line, starts the router.
import { start, route, lazy, app, $, isStatic } from './core.js';

const PAGES = {
  'home.js': ['home', 'about'], 'ledger.js': ['ledger'], 'census.js': ['census'], 'trees.js': ['trees', 'tree'],
  'figures.js': ['figures', 'figure'], 'tables.js': ['tables', 'table'], 'docs.js': ['reports', 'report', 'findings'],
};
lazy(Object.fromEntries(Object.entries(PAGES).flatMap(([file, names]) => names.map(n => [n, () => import('./' + file)]))));

route('notfound', () => { app().innerHTML = '<h1>Not found</h1><p><a href="#/">Back to the start</a></p>'; });

$('#mode').textContent = isStatic() ? 'static site, built from the committed results' : 'local app, reading the working tree';
start();
