// Static-site data layer: answers api(path) from the JSON files scripts/build_static.py exported (api/<path>.json).
// The site is read-only, so every endpoint is a plain file; nothing is computed here.
const files = new Map();
export async function staticApi(path, body) {
  if (body) throw new Error('The static site is read-only.');
  const clean = path.split('?')[0].replace(/^\/+|\/+$/g, '');
  if (!files.has(clean)) {
    files.set(clean, fetch('api/' + clean.split('/').map(encodeURIComponent).join('/') + '.json')
      .then(r => { if (!r.ok) throw new Error(`${r.status} ${clean}`); return r.json(); }));
  }
  return files.get(clean);
}
