// Usage: node scripts/credit-bundle.mjs <entry> [out]
// Bundles <entry> (@jbrowse/tables' src/bundle.ts) as build-bundle.mjs does
// and writes to [out], default inst/COPYRIGHTS, the name, version, license
// and copyright of each package with code in the bundle.
import { existsSync, readFileSync, readdirSync, writeFileSync } from 'node:fs'
import { builtinModules, createRequire } from 'node:module'
import path from 'node:path'

const entry = path.resolve(process.argv[2])
const out = process.argv[3] ?? 'inst/COPYRIGHTS'
const { build } = createRequire(entry)('esbuild')
const builtins = new Set([
  ...builtinModules,
  ...builtinModules.map(m => `node:${m}`),
])

const { metafile } = await build({
  entryPoints: [entry],
  bundle: true,
  format: 'iife',
  platform: 'neutral',
  mainFields: ['module', 'main'],
  conditions: ['import', 'default'],
  target: 'es2020',
  minify: true,
  metafile: true,
  write: false,
  logLevel: 'warning',
  absWorkingDir: path.dirname(entry),
  define: { 'process.env.NODE_ENV': '"production"' },
  plugins: [
    {
      name: 'node-builtins',
      setup(b) {
        b.onResolve({ filter: /.*/ }, args =>
          builtins.has(args.path)
            ? { path: args.path, namespace: 'empty' }
            : undefined,
        )
        b.onLoad({ filter: /.*/, namespace: 'empty' }, () => ({
          contents: 'module.exports = {}',
        }))
      },
    },
  ],
})

function packageDir(file) {
  let dir = path.dirname(path.resolve(path.dirname(entry), file))
  while (dir !== path.dirname(dir)) {
    const json = path.join(dir, 'package.json')
    if (existsSync(json) && JSON.parse(readFileSync(json, 'utf8')).name) {
      return dir
    }
    dir = path.dirname(dir)
  }
}

function licenseOf(pkg) {
  const l = pkg.license ?? pkg.licenses
  if (typeof l === 'string') {
    return l
  }
  return [l]
    .flat()
    .map(x => x?.type ?? x)
    .filter(Boolean)
    .join(' OR ')
}

function authorOf(pkg) {
  const a = pkg.author
  return typeof a === 'string' ? a.replace(/\s*[<(].*$/, '') : a?.name
}

function licenseFile(dir) {
  const file = readdirSync(dir).find(f => /^(licen[cs]e|copying)/i.test(f))
  return file ? readFileSync(path.join(dir, file), 'utf8') : ''
}

const isCopyright = l =>
  /^(Copyright\b|\(c\) \d{4})/.test(l) && !l.includes('[yyyy]')

function copyrightOf(dir, pkg) {
  const lines = licenseFile(dir)
    .split('\n')
    .map(l => l.trim())
    .filter(isCopyright)
  return [...new Set(lines)].join('; ') || authorOf(pkg) || 'its authors'
}

const packages = new Map()
for (const [file, { bytesInOutput }] of Object.entries(
  Object.values(metafile.outputs)[0].inputs,
)) {
  const dir = bytesInOutput > 0 && packageDir(file.replace(/^[^:]+:/, ''))
  if (dir) {
    const pkg = JSON.parse(readFileSync(path.join(dir, 'package.json'), 'utf8'))
    packages.set(`${pkg.name}@${pkg.version}`, { dir, pkg })
  }
}
const sorted = [...packages.values()].sort(
  (a, b) =>
    a.pkg.name.localeCompare(b.pkg.name) ||
    a.pkg.version.localeCompare(b.pkg.version),
)

const texts = new Map()
for (const { dir, pkg } of sorted) {
  const id = licenseOf(pkg)
  const text = licenseFile(dir)
    .split('\n')
    .filter(l => !isCopyright(l.trim()))
    .join('\n')
    .replaceAll(/\n{3,}/g, '\n\n')
    .trim()
  if (!id.includes('Apache') && text && !texts.has(id)) {
    texts.set(id, `${id}, as in ${pkg.name}'s LICENSE:\n\n${text}`)
  }
}

const ascii = s =>
  s
    .replaceAll('©', '(c)')
    .normalize('NFKD')
    .replaceAll(/[^\n\x20-\x7E]/g, '')
const entries = sorted.map(
  ({ dir, pkg }) =>
    `- ${pkg.name} ${pkg.version}\n  License: ${licenseOf(pkg)}\n  Copyright: ${copyrightOf(dir, pkg)}`,
)
const rule = `\n\n${'-'.repeat(72)}\n\n`

writeFileSync(
  out,
  ascii(`inst/js/jbrowse-tables.js bundles @jbrowse/tables (Apache-2.0, copyright
the JBrowse authors) with the packages below, minified.

- Source: https://github.com/GMOD/jbrowse-components/tree/main/packages/tables
- Build: packages/tables/scripts/build-bundle.mjs there, at the commit or
  release the bundle's first line names
- Each package's source is on npm under its name and version
- The Apache-2.0 text: https://www.apache.org/licenses/LICENSE-2.0

${entries.join('\n')}${rule}${[...texts.values()].join(rule)}
`),
)
console.log(`credited ${entries.length} packages in ${out}`)
