#!/bin/sh
# Replaces inst/js/jbrowse-tables.js with the bundle a published @jbrowse/tables
# release carries (its dist/jbrowse-tables.js), named on the first line.
# Usage: scripts/take-release.sh [version], the latest by default.
# The npm release names no commit, so inst/COPYRIGHTS stays as it was.
set -eu
cd "$(dirname "$0")/.."
version=${1:-$(npm view @jbrowse/tables version)}
tmp=$(mktemp -d)
npm pack --silent --pack-destination "$tmp" "@jbrowse/tables@$version" >/dev/null
tar -xzf "$tmp"/*.tgz -C "$tmp" package/dist/jbrowse-tables.js
{
  printf '// @jbrowse/tables %s (npm)\n' "$version"
  cat "$tmp/package/dist/jbrowse-tables.js"
} > inst/js/jbrowse-tables.js
rm -rf "$tmp"
echo "took @jbrowse/tables $version"
echo "inst/COPYRIGHTS still credits the last build-js.sh bundle" >&2
