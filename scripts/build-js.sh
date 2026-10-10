#!/bin/sh
# Rebuilds inst/js/jbrowse-tables.js from a jbrowse-components checkout
# (../jbrowse-components, or $JBROWSE_COMPONENTS), and names on its first line
# the commit it was built from. Run `pnpm install` in that checkout first.
set -eu
jbc=${JBROWSE_COMPONENTS:-../jbrowse-components}
out=inst/js/jbrowse-tables.js
node "$jbc/packages/tables/scripts/build-bundle.mjs" "$out.tmp"
commit=$(git -C "$jbc" rev-parse HEAD)
ref=${JBROWSE_REF:-$(git -C "$jbc" describe --tags --exact-match 2>/dev/null || echo unreleased)}
version=$(node -p "require('$jbc/packages/tables/package.json').version")
{
  printf '// @jbrowse/tables %s (%s), GMOD/jbrowse-components %s\n' "$version" "$ref" "$commit"
  cat "$out.tmp"
} > "$out"
rm "$out.tmp"
echo "built $out from $commit"
