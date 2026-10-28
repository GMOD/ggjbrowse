#!/bin/sh
# Prints the newer of @jbrowse/tables' `latest` and `next` releases on npm, by
# publish date, or nothing when neither is published. Prereleases go out
# under `next`, so `latest` alone can name one older than the bundle.
set -eu
curl -fsSL https://registry.npmjs.org/@jbrowse%2ftables | jq -r '
  .["dist-tags"] as $tags | .time as $time
  | [$tags.latest, $tags.next] | map(select(. != null))
  | if length > 0 then max_by($time[.]) else empty end'
