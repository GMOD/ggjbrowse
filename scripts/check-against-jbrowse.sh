#!/bin/sh
# The one command jbrowse-components' downstream canary and this repo's own
# bundle job both run: rebuild the engine from a jbrowse-components checkout
# ($JBROWSE_COMPONENTS, default ../jbrowse-components, already `pnpm
# install`ed), install the package with it, and run the tests on it.
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
cd "$here"
scripts/build-js.sh
lib=$(mktemp -d)
R CMD INSTALL -l "$lib" .
R_LIBS="$lib" Rscript -e 'testthat::test_dir("tests/testthat",
  package = "ggjbrowse", load_package = "installed", stop_on_failure = TRUE)'
