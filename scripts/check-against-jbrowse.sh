#!/bin/sh
# The one command jbrowse-components' downstream canary runs: rebuild the
# engine from a jbrowse-components checkout
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
# The R tutorials pull sections of the website tutorials, and a website edit
# that removes text a directive matches fails their expansion. Expanding them
# against this checkout's own pages fails the night of that edit, rather than
# at the next pin bump.
JBROWSE_COMPONENTS=${JBROWSE_COMPONENTS:-../jbrowse-components} \
  TUTORIALS_REF=HEAD Rscript scripts/tutorials.R --expand
