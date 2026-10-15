# Contributing

## Develop

- Install the dependencies in `DESCRIPTION`, then `devtools::load_all()`.
- Run the tests with `devtools::test()`. They use the data in `inst/extdata`.
  `scripts/extdata-*.sh` make the human clips.
- Run `devtools::document()` after editing a roxygen comment. `man/` and
  `NAMESPACE` are generated.
- `scripts/figures.R` redraws the README's quick-start figure.
- `scripts/gallery.R` redraws the README's JBrowse-and-R gallery and writes
  `docs/gallery.md`, the configs in `gallery/jb2export` and the README's
  table. Each entry holds a caption, the `jb2export` arguments and the
  ggjbrowse recipe, so adding an example is adding one entry.
  - Drawing needs a built jbrowse-components checkout (`JBROWSE_COMPONENTS`)
    at the version [docs/gallery.md](docs/gallery.md) names, the installed
    package and ImageMagick.
  - `--docs` rewrites the table and `docs/gallery.md` without drawing.
  - `--r` draws only the R half, which needs no jbrowse-components checkout:
    `Rscript scripts/gallery.R --r inversion`.
  - `--check` runs each R recipe, then rewrites as `--docs` does. CI runs it
    and `scripts/check-readme.R`, which runs the R blocks of the README and the
    guides in `docs/`, and fails
    if either leaves the docs changed.
- `scripts/tutorials.R` builds `docs/tutorials/` from `tutorials/*.Rmd`.
  - A source is R Markdown whose `<!-- jbrowse ... -->` comments pull
    sections of the matching JBrowse website tutorial, at the
    jbrowse-components commit `tutorials/JBROWSE_COMPONENTS` names, with the
    options the script's header lists.
  - Each `drop`, `keep`, `only` or `replace` has to match the website text,
    so a website edit that removes it fails the build.
  - A directive with `assign` writes an R chunk from a website track config
    instead, so `display <- list(...)` is the track's own `displayDefaults`.
  - `--status` lists the website commits since the pin. Bump the pin, rerun,
    and fix whatever directive stops matching.
  - `--check`, which CI runs, knits each tutorial again from the data on
    jbrowse.org and fails if its markdown changed.

## Updating the engine

- `inst/js/jbrowse-tables.js` is `@jbrowse/tables` bundled for a bare JavaScript
  host, committed so installing needs no JavaScript toolchain.
- Its first line names the release, or the jbrowse-components commit.
- To rebuild it from a checkout beside this repo:

```console
cd ../jbrowse-components && pnpm install && cd -
scripts/build-js.sh
```

- `build-js.sh` reads `$JBROWSE_COMPONENTS`, default `../jbrowse-components`.
  It then runs `scripts/credit-bundle.mjs`, which reads the bundle's esbuild
  metafile and writes each bundled package's version, license and copyright
  to `inst/COPYRIGHTS`, which DESCRIPTION's `Copyright` field names.
- `scripts/take-release.sh [version]` takes a published release from npm
  instead, the latest by default. The npm release names no commit to credit
  from, so `inst/COPYRIGHTS` keeps the last `build-js.sh` list. Before a CRAN
  submission, rebuild with `build-js.sh` so one commit gives the bundle and
  its credits.
- `scripts/check-against-jbrowse.sh` rebuilds, installs and runs the tests in
  one step. jbrowse-components' downstream canary runs it against its main
  nightly.

## CI

- **R-CMD-check:** runs on push, pull request and by hand, on Linux, macOS and
  Windows with R release, and Linux with oldrel.
- **Bundle:** runs nightly. It takes each new `@jbrowse/tables` release from
  npm, runs R CMD check on it, then commits the bundle to `main` and starts
  R-CMD-check for every platform. A failed check leaves `main` alone and the
  next run tries again.

## Where decisions live

- R draws and the engine decides. See
  [docs/faq.md](docs/faq.md#who-decides-what).
- [agent-docs/engine-owns-display-decisions.md](agent-docs/engine-owns-display-decisions.md)
  lists the rules R still copies.
