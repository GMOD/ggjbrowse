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
  table. Each entry holds the `jb2export` arguments and the ggjbrowse recipe,
  so adding an example is adding one entry. It needs a built jbrowse-components checkout
  (`JBROWSE_COMPONENTS`), the installed package and ImageMagick. `@jbrowse/img`
  4.3.0 is too old: use 5.0.0-beta or the checkout, which the two clustered
  rows need until a release carries commit `e2cef5c883`.

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
- `scripts/take-release.sh [version]` takes a published release from npm
  instead, the latest by default.
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
  [docs/how-does-this-work.md](docs/how-does-this-work.md#who-decides-what).
- [agent-docs/engine-owns-display-decisions.md](agent-docs/engine-owns-display-decisions.md)
  lists the rules R still copies.
