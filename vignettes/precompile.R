# Precompile the vignette and the pkgdown articles that fit models.
#
# Both vignettes/bicountbrms.Rmd and
# vignettes/articles/paired-count-anatomy.Rmd fit real brms models. Building
# those on every R CMD check and every CI run would need a full Stan toolchain
# plus minutes of compilation, so they are precompiled instead: each .Rmd.orig
# source is knitted HERE, on a machine with Stan, and the resulting .Rmd carries
# the output as static text. Downstream builds only render markdown.
#
# This is what lets .github/workflows/pkgdown.yaml use dependencies: '"hard"'
# and install no Stan backend at all.
#
# Run this by hand after changing either .Rmd.orig, then commit BOTH files of
# each pair. The knitr::knit() calls execute every chunk, so the numbers in the
# committed .Rmd files are real output, not transcribed.
#
#   Rscript vignettes/precompile.R
#
# EDIT THE .Rmd.orig, NEVER THE .Rmd -- the .Rmd is generated and the next run
# of this script will overwrite it.
#
# FIGURES. Downstream builds do not re-run the chunks, so any figure must
# already exist on disk. paired-count-anatomy sets dev = "svg" and
# fig.path = "figure/" in its setup chunk, so its plots land in
# vignettes/articles/figure/ and must be COMMITTED alongside the .Rmd. The
# vignette currently produces no figures; if you add a plotting chunk there, do
# the same.

stopifnot(
  "run from the package root" = file.exists("DESCRIPTION"),
  "needs a Stan backend" =
    requireNamespace("cmdstanr", quietly = TRUE) ||
    requireNamespace("rstan", quietly = TRUE)
)

# STALE INSTALL. A knit against an out-of-date installed bicountbrms does not
# fail. knitr captures the errors as chunk output and exits 0, so the result is
# an article carrying `#> Error` where the fitted numbers should be. That
# happened at 0.10.0 against an installed 0.9.0, whose binegbin() still declared
# five dpars and whose binegbin_mfd_to_dpars() still emitted `shapex`; the knit
# reported success and the article was wrong. Compare the versions rather than
# trusting the exit code.
if (!requireNamespace("bicountbrms", quietly = TRUE)) {
  stop("bicountbrms is not installed. Run R CMD INSTALL . first.", call. = FALSE)
}
local({
  installed <- as.character(utils::packageVersion("bicountbrms"))
  declared  <- as.character(read.dcf("DESCRIPTION")[1, "Version"])
  if (!identical(installed, declared)) {
    stop(sprintf(
      "installed bicountbrms is %s but DESCRIPTION declares %s.\n  Run R CMD INSTALL . before precompiling.",
      installed, declared), call. = FALSE)
  }
})

# input/output are relative to the directory each pair lives in, so knitr
# resolves any relative paths inside the document the same way a build would.
# Every target here fits at least one model, which is why it is precompiled.
# The two articles that do not fit anything -- choosing-priors and
# migration-and-errata -- are plain .Rmd with no .orig, and are absent from
# this list on purpose.
targets <- list(
  list(dir = "vignettes",          stem = "bicountbrms"),
  list(dir = "vignettes/articles", stem = "paired-count-anatomy"),
  list(dir = "vignettes/articles", stem = "families-and-parameters"),
  list(dir = "vignettes/articles", stem = "partially-observed-fit")
)

for (tg in targets) {
  message("Knitting ", tg$dir, "/", tg$stem,
          ".Rmd.orig -- this fits a model, allow a few minutes.")
  old <- setwd(tg$dir)
  knitr::knit(
    input  = paste0(tg$stem, ".Rmd.orig"),
    output = paste0(tg$stem, ".Rmd")
  )
  setwd(old)
}

message("Done. Commit both the .Rmd.orig and the .Rmd of each pair.")
