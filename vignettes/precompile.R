# Precompile the vignette and the pkgdown articles.
#
# Two documents fit brms models: vignettes/bicountbrms.Rmd and the
# partially-paired-designs article. Building them on every R CMD check and CI run
# would need a Stan toolchain and several minutes of compilation, so each
# _<name>.Rmd source is knitted here, on a machine with Stan, and the resulting
# <name>.Rmd holds the output as static text. Downstream builds only render
# markdown.
#
# Because the output is static text, .github/workflows/pkgdown.yaml can use
# dependencies: '"hard"' and install no Stan backend.
#
# Run this by hand after changing either _<name>.Rmd, then commit BOTH files of
# each pair. The knitr::knit() calls execute every chunk, so the numbers in the
# committed .Rmd files are real output, not transcribed.
#
#   Rscript vignettes/precompile.R
#
# EDIT THE _<name>.Rmd, NEVER THE <name>.Rmd -- the latter is generated and the
# next run of this script will overwrite it.
#
# FIGURES. Downstream builds do not re-run the chunks, so any figure must
# already exist on disk. The two articles set dev = "svg" and
# fig.path = "figure/" in their setup chunks, so their plots are written to
# vignettes/articles/figure/ and must be committed with the .Rmd files. The
# vignette produces no figures; a plotting chunk added there needs the same
# setup.

stopifnot(
  "run from the package root" = file.exists("DESCRIPTION"),
  "needs a Stan backend" =
    requireNamespace("cmdstanr", quietly = TRUE) ||
    requireNamespace("rstan", quietly = TRUE)
)

# STALE INSTALL. A knit against an out-of-date installed bicountbrms does not
# fail: knitr captures the errors as chunk output and exits 0, so the article
# shows `#> Error` where the fitted numbers should be. The check below compares
# the installed version with DESCRIPTION. It cannot detect an install of the
# same version made before later source edits, so run R CMD INSTALL . before
# every precompile.
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
# The reparameterisation article fits no model, but it is precompiled with the
# others so that every _<name>.Rmd source is built the same way. The prediction
# vignette runs no code, is a plain .Rmd file with no _ source, and is not
# listed.
targets <- list(
  list(dir = "vignettes",          stem = "bicountbrms"),
  list(dir = "vignettes/articles", stem = "reparameterisation"),
  list(dir = "vignettes/articles", stem = "partially-paired-designs")
)

for (tg in targets) {
  message("Knitting ", tg$dir, "/_", tg$stem,
          ".Rmd -- a document that fits a model takes several minutes.")
  old <- setwd(tg$dir)
  knitr::knit(
    input  = paste0("_", tg$stem, ".Rmd"),
    output = paste0(tg$stem, ".Rmd")
  )
  setwd(old)
}

message("Done. Commit both the _<name>.Rmd and the <name>.Rmd of each pair.")
