# Unifying the symmetric and asymmetric joint families

**Written for pairedcountbrms 0.8.0; the package is now bicountbrms 0.9.0.**
This note records why `binegbin_joint()` was generalised to six distributional
parameters rather than joined by a second family, what was measured to establish
that the change is safe, and what the Belize EM project (`tnc001-belize-em`)
must do to adopt it.

The generalisation itself is unchanged by the 0.9.0 rename and split, and the
measurements below stand as recorded — they were made under the name
`pairedcountbrms`, and are left saying so. Only the forward-looking adoption
steps in §5 have been updated to the current package name; install
`bicountbrms` where they say to install the package.

---

## 1. What changed

Up to 0.7.0, `binegbin_joint()` declared five dpars:

```
mu, lambdaone, lambdatwo, shapes, shapex
```

with the single `shapex` governing **both** source-only excess components. In
parallel, the Belize project carried a local family, `binegbin_joint_ax()`
(`R/binegbin_joint_ax.R`), identical in every respect except that the excess
dispersion was split by margin — six dpars, and the family every shipping fit
in that project uses.

0.8.0 makes that split the package's own. `binegbin_joint()` now declares:

```
mu, lambdaone, lambdatwo, shapes, shapexone, shapextwo
```

all log-linked, `lb = 0`, with the response still carrying
`vint(y2, y1_obs)`. The symmetric model is the constraint
`shapexone == shapextwo`, expressed as a formula rather than as a separate
family.

## 2. Design: two families (A) or one generalised family (B)

**Option A** — keep `binegbin_joint()` at five dpars and add
`binegbin_joint_ax()` alongside it. No breaking change; but the package would
then ship two families that differ by a single equality constraint, with
duplicated Stan code, duplicated `log_lik_*`/`posterior_predict_*` methods, a
doubled test surface, and a MethodsX paper obliged to explain why the
distinction exists.

**Option B** — generalise to six dpars. One family, one Stan function, one
likelihood to describe. The cost is a breaking change to the formula API.

**B was chosen.** The argument that decided it is not aesthetic; it is that the
cost of B turned out to be much smaller than it looks, and that was measured
rather than assumed. The decisive question was what happens to **existing
`brmsfit` objects**, and it has two independent parts.

### 2a. Does a stored fit still find the right post-processing function?

Yes, and it never could have failed. brms resolves a custom family's methods
through `brms:::custom_family_method()`:

```r
function (family, name) {
    out <- family[[name]]
    if (!is.function(out)) {
        out <- paste0(name, "_", family$name)
        out <- get(out, family$env)
    }
    out
}
```

The lookup is `get(name, family$env)`. Inspecting a stored fit shows
`family$env` is an **empty** environment whose parent is `R_GlobalEnv`, so the
`get()` walks straight out to the live search path. Resolution is therefore
**late-bound**: a stored fit runs whatever code is attached *now*, not
anything frozen at fit time.

Measured directly (scratch script `02-late-binding.R`):

| Test | Result |
|---|---|
| `log_lik()` on a stored fit with the family file **not** sourced | errors: `object 'log_lik_binegbin_joint_ax' not found` |
| same, after sourcing the family | returns a 3 × 516 matrix of finite values |
| redefine the function to return `-999`, call `log_lik()` again | returns `-999` |

The last row is the proof: the fit picks up the redefinition. So changing the
package cannot break function *resolution* for a fit whose family **name** is
unchanged — and under B the name `binegbin_joint` is unchanged.

### 2b. Does the function body still find the dpars it wants?

This is where B could genuinely have broken, and it is the only real risk. A
`brmsfit` stores its own family object, so `prepare_predictions()` builds
`prep$dpars` from *that* family. A five-dpar fit arrives with `shapex`
populated and `shapexone`/`shapextwo` absent; a six-dpar body asking for
`shapexone` would fail on a bare `!is.null(x) is not TRUE` from
`brms::get_dpar()`.

The package already had the idiom for this, introduced at 0.7.0 for the
`lambdaem` → `lambdaone` rename: resolve the name **against the fit** and hand
that name to brms (`.get_rate()` in `R/utils.R`). 0.8.0 generalises it to
`.get_dpar_any()`, an ordered candidate list:

```
shapexone  <-  shapexone | shapexem | shapex
shapextwo  <-  shapextwo | shapexlb | shapex
```

Falling back to `shapex` **last** is what makes the symmetric case work with no
shim: a five-dpar fit resolves both excess dispersions to its single `shapex`,
which is exactly the constraint it was estimated under.

Measured on seven stored five-dpar fits (`05-symmetric-nonregression.R`),
comparing the 0.8.0 code against the 0.7.0 code taken from git HEAD, on the
same `prep`:

| Fit | dpars stored | obs × draws | max abs diff | bitwise identical | `posterior_predict` identical |
|---|---|---|---|---|---|
| `mfd3/fit_alb_p2_vMc` | 5 | 516 × 3000 | 0 | yes | yes |
| `mfd3/fit_yft_p2_vMc` | 5 | 516 × 3000 | 0 | yes | yes |
| `mfd3/fit_bsh_p2_vMc_sym` | 5 | 516 × 3000 | 0 | yes | yes |
| `mfd3/cv_sai_conRE` | 5 | 516 × 3000 | 0 | yes | yes |
| `onemodel/S_alb` | 5 | 516 × 3000 | 0 | yes | yes |
| `onemodel/U_yft` | 5 | 516 × 3000 | 0 | yes | yes |
| `onemodel/Anp3_adopted` | 5 | 1548 × 3000 | 0 | yes | yes |

and end-to-end through brms's own dispatch, with only the 0.8.0 code visible:

```
log_lik():           OK, dim 20x516, all finite
posterior_predict(): OK, dim 20x516, all >= 0
loo():               OK, elpd_loo = -1747.247
```

**Conclusion: stored fits survive, and no compatibility shim is required.**
Retaining the old family name as a thin wrapper — the shim contemplated in
advance — is unnecessary, because under B the family name does not change and
the dpar fallback covers the only thing that does. Had B instead *renamed* the
family (say to `binegbin_joint_asym`), a shim would have been needed, and it
would have had to be a real family constructor named `binegbin_joint()`
returning the old five-dpar object plus `log_lik_binegbin_joint()` /
`posterior_predict_binegbin_joint()` methods — i.e. keeping the whole old
family alive, which is option A wearing a different hat. That is a further
argument for not renaming.

## 3. The equivalence check (run once, 2026-08-03)

**This is a record, not a test.** It was run once during the migration and is
deliberately *not* part of `tests/testthat/`. Two reasons: it depends on fits
stored outside this repository, so as a test it could only ever skip on CI, on
CRAN and on any other checkout — and a permanently-skipped test reads as
coverage that is not there; and its subject is frozen, since those ten fits are
fixed artefacts, so once the comparison passes it can only change if *this
package* changes. The package-side property is pinned by portable tests that
need no external repository (see the list at the end of this section).

The authoritative confirmation will come from `tnc001-belize-em` itself, which
refits its models against the packaged family after adoption.

Environment for the run below: pairedcountbrms 0.8.0 (working tree), brms
2.23.0, R 4.6.0 (2026-04-24 ucrt), Windows 11 x64.

The packaged six-dpar family must be **likelihood-identical** to
`binegbin_joint_ax()`, or adopting it silently changes ten shipping fits.

Method: read each stored fit, build one `prep` from it, and evaluate the
pointwise log-likelihood twice — once through the project-local
`log_lik_binegbin_joint_ax()` and once through the packaged
`log_lik_binegbin_joint()`, which resolves that fit's `lambdaem`/`lambdalb`/
`shapexem`/`shapexlb` spellings via `.get_dpar_any()`. Both functions see the
same numbers, so any disagreement is a real difference in the likelihood and
not a naming artefact. Compared **element by element** over the full
`draws × observations` array — not `elpd`, which can agree while individual
points do not.

| Fit | obs | draws | elements | matched rows | y2-only rows | max abs diff | bitwise identical |
|---|---|---|---|---|---|---|---|
| `Lax_alb_vMc`  | 516 | 3000 | 1,548,000 | 98 | 418 | 0 | 100% |
| `Lax_bet_vMc`  | 516 | 3000 | 1,548,000 | 98 | 418 | 0 | 100% |
| `Lax_bil_vMc`  | 516 | 3000 | 1,548,000 | 98 | 418 | 0 | 100% |
| `Lax_bsh_vMc`  | 516 | 3000 | 1,548,000 | 98 | 418 | 0 | 100% |
| `Lax_eskh_vMc` | 516 | 3000 | 1,548,000 | 98 | 418 | 0 | 100% |
| `Lax_sai_vMc`  | 516 | 3000 | 1,548,000 | 98 | 418 | 0 | 100% |
| `Lax_swo_vMc`  | 516 | 3000 | 1,548,000 | 98 | 418 | 0 | 100% |
| `Lax_ttus_vMc` | 516 | 3000 | 1,548,000 | 98 | 418 | 0 | 100% |
| `Lax_wah_vMc`  | 516 | 3000 | 1,548,000 | 98 | 418 | 0 | 100% |
| `Lax_yft_vMc`  | 516 | 3000 | 1,548,000 | 98 | 418 | 0 | 100% |

**15,480,000 values, every one bitwise identical** (max absolute difference
exactly `0`, not merely within tolerance). All values finite.

To reproduce: read each fit, build one `prep`, and evaluate
`log_lik_binegbin_joint(i, prep)` against `log_lik_binegbin_joint_ax(i, prep)`
for every observation, comparing the full draws vector each time. The whole
check is about thirty lines and takes roughly six minutes over the ten fits.

Separately, and **this is** in the suite
(`tests/testthat/test-binegbin_joint_asym.R`), these self-contained checks pin
the same properties without any external dependency:

- the asymmetric likelihood is a proper pmf on both branches, and the marginal
  identity (matched branch integrated over `y1` == the `y2`-only branch) still
  holds with `shapexone != shapextwo`;
- `shapexone` does **not** enter the `y2`-only branch, and does enter the
  matched one (so the first assertion is not vacuous);
- transposing both rates and both dispersions transposes the joint, which pins
  `shapexone` to the `y1` term and `shapextwo` to the `y2` term against a
  silent swap in the Stan signature;
- `shapexone == shapextwo` reproduces `binegbin_lpmf_r` exactly over a
  108-point parameter grid (symmetric non-regression);
- a five-dpar `prep` and a `shapexem`/`shapexlb` `prep` both resolve to the
  packaged names and give identical results;
- the Stan lpmf against the R reference on an asymmetric grid (`< 1e-8`);
- **simulation recovery** with `shapexone = 0.7` and `shapextwo = 6.0` — an
  order of magnitude apart — asserting not only that each truth is bracketed
  but that the 95% interval for `log(shapexone) - log(shapextwo)` excludes
  zero, i.e. that the data distinguished the two rather than the priors doing
  it.

## 4. The dpar mapping table

The package is deliberately **source-agnostic**; the project is
domain-specific. The package keeps the generic names.

| Belize project (`binegbin_joint_ax`) | Package 0.8.0 (`binegbin_joint`) | Role |
|---|---|---|
| `mu`         | `mu`        | rate of the shared component |
| `lambdaem`   | `lambdaone` | rate of the source-1-only excess (EM) |
| `lambdalb`   | `lambdatwo` | rate of the source-2-only excess (LB) |
| `shapes`     | `shapes`    | NB2 dispersion of the shared component |
| `shapexem`   | `shapexone` | NB2 dispersion of the source-1-only excess |
| `shapexlb`   | `shapextwo` | NB2 dispersion of the source-2-only excess |

and the two supplementary integers:

| Belize project | Package convention | Binding | Role |
|---|---|---|---|
| `vint(y_lb, ...)`  | `vint(y2, ...)`     | `vint1[n]` | the always-observed second count |
| `vint(..., em_obs)` | `vint(..., y1_obs)` | `vint2[n]` | 0/1 flag: was the first count observed? |

**The `vint()` terms need no renaming.** They are matched *positionally* —
`vars = c("vint1[n]", "vint2[n]")` — so the data column names are free. `y_lb`
and `em_obs` may stay exactly as they are; only the documentation convention
differs. The right-hand column is what the package's docs call them, not a
requirement.

### 4a. A correction to the expected migration cost

The renaming was expected to break read-outs that grep posterior draw columns —
`b_lambdaem_Intercept` becoming `b_lambdaone_Intercept`. **It does not, because
no such column exists.** Every dpar in the project's models is defined through
`nlf()` in terms of non-linear parameters, and it is the *nlpar* names that
reach the draws. `Lax_alb_vMc`'s formula:

```
y_em | vint(y_lb, em_obs) ~ 1
mu       ~ eta + log_inv_logit(con)
lambdaem ~ log(2) + eta + log_inv_logit(-con) + log_inv_logit( 2 * methd)
lambdalb ~ log(2) + eta + log_inv_logit(-con) + log_inv_logit(-2 * methd)
shapes   ~ log(1/kappas^2)
shapexem ~ log(1/kappaxem^2)
shapexlb ~ log(1/kappaxlb^2)
eta ~ 1 + (1 | vessel) + (1 | vessel:trip_id)
con ~ 1 + (1 | vessel)
methd ~ 1
kappas ~ 1
kappaxem ~ 1
kappaxlb ~ 1
```

and its entire non-`r_`/`z_`/`L_` variable set:

```
b_eta_Intercept    b_con_Intercept    b_methd_Intercept
b_kappas_Intercept b_kappaxem_Intercept b_kappaxlb_Intercept
sd_vessel__eta_Intercept  sd_vessel:trip_id__eta_Intercept  sd_vessel__con_Intercept
```

Scanning all 54 stored `binegbin_joint*` fits across `fits/fold`, `fits/mfd3`
and `fits/onemodel`: **zero** have a dpar name in any draw column.

So the dpar rename touches **model source only** — the left-hand sides of
`nlf()` terms — and touches **no stored draws, no read-out, and no saved
object**. The nlpar names `kappaxem`/`kappaxlb` are the project's own choice
and the package does not require renaming them; if the project prefers to keep
its EM/LB vocabulary in the parameters it actually reads, it can, and the
mapping above stays confined to six `nlf()` left-hand sides.

## 5. Adoption steps for `tnc001-belize-em`

Not performed here — this is the other session's work. In order:

1. **Delete nothing yet.** Keep `R/binegbin_joint_ax.R` on disk until step 5.

2. **Install bicountbrms** (0.9.0 or later) and attach it where the local family
   was sourced. The family, its dpars and its likelihood are unchanged from the
   `pairedcountbrms` 0.8.0 under which the checks below were run; only the
   package name differs.

3. **Rewrite the `nlf()` left-hand sides** in the model-building code:

   ```r
   nlf(lambdaem ~ ...)   ->  nlf(lambdaone ~ ...)
   nlf(lambdalb ~ ...)   ->  nlf(lambdatwo ~ ...)
   nlf(shapexem ~ ...)   ->  nlf(shapexone ~ ...)
   nlf(shapexlb ~ ...)   ->  nlf(shapextwo ~ ...)
   ```

   Measured over the Belize repo (`*.R` and `*.qmd`), this is **142 `nlf()`
   sites across 37 files**: 55 `lambdaem`, 54 `lambdalb`, 17 `shapexem`, 16
   `shapexlb`. They sit on **128 lines** — some lines carry two `nlf()` terms —
   so a line-counting `grep -c` reports 128 and an occurrence-counting
   `grep -oh ... | wc -l` reports 142. Both are right; a substitution pass
   should change 142 things. The pattern is

   ```bash
   grep -rEoh "nlf\( *(lambdaem|lambdalb|shapexem|shapexlb) *~" --include=*.R --include=*.qmd .
   ```

   Re-derive the counts before and after as a check.

   The right-hand sides, the nlpar names (`eta`, `con`, `methd`, `kappas`,
   `kappaxem`, `kappaxlb`), and every `prior()` call are **unchanged** — the
   priors attach to nlpars, not to dpars. `R/mfd3/mfd3_priors.R` mentions
   `binegbin_joint_ax` only in a comment and contains no `dpar =` argument, so
   it needs no edit. Nothing in this note asks for a prior to change.

4. **Swap the family and stanvars:**

   ```r
   family   = binegbin_joint_ax(),            ->  family   = binegbin_joint(),
   stanvars = binegbin_joint_ax_stanvars()    ->  stanvars = binegbin_joint_stanvars()
   ```

5. **Decide what to do about the existing fits.** They are stored under the
   family name `binegbin_joint_ax`, so brms will look up
   `log_lik_binegbin_joint_ax` / `posterior_predict_binegbin_joint_ax` on the
   search path. Two options:

   - **Keep the file.** `R/binegbin_joint_ax.R` continues to work untouched and
     nothing needs refitting. Simplest, and correct: the equivalence test in
     §3 proves the two give bitwise-identical likelihoods.
   - **Replace it with a three-line alias**, if the intent is that only the
     packaged likelihood is ever evaluated:

     ```r
     log_lik_binegbin_joint_ax <-
       function(i, prep) bicountbrms::log_lik_binegbin_joint(i, prep)
     posterior_predict_binegbin_joint_ax <-
       function(i, prep, ...) bicountbrms::posterior_predict_binegbin_joint(i, prep, ...)
     ```

     These work because the packaged functions resolve `lambdaem`/`lambdalb`/
     `shapexem`/`shapexlb` directly (§2b). Note the *aliases still need the
     `_ax` names* — that requirement comes from the stored fits' family name,
     not from the package.

   Refitting is not required either way. Fits produced *after* step 4 will
   carry the family name `binegbin_joint` and need neither.

6. **Verify there, not here.** The models are being refitted against the
   packaged family anyway, so the authoritative equivalence check belongs in
   `tnc001-belize-em` at that point — which is why it is not a test in this
   package (§3).

   When comparing a refit against its stored predecessor, note that the two
   are **not** expected to agree draw for draw: refitting re-runs the sampler,
   so draws differ by Monte Carlo error even under an identical likelihood.
   Compare the posteriors (e.g. overlapping intervals for `b_eta_Intercept`,
   `b_con_Intercept`, `b_methd_Intercept` and the three kappas), or compare
   `log_lik()` evaluated at a *common* parameter vector, which is the exact
   comparison §3 made and is the one that can be identical.

## 6. Things that will bite

Ordered by how easily each is misdiagnosed as a packaging fault when it is not.

**Sampling diagnostics on the asymmetric model.** Freeing the second dispersion
adds a weakly-identified parameter, and the likelihood for a Negative-Binomial
shape is flat towards large values — a big shape is nearly Poisson, so the data
cannot separate 30 from 300. brms leaves a custom family's non-`mu` dpars flat
and improper. Run unregularised, this package's own asymmetric recovery test
produced **680 divergent transitions and a max Rhat of 1.54**; with weakly
informative priors on the dispersions it is clean. If adoption produces
divergences, look here first — it is a property of the model, not of the
packaging. The Belize models already put priors on the dispersions through the
`kappa` parameterisation, so they should not hit this, but the failure
signature is worth recognising.

**`binegbin_mfd_to_dpars()` now takes `kappaxone`/`kappaxtwo`**, and
`binegbin_dpars_to_mfd()` takes `shapexone`/`shapextwo`, returning the
per-margin names. Use these rather than re-deriving `shape = 1/kappa^2` by
hand. `kappax` still works and still returns `shapex`; supplying it alongside
either new argument is an error.

**A refit will not match its stored predecessor draw for draw.** See step 6
above. Monte Carlo error is expected; the likelihood is what is identical.

**`vint()` needs no renaming.** `y_lb` and `em_obs` are matched positionally
and can stay exactly as they are. Only the four `nlf()` left-hand sides change.

**Priors need no editing.** They attach to nlpars (`eta`, `con`, `methd`,
`kappas`, `kappaxem`, `kappaxlb`), not to dpars. `R/mfd3/mfd3_priors.R`
contains no `dpar =` argument and mentions `binegbin_joint_ax` only in a
comment. Nothing in this migration requires a prior to change.

**Read-outs are safe.** No dpar name appears in any stored draw column — see
§4a. Scanned across all 54 stored `binegbin_joint*` fits: zero hits.

## 7. Identifiability note carried into the docs

`shapextwo` governs the always-observed margin and appears on **both**
likelihood branches, so the `y2`-only rows inform it. `shapexone` appears only
on the matched branch — it is integrated out of the `y2`-only branch along with
the `y1` margin — so it is identified **solely by the matched rows**, exactly as
`lambdaone` is.

This is not incidental for the Belize fits: they carry **98 matched rows against
418 `y2`-only rows**, so `shapexone` is learned from under a fifth of the data
while `shapextwo` uses all of it. Any comparison of the two dispersions should
be read against that asymmetry, and the prior on `shapexone` does
correspondingly more of the work. The same caveat now appears in the family's
roxygen block and in the README's Limitations section.
