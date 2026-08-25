# Work order: unify the joint families — release 0.10.0

*Written 25 August 2026. Companion to `family-unification.md`, which records the same
decision taken at 0.8.0 for one family; this finishes the job for all four.*

**This ships as 0.10.0, not 1.0.0.** While the major version is 0 the API is not held out as
stable, so a breaking change is in order here. 1.0.0 comes later — the release that declares
the surface settled and carries the archived DOI a manuscript can cite — and may contain no
code change at all beyond what lands here.

This describes what the release needs to achieve and why. **The route is yours to work
out** — read the code, and if you find a better way to reach these outcomes than the one
sketched here, take it and say so. Where this note asserts something about the code,
verify it rather than trusting it; some of it was read quickly.

---

## 1. What is wrong now

The package ships four custom families:

| | dpars | `vars` |
|---|---|---|
| `bipois()` | 3 | `vint1` |
| `bipois_cens()` | 3 | `vint1`, `vint2` |
| `binegbin()` | 5 (one shared `shapex`) | `vint1` |
| `binegbin_cens()` | 6 (`shapexone`/`shapextwo`) | `vint1`, `vint2` |

Three problems, and they are one problem.

**The dispersion asymmetry is backwards.** `binegbin_cens()` got two source-specific
dispersions at 0.8.0; `binegbin()` never did. But in `binegbin()` every row is matched, so
both excess dispersions are informed by every row — it is the family where they are *easiest*
to identify. In `binegbin_cens()` the unpaired rows inform only `shapextwo`, so `shapexone`
rests on matched rows alone. The capability sits with the family that can least support it.

**`_cens` is the wrong word, and wrong inside brms specifically.** Censoring means the value
is known to lie in a set. On these rows the first count is not observed at all, and the
likelihood marginalises over its whole support. brms already uses `cens()` as an addition
term meaning exactly the bounded thing — `left`, `right`, `interval`, no code for "not
observed" — so a brms user meets a familiar word attached to an incompatible mechanism.

**The two families differ by almost nothing.** `binegbin_cens_lpmf`'s `y1_obs == 1` branch
*is* `binegbin_lpmf` — the roxygen says so and a test pins it. Once the dispersions unify,
the only real difference is whether a second `vint` is present. `family-unification.md` §2
rejected shipping "two families that differ by a single equality constraint, with duplicated
Stan code, duplicated `log_lik_*`/`posterior_predict_*` methods, a doubled test surface, and
a MethodsX paper obliged to explain why the distinction exists". That is what is still here.

---

## 2. What 0.10.0 should look like

**Two `custom_family` names — `bipois` and `binegbin` — with two constructors each:**

```
bipois()     bipois_partialobs()        binegbin()   binegbin_partialobs()
```

Both constructors in a pair return the **same** `name`, so there is one set of
post-processing methods per component distribution and one lpmf. They differ in `vars`:
one `vint` or two. The `_partialobs` variants take the observation flag; the plain ones
don't, and a user fitting fully paired data never has to know the flag exists. That last
point is the reason for two constructors rather than one mandatory flag — the simple case
should not pay for the complex one in comprehension, not just in keystrokes.

Both negative-binomial constructors carry the full six dpars.

**Why `partialobs`.** What is partially observed is the *pair*, which is a property of the
user's data. `unobs` and `unmatched` name only the minority row class; `partial` alone
collides with partial likelihood and partial pooling; `cens` is taken and means something
else.

**One likelihood implementation**, so the matched branch is not a second copy that can
drift. Two ways to get there, and **the second is worth testing first because it is
simpler if it works**.

*Delegation.* A single lpmf carrying the flag, with the one-`vint` signature delegating to
it with the flag set to 1. **This requires Stan function overloading** — brms derives the
Stan function name from the family name, so both constructors generate calls to
`binegbin_lpmf` with different arities, and that is two user-defined functions sharing one
name. Stan allows it only from **2.29** (February 2022). Current rstan and cmdstanr are
past that, but rstan sat at 2.21 on CRAN until 2023 and `DESCRIPTION` sets no floor on
either. If you take this route, add version floors and say why in `NEWS.md`, because the
failure mode otherwise is an opaque compile error.

*A literal in `vars`.* `vars` entries are pasted into the generated Stan call, so the
matched-only constructor may be able to pass the flag as a constant —
`vars = c("vint1[n]", "1")` — giving one Stan function, one signature, no overloading and
no version floor. Whether brms tolerates a non-variable entry there is a short test.
**Run it before committing to delegation.** If it works, prefer it and note in `NEWS.md`
why the design looks the way it does.

**Methods that work out which shape they have.** They already read `prep$data$vint1` and
`prep$data$vint2` directly, and `.get_dpar_any()` already resolves dpar spellings across
versions, so the machinery mostly exists.

---

## 3. Things you need that are not in the code

**Remove `_cens` and `_joint` completely — constructors and methods alike.** Nobody outside
this account has used the package, and the one consumer with stored fits handles its own
compatibility (below). 0.10.0 should contain no deprecated code, no forwarders and no
obligation to remove anything later.

`binegbin_joint`'s own documentation promises removal "in the next major version", which
0.10.0 is not. Remove it anyway, and **amend that sentence rather than leaving a documented
promise quietly unmet**.

**One consumer has stored fits that depend on the old names, and it is handling that
itself.** In the sibling project `tnc001-belize-em`, every fit in `fits/final` has
`family$name == "binegbin_cens"`. brms resolves post-processing late, by `get()` on that
stored string, so those fits need the 0.9.1 methods on the search path. That project pins
bicountbrms 0.9.1 in its own `renv.lock` until its analysis is re-run — which is a stronger
guarantee than anything this package could offer, because it freezes the whole package
rather than a hand-picked subset of functions. **Nothing is required of 0.10.0 on their
behalf.** Do not add a compatibility layer for them.

**The test suite has a demonstrated blind spot, in this exact code.** `NEWS.md` for 0.9.1
records that exchanging `shapexone` and `shapextwo` inside
`posterior_predict_binegbin_cens()` "passed 412 assertions across six test files with no
failures", and `binegbin.R`'s roxygen warns that the dpar order and the Stan signature can
be reordered independently, "silently swapping which rate or dispersion governs which
component". Treat the existing green suite as weak evidence for anything touching dispersion
routing. `tests/testthat/test-cens-predict.R` is the model for a test that would have caught
it — write the equivalent for the unified family **before** changing the family, and satisfy
yourself it fails on today's code for the right reason.

**Most users will take the one-`vint` path, and nothing currently tests it.** Whatever shape
the method branching takes, the untested path is the common one.

**Stored plain `binegbin` fits collide with the new name, and this is the one compatibility
question 0.10.0 cannot delegate.** The `_cens` fits belong to a project that pins its own
version. But stored *five-dpar* `binegbin` fits carry the name this release keeps, so they
land on the new unified methods whatever anyone pins, carrying a prep with a single
`shapex` and no `vint2`. That
should work — `.get_dpar_any()` falls back to `shapex` for both dispersion names, and the
absent `vint2` selects the matched branch — but it works because two independent mechanisms
happen to line up, which is the kind of thing that holds until it doesn't. **Test it
explicitly rather than assuming it**, and check first whether any such fits exist; the
Belize `fits/final` set is entirely `binegbin_cens`, but other directories were not looked
at.

**Do not regenerate the Belize fits, and do not try to keep them working.** They are being
written up, and that project pins 0.9.1 for as long as it needs to. Neither is this
release's problem.

---

## 4. What "done" looks like

- Four constructors, two family names, one lpmf and one set of post-processing methods per
  component distribution.
- A fully paired model fits without the user supplying or knowing about an observation flag.
- A partially paired model fits with `vint(y2, y1_obs)` and gives the same answers it does
  today.
- One-`vint` and two-`vint` results agree: for each of `log_lik`, `posterior_predict` and
  `posterior_epred`, the one-`vint` path on matched data equals the two-`vint` path on the
  same data with the flag set to 1.
- A test that fails if the two source-specific dispersions are exchanged anywhere they are
  routed.
- No deprecated code ships. `_cens` and `_joint` are gone entirely, and nothing is left
  behind to remove at some later version.
- `R CMD check` clean, and the Stan-fitting tests demonstrably *ran* rather than skipping —
  `helper-stan.R` notes that r-lib's action sets `NOT_CRAN=true`, so a run where every
  fitting test skipped for want of a toolchain looks identical to one where they passed.
- Every description of the family surface, anywhere in the repository, matches what ships
  (§5).

---

## 5. Documentation

Two phases. **Phase A makes everything true; phase B makes it findable.** Do A first and
verify it, so that if you run out of road the package is at least coherent. But do B in the
same pass rather than later — A already requires reading every word of the README, the
vignette and the articles, and B mostly moves that same text around. Splitting the two jobs
in time means reading ~4,900 words twice.

### Phase A — nothing in the repository describes a surface that no longer exists

**Everything that describes the four families goes stale at once**, and in this repository
that is more places than it looks. An inventory from a quick pass — **treat it as a starting
point, not a complete list**, and search for others:

| Where | What is in it |
|---|---|
| `DESCRIPTION` | The `Description` field names all four families in prose — "Four families are supplied… and a censoring-aware extension of each (bipois_cens, binegbin_cens)" — and asserts "The Negative-Binomial families let the two source-specific components differ in overdispersion". Both statements change |
| Roxygen in `R/*.R` | Family docs, `@details`, and worked `brm()` calls. The `binegbin()` example uses `shapex ~ 1`, a dpar that stops existing |
| Non-roxygen file headers in `R/*.R` | Long `#` comment blocks explaining the construction. These do not regenerate and are easy to miss |
| `man/*.Rd` | Regenerated — but only if you remember to roxygenise |
| `_pkgdown.yml` | Reference sections by family name; the section descriptions; the "Deprecated" section |
| `README.md` | ~4,900 words, including the five-vs-six dpar distinction and the notation bridge |
| `vignettes/bicountbrms.Rmd` | **Precompiled.** See the trap below |
| `vignettes/articles/paired-count-anatomy.Rmd` | Same, and it discusses parameters throughout |
| `NEWS.md` | The breaking change, and what someone holding fits made with `_cens` or `_joint` should do — install the last version that defined them, since brms resolves post-processing off the attached search path |
| `CLAUDE.md` | Repository conventions, which reference the family surface |
| `tests/` | Test names and comments referring to the old families |

**Correcting what is stale is only half of it — the `_partialobs` constructors are new and
need documentation written for them.** Someone arriving at the package should be able to
learn, without reading the source:

- **what the flag is and how to build it** — a 0/1 integer column, 1 where both counts were
  recorded, and where it goes in the formula;
- **what happens to each kind of row** — a matched row uses the full joint; a row whose
  first count was never recorded contributes the second count's marginal from the same
  model, so it still informs the shared component, the second source's rate and dispersion,
  and any group-level effects. It is not dropped, and it is not given a different model;
- **what you get afterwards** — the fitted model can impute the unobserved first count
  conditional on the observed second one, which is usually why someone wanted this;
- **the design consequence worth knowing before they collect data** — the first source's
  rate and excess dispersion are identified by the *matched* rows alone, so a design with
  few of those learns them weakly and leans on their priors. Someone with 20 matched rows in
  500 should meet that fact in the documentation rather than discover it in a posterior;
- **how to tell which shape a stored fit used**, since `family$name` is now the same either
  way and only the presence of `vint2` distinguishes them.

This is the package's distinguishing capability. It should be prominent rather than a note
at the bottom of a reference page.

**Three traps.**

*The vignettes are precompiled.* `bicountbrms.Rmd` is generated from `bicountbrms.Rmd.orig`
by `vignettes/precompile.R`. Edits made to the `.Rmd` alone are silently lost at the next
precompile, and a stale `.Rmd` alongside an updated `.Rmd.orig` is worse than either.

*`migration/family-unification.md` is a historical record, not documentation.* It describes
a decision taken at 0.8.0 and should keep saying what it says. Supersede it with a pointer
rather than rewriting it — the same way it handles the 0.9.0 rename that post-dates it.

*The word matters more than the mechanism here.* Wherever the docs currently say censoring,
what is true is that the first count was not observed at all and the likelihood marginalises
over its whole support. Somewhere prominent it is worth saying that this is unrelated to
brms's own `cens()` addition term, which means bounded values — a brms user will otherwise
map one onto the other.

### Phase B — the documentation is findable, and demonstrates the thing it is about

**The README is doing an article's job.** At ~4,900 words it sets out the construction, the
moments, the reparameterisation, prior guidance, validation and migration history. A README
that long is not read. Target 800–1,000 words: what the package is, how to install it, the
families in a table, the notation bridge from symbol to dpar name, one minimal example,
links. Everything else moves into articles — the sketch already agreed is three of them,
covering the families and their parameters, choosing priors, and migration and errata.

**`_pkgdown.yml` has no `articles:` section at all.** It carries `url`, `template` and four
`reference:` sections and nothing else, so the site currently has a reference index and no
articles — neither the get-started vignette nor `paired-count-anatomy` is deliberately
placed. Whatever comes out of the README needs a home on the site, or the split just hides
it somewhere else.

**No worked partially observed fit exists anywhere in the package.** Every fitted example in
every document uses fully matched data. That branch is the package's distinguishing
capability and the entire basis of the design it was written for, and it is never simulated,
never fitted, and its predictions never shown. It is also the first thing a reader who
installs the package will try. Write one: simulate a partially paired dataset, fit it,
check it, and show the imputed first margin against what was withheld. This is the reason
phase B belongs here rather than later — the example must demonstrate the new surface, so
it cannot be written before the API lands, and writing it is the best test that the new
surface is usable.

There is a second-order benefit to phase B worth knowing. Several passages currently run
word-for-word between the public README and a manuscript describing this method — runs of
30 to 45 consecutive identical words. Moving that material into articles relocates the
overlap to a place where overlap with a paper is expected, and any rewriting it prompts
reduces it further.

---

## 6. Deliberately not in this work order

Release mechanics — version bump, `.zenodo.json`, tagging, the Zenodo webhook, the DOI —
are not yours. Nor is 1.0.0: that is a later, separate decision about when the surface has
proven itself, tracked in `tnc001-belize-em/docs/bicountbrms-release-checklist.md`.

If something here turns out to be wrong about the code, or a decision looks mistaken once
you are inside it, say so rather than working around it.

---

## 7. Outcome, recorded 25 August 2026

Added after the release was built. §6 asks for anything that turned out wrong
about the code to be said rather than worked around; this is that.

**The literal in `vars` works, so delegation was not needed.** §2 asked for it
to be tested before committing to Stan function overloading. It was, and it
holds: `brms:::stan_log_lik_custom()` is the only consumer of `family$vars`,
`custom_family()` validates the entries no further than `as.character()`, and a
non-variable entry is pasted into the generated call verbatim. `binegbin()`
declares `vars = c("vint1[n]", "1")` and reaches the same `binegbin_lpmf` with
the flag fixed at `1`; `standata()` carries no `vint2`. **No floor on the Stan
version was added, and none is needed.** Because that pasting is not a
documented brms guarantee, `tests/testthat/test-stancode-shape.R` pins the
generated call for all four constructors and runs in the fast suite.

**Stored plain `binegbin` fits exist, and they are older than §3 assumes.** §3
asks whether any do, and notes that `fits/final` is entirely `binegbin_cens` —
correct, all ten of them. But elsewhere under `tnc001-belize-em/fits/` there
are **103** fits with `family$name == "binegbin"`, and every one carries
pre-0.7.0 rate names (`lambdaem`/`lambdalb`) *as well as* the single `shapex`.
So the compatibility path is three independent mechanisms lining up, not the
two §3 describes: `.get_rate()` on the rates, `.SHAPEX*_NAMES` falling through
to `shapex`, and an absent `vint2` selecting the matched branch. All three are
now pinned in `test-dpar-compat.R`, and one real fit
(`fits/fit_alb_binegbin.rds`) was post-processed under 0.10.0 as a check:
`log_lik()`, `posterior_epred()`, `posterior_predict()` and `loo()` all return
finite, sensible values. Nothing was regenerated.

**One thing §5 does not mention needed a code change.**
`binegbin_mfd_to_dpars()` emitted a dpar named `shapex` when given `kappax`.
After unification no family accepts that name, so the converter's output could
no longer be passed to `brm()` — a functional defect, not a stale sentence.
`kappax` now writes `shapexone` and `shapextwo` at a common value. The inverse
still *accepts* `shapex`, because its input is a stored fit and that is the
genuine dpar name on all 103 of them.

**The Stan-side blind spot §3 warns about was still open, and is now closed.**
Unification made the Stan-vs-Stan equivalence checks vacuous — they compared a
function to itself. Deleting them would have left the Stan lpmf's dispersion
routing tested only by a grid that passes the *same* value for both
dispersions, which a swap cannot disturb: injecting one into
`binegbin_stan_funs` leaves that grid byte-identical. They were replaced by
asymmetric routing checks, which catch the same swap by 45 log units. All five
mutation sites — the three R methods, the R reference and the Stan lpmf — are
now caught.

The `_partialobs` naming, the removal of every `_cens`/`_joint` name with no
deprecation layer, and the decision not to regenerate the Belize fits were all
adopted as written.
