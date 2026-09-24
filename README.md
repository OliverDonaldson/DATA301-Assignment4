# DATA 301 Assignment IV

Oliver Donaldson, 300659923.

**Exercise 1:** three estimators for the rate of an exponential distribution (maximum likelihood,
an IQR-based analogy estimator, and a first-order bootstrap bias-corrected version of the latter),
compared by Monte Carlo over `lambda` ∈ {0.1, 1, 10} × `n` ∈ {5, 30, 100}.
**Exercise 2:** exploratory visual analysis of quarterly Australian beer production
(`ausbeer`, `fpp2`).

## Reproducing this

```
Rscript install-deps.R
quarto render "DATA301-A4-Oliver-Donaldson-300659923.qmd"
```

Produces `DATA301-A4-Oliver-Donaldson-300659923.pdf`. Needs R, Quarto and a LaTeX installation
(`quarto install tinytex` if there is none). In RStudio, opening the `.qmd` and clicking Render
does the same.

Every random quantity is generated under a seed fixed in the setup chunk, so a clean clone
reproduces the tables and figures exactly. The render takes about 4 minutes: the Monte Carlo
in Exercise 1 is 50,000 replications × 9 cells, each with 200 bootstrap resamples (90 million
resamples in total). Reduce `R_REPS` in the setup chunk to render faster; the conclusions are
unchanged from about 5,000 replications onward, but the reported Monte Carlo standard errors grow.

Built with R 4.5.2, Quarto 1.10.18, fpp2 2.5.1, forecast 9.0.1.

## Layout

| Path | What it is |
|---|---|
| `DATA301-A4-Oliver-Donaldson-300659923.qmd` | the submission; self-contained and renders on its own; all code is in here, and none is echoed in the PDF |
| `working/` | the exploratory scripts behind the design decisions in the submission |
| `install-deps.R` | installs the five packages needed |

### `working/`

These are not needed to render the submission. They are the scripts that settled the choices the
submission then states, kept so the reasoning is auditable.

| Script | Question it answered |
|---|---|
| `01-pilot-mc-error.R` | How many Monte Carlo replications are needed? Pilot at R = 2,000, measuring the relative MC standard error per cell. |
| `02-benchmark-bootstrap.R` | Three ways to compute 200 bootstrap IQRs; the counting-sort version is 5–13× faster than calling `quantile()` per resample. |
| `03-verify-counting-sort.R` | Proof the fast version is exact: counting sort, row sort and `stats::quantile(type = 8)` agree to machine precision on identical resamples. |
| `04-diagnose-n5-failure.R` | Why the bias correction breaks at n = 5: anatomy of the bootstrap mean, share of MSE from the worst replication, and MSE drift as R grows. |
| `05-full-simulation.R` | The full grid as a standalone run: `Rscript 05-full-simulation.R 50000`. |
| `06-explore-ausbeer.R` | Exercise 2 exploration: the plots, the swing-vs-level check (with its residual autocorrelation) and the swing on several Box–Cox scales. |

Scripts 05 and 06 write their outputs (`res_*.rds`, `e2_*.png`) next to themselves; these are
git-ignored.
