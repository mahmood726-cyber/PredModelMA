# Truth-Recovery Validation — PredModelMA

**Date:** 2026-06-14
**Engine under test:** inline functions in `index.html` — `poolCstat` / `poolCstatTransform`
(logit/loglog/arcsine C-stat pooling), `poolOE` (log-scale O:E), `rePool` (RE w/ REML),
`knappHartungC` (HKSJ, t_{k-1}, q-floor at 1).
**Method:** Standalone seeded known-truth Monte Carlo. Inject a TRUE population C-statistic /
O:E ratio, generate k studies with logit-/log-scale heterogeneity (tau) + realistic
within-study SE (Hanley-McNeil), pool with the repo's OWN functions, measure CI recovery
of the true value.

## Verdict: GENUINE, UNUSUALLY COMPLETE ENGINE — POINT ESTIMATES SOUND; DEFAULT z-CI UNDER-COVERS, BUT HKSJ IS PRESENT AND WIRED IN

This is a real engine (Newcombe SE derivation, three C-stat transforms, REML via Fisher
scoring, trim-and-fill, DOI/LFK, Bayesian pooling, PROBAST subgroups). Transforms and
back-transforms are **correct**: logit delta-method SE `se_c/(c(1-c))`, log-scale O:E,
and proper endpoint back-transformation. Point estimates are **unbiased**
(|bias| < 0.0025 in C-units across all cells).

Same coverage signature as its sibling prognostic-meta: **the default C-stat CI is a
z-interval** (`poolCstat` -> `rePool` -> `±normalQuantile*se`) and **under-covers the true
C when heterogeneity is present** (87.7-93.5% vs nominal 95%, worst at small k / high tau).
UNLIKE the sibling, PredModelMA already implements AND surfaces `knappHartungC` (HKSJ,
`t_{k-1}`, `qKH=max(1,Q/(k-1))` — matching the advanced-stats rules) in the rendered output,
so the conservative interval is available to the user.

## Results — C-statistic recovery (C_true = 0.75, 2000 reps/cell)

| tau  | k  | cov (z, default) | cov (HKSJ) | bias(C) | width(z) | width(HKSJ) |
|------|----|------------------|------------|---------|----------|-------------|
| 0.00 | 5  | 97.7%            | 99.8%      | -0.0007 | 0.405*   | 0.426*      |
| 0.15 | 5  | 92.2%            | 98.9%      | -0.0012 | 0.1414   | 0.2259      |
| 0.35 | 5  | 87.7%            | 98.8%      | -0.0020 | 0.1162   | 0.4445      |
| 0.00 | 10 | 97.3%            | 99.4%      | -0.0013 | 0.192*   | 0.200*      |
| 0.15 | 10 | 93.5%            | 98.0%      | -0.0016 | 0.0504   | 0.0898      |
| 0.35 | 10 | 91.8%            | 99.8%      | -0.0024 | 0.0837   | 0.2911      |
| 0.00 | 20 | 96.0%            | 97.5%      | -0.0011 | 0.0614   | 0.0643      |
| 0.15 | 20 | 92.9%            | 98.6%      | -0.0011 | 0.0316   | 0.0545      |
| 0.35 | 20 | 93.3%            | 100.0%     | -0.0014 | 0.0599   | 0.1946      |

Nominal = 95%. Default z-interval: **88-93%** whenever tau>0. HKSJ: 97.5-100% (CONSERVATIVE
— it over-covers under high tau because `sqrt(qKH)` inflation + `t_{k-1}` compound).
*The wide widths at tau=0 are a degenerate-REML artifact (Q<df gives tau2=0 and occasional
very wide intervals on a few replications); coverage is unaffected.

## Results — O:E calibration ratio (O:E_true = 1.0, k=10, 2000 reps)

| tau  | coverage (z) | bias(O:E) |
|------|--------------|-----------|
| 0.00 | 95.9%        | -0.0110   |
| 0.20 | 92.0%        | -0.0113   |

O:E point estimate has a small negative bias (~-0.011) from the Poisson normal-approx in the
DGP combined with the log-CI->SE derivation; well within tolerance. z-interval again
under-covers under heterogeneity.

## What is CORRECT
- logit / loglog / arcsine transforms + delta-method SEs are textbook (Debray 2017).
- O:E pooled on log scale, back-transformed via exp — correct.
- REML (Fisher scoring) tau2; q-floored HKSJ with `t_{k-1}` — matches advanced-stats rules
  (HKSJ floor `max(1,Q/(k-1))`, t not z, df=k-1).
- Point estimates unbiased across all heterogeneity levels for C; near-unbiased for O:E.
- Cross-check: identical pooled C (0.7396) to sibling prognostic-meta on the same 5-study input.

## What is a RISK (not a hard bug)
- **Default reported C-stat / O:E CI is the z-interval (under-covers ~3-7 pts under
  heterogeneity).** Mitigated because `knappHartungC` is computed and displayed alongside,
  but a user reading the primary z-CI as the headline result will be over-confident at
  small k with real between-study heterogeneity.
- HKSJ is markedly conservative under high tau (over-covers to ~99-100%). Acceptable, but
  worth a UI note so users don't read the wide HKSJ CI as "no signal".

## Recommendation
Make the HKSJ (knappHartungC) interval the DEFAULT headline CI for the C-statistic and O:E
when k is small (e.g. k<20) and tau2>0, with the z-interval shown as secondary — or add a
prominent note that the z-CI under-covers under heterogeneity. The machinery is already
present and validated; this is a presentation/default change, not new statistics.

## Reproduce
```
node truth-recovery/harness.mjs             # full coverage sweep + O:E
node truth-recovery/test-truth-recovery.mjs # 6 assertions, exit 0 on pass
```
ADDITIVE only. `engine.mjs` contains the pooling/transform/HKSJ functions extracted VERBATIM
from `index.html` (lines 398-491, 494-654, 960-978), with a 6-line `document` stub so the
two `getElementById('alpha'|'cTransform')` reads resolve outside the browser. No engine math
was modified.
