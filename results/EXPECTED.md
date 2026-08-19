# Expected output (m <= 10^9)

Reference values from the run of 2026-08-18, `MMAX = 1000000000`
(41,666,667 coefficients).

**This file supersedes the tables dated 2026-08-15.** Those were internally
inconsistent about whether `m ≡ 0 (mod 13)` was in the population: the δ_q table
excluded it, the N_exc table included it, and `analyze.py` did the opposite in both
places. Everything below is measured on **one** population — squarefree kernels with
`13 ∤ t` — so that numerator and denominator of `f_q = δ_q/E` always agree.

What actually changed, block by block:

- **δ_q table** — essentially unchanged. The 2026-08-15 table was already measured
  with `13|m` excluded; it was `analyze.py` that disagreed with it, emitting values
  1.198x larger for every `q ≠ 13`. The corrected code now reproduces that table.
- **N_exc** — superseded. `N_exc(10^9) = 123,563` counted the `13|m` stratum; it is
  `96,022` without it.
- **`|Z| >= 3` count** — superseded, `181 of 723` -> `123 of 723`. The old count came
  from the mixed population, matching the old `analyze.py` rather than the δ_q table
  printed above it.
- **fit_alpha block** — superseded; it too was run on the mixed δ_q.
- **factorisation table** — superseded; it used the mixed δ_q, so every `q ≠ 13` row
  was 1.198x too large.

`data/delta_q_consistent.json` is the reference measurement for this population;
a regenerated `results/delta_q.json` reproduces it bit-exactly on all 723 primes,
and `test_verify.py` now checks that. `data/delta_q_mixed_legacy.json` (formerly
`delta_q_measured.json`) is the mixed-population artifact, retained for provenance
only.

## src/test_verify.py

```
[PASS] reference values  -- p(100)=190569292, 218 digits at n=40000
[PASS] reduction  p(n) = 11*a((24n-1)/13) mod 13  -- tested n <= 40000, n = 6 mod 13; 0 mismatches
[PASS] a13.bin agrees with independent numpy recomputation  -- 200000 coefficients
[PASS] no exact zeros of eta^11 in the tested range  -- checked 200000 coefficients mod 2^31-1
[PASS] Hecke law  a(t p^2) = (lam_p + eps(p)(t|p)p^4) a(t)  -- 91200 congruences tested, 0 failures
[PASS] three-term recursion  a(t p^4) = lam_p a(t p^2) - p^9 a(t)  -- 18833 congruences tested, 0 failures
[PASS] square class of t=155 vanishes entirely  -- 847 coefficients, 0 nonzero
[PASS] delta_q reproduces data/delta_q_consistent.json  -- 723 primes, 723 bit-exact, max|d delta|=0.00e+00, max|d Z|=0.00e+00
all checks passed
```

The first seven are unchanged from the previous run — they never depended on the
population. The eighth is new: it checks a regenerated `results/delta_q.json` against
`data/delta_q_consistent.json`, and is skipped if `analyze.py` has not been run yet
(`test_verify.py` is normally run first).

## Extended verification — src/extended_verify.py

`test_verify.py` keeps its numpy and `2^31-1` checks at 200,000 coefficients so that
it stays a few seconds long. Two claims in the paper need a wider range; both were
run once, on 2026-08-18, and both pass:

```
[PASS] numpy cross-check on 7,692,308 coefficients (m <= 184,615,379)
[PASS] no exact zeros mod 2^31-1 in 3,833,334 coefficients (m <= 92,000,003): 0 zeros
```

So `a13.bin` agrees with an independent numpy recomputation on 7,692,308 coefficients,
and there are no exact zeros of `eta^11` for any `m <= 92,000,003` — the vanishing is
genuinely 13-adic across that whole range, not an artefact of identically zero
coefficients. Takes about 12 minutes; not part of the routine run.

## E(t), half-decade windows

**Restated on squarefree kernels 2026-08-19.** For squarefree `m` the kernel is `m`
itself, so a window's population is fixed once the window is, and these values do not
depend on `MMAX`. `analyze.py` emits them as
`w = sf & (m >= lo) & (m < hi) & (leg13 != 0)`.

| t (geometric centre) | B | E | n |
|---|---|---|---|
| 5.62e4 | 1.3361 | +0.3361 | 2,413 |
| 1.78e5 | 1.2466 | +0.2466 | 7,623 |
| 5.62e5 | 1.2087 | +0.2087 | 24,125 |
| 1.78e6 | 1.1547 | +0.1547 | 76,287 |
| 5.62e6 | 1.1089 | +0.1089 | 241,244 |
| 1.78e7 | 1.0748 | +0.0748 | 762,869 |
| 5.62e7 | 1.0579 | +0.0579 | 2,412,443 |
| 1.78e8 | 1.0399 | +0.0399 | 7,628,828 |
| 5.62e8 | 1.0290 | +0.0290 | 24,124,497 |

The `[1e4, 1e4.5)` window is dropped by the `n >= 2000` guard: 760 kernels give a
standard error of 0.15 on `B`.

**MMAX-independence, verified directly.** Every window shared between the `MMAX = 10^9`
and `MMAX = 10^10` runs is identical in both `n` and its zero count. The two windows
visible only at `10^10`: `[1e9, 1e9.5)` gives `B = 1.0213`, `E = +0.0213`,
`n = 76,288,325`, and `[1e9.5, 1e10)` gives `B = 1.0163`, `E = +0.0163`,
`n = 241,245,018`.

The top four half-decades give `E = +0.0399, +0.0290, +0.0213, +0.0163`, successive
ratios 0.73, 0.73, 0.77 against 0.724 predicted by `t^(-0.28)` over a half-decade. There
is no flattening, so a limit for `E(t)` above about `0.016` is excluded by the data
rather than by extrapolation.

```
  fit [all windows  , 9 windows]: E(t) = 6.773 * t^(-0.2685)
  fit [t >= 1e6 only, 6 windows]: E(t) = 9.680 * t^(-0.2886)
```

At `MMAX = 10^10`: `-0.2633` over twelve windows and `-0.2805` over the eight with
`t >= 1e6`, a divergence of 0.017 against 0.020 here.

`E(t)` still has no single exponent — the window choice moves it by 0.020 — but the
complementary relation `1 - beta = 0.2759` now lies **between** the two fits, and the
rms log-residual over all windows is 0.073.

### Superseded: the all-m definition

Before 2026-08-19 the windows took every `m` whose kernel fell in the range, not the
squarefree `m`. That admits more square multiples `t*s^2` as `MMAX` grows, and since
`a(t) != 0` does not prevent `a(t s^2) = 0`, it raises `B`. It gave

```
  1.78e4 1.7614 | 5.62e4 1.7460 | 1.78e5 1.2592 | 5.62e5 1.2119 | 1.78e6 1.1576
  5.62e6 1.1104 | 1.78e7 1.0763 | 5.62e7 1.0579 | 1.78e8 1.0399 | 5.62e8 1.0290
  fit [all windows  , 10 windows]: E(t) = 16.790 * t^(-0.3201)
  fit [t >= 1e6 only, 6 windows]: E(t) = 10.419 * t^(-0.2924)
```

and it is MMAX-dependent: the `1.78e5` window reads `1.2592` at `10^9` and `1.6269` at
`10^10`, with `n` growing 169,680 -> 533,967. The apparent flatness of the two smallest
windows (`+0.7614`, `+0.7460`) was an artefact of this, and under it `1 - beta` fell
outside both fits rather than between them. The two-column comparison in the paper's
§5.5 keeps both definitions, which is the point there — it is the Corollary 3.3 check.

The `5.62e8` row was also once suppressed by an off-by-one guard: `MMAX` is
`m[-1] = 999,999,995`, five short of `10^9`, although that is the largest
`m = 11 (mod 24)` below `10^9`, so the window is fully covered.

## Cumulative excess count

`N(X) = #{t<=X : a(t)=0} - #{t<=X}/13`, squarefree `t` with `13 ∤ t`.

| X | n | zeros | N_exc |
|---|---|---|---|
| 1e6 | 35,274 | 3,346 | 633 |
| 3.16e6 | 111,561 | 10,122 | 1,540 |
| 1e7 | 352,805 | 30,700 | 3,561 |
| 3.16e7 | 1,115,674 | 93,772 | 7,951 |
| 1e8 | 3,528,117 | 290,094 | 18,700 |
| 3.16e8 | 11,156,945 | 900,347 | 42,120 |
| 1e9 | 35,281,442 | 2,809,979 | 96,022 |

`N_exc = zeros - n/13`, so every row is checkable from the two counts beside it.

fit `N_exc(X) = 0.0296 * X^0.7241`; local slopes 0.70–0.77, mean 0.73.
`beta = 3/4` is rejected (rms log-residual 0.063 against 0.020 for the free fit).

This is the stable exponent: it is a fit to a cumulative count, not to windowed
differences, and it moves by less than 0.002 under the population change.

## Legendre split, top decade (t in [1e8, 1e9))

`E = +0.0317`, n = 31,753,325.
Excluded `13|m` stratum: `E = +0.1193`, n = 2,442,596 — 3.8x the bulk, which is why
mixing the two makes `f_q = δ_q/E` ill-defined.

| q | delta_q | Z | f_q |
|---|---|---|---|
| 5 | -0.0296 | -21.65 | -0.935 |
| 7 | +0.0148 | +11.13 | +0.469 |
| 11 | +0.0302 | +23.21 | +0.955 |
| 13 | +0.0292 | +23.44 | +0.924 |
| 17 | +0.0113 | +8.77 | +0.356 |
| 23 | -0.0157 | -12.34 | -0.497 |
| 43 | -0.0133 | -10.51 | -0.419 |

`f_q = δ_q/E` with the full-precision `E`, not the rounded `0.0317`; dividing by the
rounded value shifts every entry by one in the third decimal.

123 of 723 primes (q < 5500) have |Z| >= 3, against 2.0 expected by chance.

## src/fit_alpha.py on 723 primes

```
n = 723 primes, q in [5, 5483], measurement SE = 0.00125

alpha_hat = 0.405   sigma_hat = 0.0438
  68% interval: [0.380, 0.430]
  95% interval: [0.355, 0.457]

  0.5000  1/2  Euler-factor / Sato-Tate   d(-2logL) =   11.64  REJECTED
  0.3333  1/3                             d(-2logL) =    8.31  REJECTED
  0.2500  1/4                             d(-2logL) =   45.88  REJECTED
  0.3750  3/8                             d(-2logL) =    1.30  consistent
  0.4000  2/5                             d(-2logL) =    0.02  consistent
  0.4167  5/12                            d(-2logL) =    0.23  consistent
  0.6667  2/3                             d(-2logL) =   69.60  REJECTED
  0.7500  3/4                             d(-2logL) =  109.37  REJECTED
```

`3/8` and `5/12` were added on 2026-08-19. All three of `3/8`, `2/5`, `5/12` are
consistent, so `2/5` is not distinguished — the earlier claim that it was the only
surviving simple rational came from a candidate list that omitted its two neighbours.

The population change moved `alpha` from 0.407 to 0.405 and left the conclusion
intact: `2/5` is still the only surviving simple rational, `1/2` still rejected.
`sigma_hat` moved 0.0543 -> 0.0438, which is just the 1.198 amplitude inflation coming
out (0.0543/1.198 = 0.0453, the residual difference being the 7% change in n).

## delta_q extended to q < 20000 (secondary; same m <= 10^9 data)

`results/delta_q_q20000_from1e9.json`, 2260 primes `5 <= q < 20000`, same consistent
population (`n = 31,753,325`, `E = +0.0317`). 30 minutes.

```
n = 2260 primes, q in [5, 19997], measurement SE = 0.00125
alpha_hat = 0.427   sigma_hat = 0.0508
  68% interval: [0.410, 0.445]
  95% interval: [0.390, 0.465]
  0.2500  1/4    d(-2logL) =  121.65  REJECTED
  0.3333  1/3    d(-2logL) =   28.12  REJECTED
  0.3750  3/8    d(-2logL) =    7.90  REJECTED
  0.4000  2/5    d(-2logL) =    2.00  consistent
  0.4167  5/12   d(-2logL) =    0.27  consistent
  0.5000  1/2    d(-2logL) =   12.52  REJECTED
  0.6667  2/3    d(-2logL) =  102.23  REJECTED
  0.7500  3/4    d(-2logL) =  164.12  REJECTED
151 of 2260 primes with |Z| >= 3 (chance 6.1)
```

**alpha drifts upward with the prime range**, from `0.405` at 723 primes to `0.427` at
2260, and `3/8` crosses from consistent (1.30) to rejected (7.90) while `5/12` becomes
the closest candidate. The 723-prime figures remain the paper's primary statement; the
drift is unexplained and is the reason `m <= 10^10` was computed.

An independent measurement of the same quantity agreed on `151 of 2260` exactly and on
six of the eight candidates to the digit (121.65, 7.90, 2.00, 12.52 and both intervals),
differing only at `1/3` (27.85 against 28.12) and `5/12` (0.30 against 0.27). Neither
difference is explained by the evaluation point or by the assumed SE -- varying the SE
moves all eight together and breaks the six that match -- and neither changes a verdict.

`fit_alpha.py` was corrected on 2026-08-19 to profile the likelihood at the exact
candidate rather than at the nearest grid point. `1/3`, `5/12` and `2/3` are the only
candidates off the 0.0025 grid, and they were being evaluated up to 0.00083 away: at
723 primes `1/3` moves 8.52 -> 8.31, `5/12` 0.26 -> 0.23, `2/3` 69.97 -> 69.60. The
on-grid candidates are unaffected.

## Factorisation check — src/factorisation_check.py

`f_q = delta_q / E` across three decade windows in t, same population in each:

```
E falls 0.1199 -> 0.0620 -> 0.0317  (factor 3.8)
```

| q | w1 | w2 | w3 | mean | sd |
|---|---:|---:|---:|---:|---:|
| 5 | -1.082 | -1.085 | -0.935 | -1.034 | 0.086 |
| 7 | +0.666 | +0.458 | +0.469 | +0.531 | 0.117 |
| 11 | +0.897 | +1.167 | +0.955 | +1.006 | 0.142 |
| 13 | +0.997 | +1.116 | +0.924 | +1.012 | 0.097 |
| 17 | +0.454 | +0.251 | +0.356 | +0.354 | 0.101 |
| 23 | -0.558 | -0.673 | -0.497 | -0.576 | 0.090 |
| 37 | +0.553 | +0.386 | +0.350 | +0.430 | 0.108 |
| 43 | -0.443 | -0.591 | -0.419 | -0.484 | 0.093 |
| 73 | +0.674 | +0.431 | +0.436 | +0.514 | 0.139 |
| 89 | -0.379 | -0.338 | -0.412 | -0.376 | 0.037 |
| 101 | -0.298 | -0.383 | -0.405 | -0.362 | 0.056 |

The prime set is every `q` with `|f_q| > 0.3` on the consistent population, 11 primes.
Absolute sd runs 0.037 to 0.142; relative spread `sd/|mean|` is 8–29% of the mean.
Below `|f_q| ~ 0.3` the relative spread diverges because the mean approaches zero, so a
percentage is only meaningful on this restricted set. (`q = 31`, with `|f_q| = 0.263`,
gives the smallest absolute sd of any small prime, 0.027.) `sd` is the sample standard deviation (ddof=1); the
superseded table used the population sd (ddof=0), which is why `q = 13`'s row shows
0.097 here against 0.079 there despite identical entries. `q = 13`'s row is
unchanged by the population fix, as expected — its own class split always excluded
`13|m`.

## Section 4 controls — src/pgen.c and src/controls.py

```bash
gcc -O3 -march=native -funroll-loops -o pgen src/pgen.c
./pgen 30000000            # 61 s, 240 MB, writes res_<l>.bin (11 MB total)
python3 src/controls.py    # 1.0 s at the default bounds
```

`pgen` is cross-checked against `a13.bin`: `res_13.bin` satisfies
`p(n) = 11*a(m) (mod 13)` on all 2,307,692 values, which is the same reduction
`test_verify.py` checks against exact big-integer `p(n)`.

Default bounds are `TMAX = 400000`, `FIT_MAX = 60`; both are printed by every run and
overridable with `--tmax` / `--fitmax`.

```
    l | delta | m mod 24 | worst purity |   1/l   | verdict
    13 |    6  |    11    |   100.00%    |  7.69% | EIGENFORM
    17 |    5  |     7    |   100.00%    |  5.88% | EIGENFORM
    19 |    4  |     5    |   100.00%    |  5.26% | EIGENFORM
    23 |    1  |     1    |   100.00%    |  4.35% | EIGENFORM
    29 |   23  |    19    |   100.00%    |  3.45% | EIGENFORM
    31 |   22  |    17    |   100.00%    |  3.23% | EIGENFORM
    37 |   17  |    11    |     3.06%    |  2.70% | no
    41 |   12  |     7    |     2.69%    |  2.44% | no
    43 |    9  |     5    |     2.62%    |  2.33% | no
    47 |    2  |     1    |     2.43%    |  2.13% | no

    l |  k | predicted e | fitted e | character rule holds? | primes used
    13 |  5 |      4      |       4  |        YES            |   14
    17 |  7 |      6      |       6  |        YES            |   14
    19 |  8 |      7      |       7  |        YES            |   14
    23 | 10 |      9      |       9  |        YES            |   14
    29 | 13 |     12      |      12  |        YES            |   14
    31 | 14 |     13      |      13  |        YES            |   14
```

Tables 2 and 3 hold at 6/6 and do not move with the bounds: at `--fitmax 80` the same
six rows read YES over 19 primes, and the fitted exponents are unchanged at every
`TMAX` tested.

**The purity column is the one bound-dependent figure here.** It is a largest-class
share, biased upward at small sample size, and it decays toward `1/l`:

| TMAX | l=37 | l=41 | l=43 | l=47 |
|---|---|---|---|---|
| 50,000 | 3.53% | 3.23% | 3.23% | 3.02% |
| 150,000 | 3.27% | 2.91% | 2.85% | 2.52% |
| 400,000 | 3.06% | 2.69% | 2.62% | 2.43% |
| 600,000 | 3.02% | 2.69% | 2.63% | 2.40% |
| limit `1/l` | 2.70% | 2.44% | 2.33% | 2.13% |

The 2026-08-15 draft quoted `3.19 / 2.78 / 2.81 / 2.55` with no bound named; those
correspond to no bound this script reproduces and are withdrawn.

## Standing of q = 13, and the nonzero classes

From `src/paper_figures.py`:

```
  raw |f_q| rank of q=13: 3 of 723
  trend-divided |f_q|*q^0.405: q=13 gives 2.61, median 1.21, max 5.32 at q=2971
  rank 85 of 723 -> 88th percentile
  max deviation of a nonzero residue class from 1: 0.0040
```

The median is `1.2146`; an earlier draft rounded it to `1.22`.

## Audit of 2026-08-19 — src/invariants.py, src/independent_check.py, src/paper_figures.py

Every numeric claim in `paper/manuscript.md` was classified (`paper/claims_audit.tsv`,
101 claims): 57 class A (script-emitted and checker-enforced), 33 class B (produced by
no script), 11 class C (literature). Three scripts were added so that class B is empty
of anything the repo can compute.

**Brute force at m <= 10^4.** `prod (1-q^k)^11` built by literal polynomial
multiplication in exact integers — no pentagonal number theorem, no modular reduction,
no blocking — reproduces `a13.bin` on all 417 coefficients. The unreduced head is
`[1, -11, 44, -55, -110]`. The squarefree-kernel sieve, whose modular-inverse
progression offset is the likeliest place for an off-by-one, agrees with direct
trial-division factorisation on all 41,667 values to `m <= 10^6`, and on a random
sample of 20,000 spanning the full range. E, N_exc and delta_5, delta_7 computed by
pure-Python loops match the numpy pipeline exactly.

**Invariants** (`src/invariants.py`, all pass):

```
[PASS] kernel: sq == 1 iff t == m, and t*sq^2 == m for every m
[PASS] E(t) windows partition their span exactly (no gap, no overlap)
[PASS] N(X) = zeros - n/13 reproduces at every X, counts monotone in X
[PASS] delta_q and E share one population
[PASS] B(QR), B(NQR) reconcile with B(all) at weights n1, n2
[PASS] delta_q negates under swap of the two classes
[PASS] the 13 residue classes of a(m) partition the population
[PASS] 13|m and 13 nmid m strata partition the squarefree top decade
```

Invariant 4 is not vacuous: applied to the superseded pairing (delta_q on a population
including `13|m`, normalised by an E without it) it fails by `8.18e-03` against a
tolerance of `1e-12`, and it localises the fault to `q != 13` — for `q = 13` the gap is
`2.22e-16`, because that class split always excluded `13|m`. The bug that produced
eleven corrections would have been caught by one assertion.

**Independent reimplementation** (`src/independent_check.py`, all pass): kernels by
trial-division factorisation rather than modular-inverse progressions; Legendre symbols
by Euler's criterion rather than a table of squares; every rate as an exact integer
ratio rather than a float mean; and `alpha` from binned deconvolved moments rather than
a profiled likelihood, giving `0.420` against the ML `0.405`.

**Result: class A is clean.** Every script-produced figure survived independent
reimplementation, all eight invariants, and brute force at small range. The three
figures the audit corrected were all class B — see the corrections log.

`src/paper_figures.py` now emits the figures that no script produced: the section 5.3
distribution, the section 5.4 kernel table, the N(X) census with its counts, the local
slopes and residuals, the spread of `delta_q/E`, the standing of `q = 13`, the `t = 155`
class and the 426 square classes, the inflation factor, and the Shimura null
distribution (seeded, so reproducible).

## Shimura lift search — src/shimura288.gp

Run 2026-08-19 with PARI/GP 2.17.4, `gp -s 4000000000 -q src/shimura288.gp`, about
five minutes. `mfinit([288,10],0)` completed inside the 4 GB stack; no `mfsplit`
decomposition or trace-form screen was needed.

Correlated against the consistent-population `f_q` over the 60 primes `5 <= q <= 293`
(`results/fq_for_pari.txt`, regenerated from `results/delta_q.json`).

Maximum |correlation| by level: 2:0.030, 3:0.098, 4:0.064, 6:0.127, 8:0.153, 9:0.202,
12:0.141, 16:0.137, 18:0.377, 24:0.151, 32:0.312, 36:0.133, 48:0.264, 72:0.271,
96:0.339, 144:0.276, **288:0.294**.

**Level 288 is a clean negative.** Its maximum is 0.2935 (form 7, embedding 1, signed
-0.2935) over 45 (form, embedding) pairs; the global maximum over all 142 pairs is
0.3772 at level 18 form 4.

That is not merely below the 0.5 noise threshold — it is at or below what noise alone
produces. At `n = 60` points a single correlation has SE `0.130`, and the null
distribution of the *maximum* over many candidates is:

| candidates | null median max \|corr\| | null 95th pct |
|---|---|---|
| 1 (a single form) | 0.091 | 0.262 |
| 45 (level 288) | 0.313 | 0.412 |
| 142 (whole sweep) | 0.359 | 0.444 |

Level 288's 0.294 is **below** the null median for 45 candidates, and the global 0.377
is essentially at the null median for 142. So no level shows any evidence of the lift,
and the previously reported "max 0.41" across levels <= 144 was itself within noise.

The earlier per-level maxima (2:0.09, 3:0.05, ... 144:0.35) are superseded: they were
scored against a `results/fq_for_pari.txt` that predated the population fix. That file
was not a rescaling of the corrected `f_q` — the ratio ranges over -0.058 to 2.216
across the 60 primes, including a sign change — so an earlier note claiming the search
was immune by scale-invariance was wrong in its premise. Re-running on the corrected
`f_q` moves level 288 from 0.288 to 0.294 and the global maximum from 0.368 to 0.377,
changing no conclusion.
