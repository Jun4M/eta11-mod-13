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
  ~1.24x larger for every `q ≠ 13`. The corrected code now reproduces that table.
- **N_exc** — superseded. `N_exc(10^9) = 123,563` counted the `13|m` stratum; it is
  `96,022` without it.
- **`|Z| >= 3` count** — superseded, `181 of 723` -> `123 of 723`. The old count came
  from the mixed population, matching the old `analyze.py` rather than the δ_q table
  printed above it.
- **fit_alpha block** — superseded; it too was run on the mixed δ_q.
- **factorisation table** — superseded; it used the mixed δ_q, so every `q ≠ 13` row
  was ~1.24x too large.

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

| t (geometric centre) | B | E | n |
|---|---|---|---|
| 1.78e4 | 1.7614 | +0.7614 | 53,249 |
| 5.62e4 | 1.7460 | +0.7460 | 95,008 |
| 1.78e5 | 1.2592 | +0.2592 | 169,680 |
| 5.62e5 | 1.2119 | +0.2119 | 301,540 |
| 1.78e6 | 1.1576 | +0.1576 | 504,556 |
| 5.62e6 | 1.1104 | +0.1104 | 914,243 |
| 1.78e7 | 1.0763 | +0.0763 | 1,892,938 |
| 5.62e7 | 1.0579 | +0.0579 | 2,707,995 |
| 1.78e8 | 1.0399 | +0.0399 | 7,628,828 |
| 5.62e8 | 1.0290 | +0.0290 | 24,124,497 |

The last row was previously suppressed by an off-by-one guard: `MMAX` is
`m[-1] = 999,999,995`, five short of `10^9`, although `999,999,995` is in fact the
largest `m ≡ 11 (mod 24)` below `10^9`, so the window is fully covered.

Two power-law fits, over different window sets:

```
  fit [all windows  , 10 windows]: E(t) = 16.790 * t^(-0.3201)
  fit [t >= 1e6 only, 6 windows]: E(t) = 10.419 * t^(-0.2924)
```

**`E(t)` is not well described by a single power law.** The two smallest windows are
essentially flat (`0.7614`, `0.7460`) and only the large-`t` windows decay, so any
single exponent is a statement about which windows were chosen. Quote the cumulative
count exponent instead — see below. (The previously quoted `t^(-0.28)` is neither of
these; it came from the large-`t` windows together with the complementary relation
`1 - 0.7225`.)

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
  0.3333  1/3                             d(-2logL) =    8.52  REJECTED
  0.2500  1/4                             d(-2logL) =   45.88  REJECTED
  0.4000  2/5                             d(-2logL) =    0.02  consistent
  0.6667  2/3                             d(-2logL) =   69.97  REJECTED
  0.7500  3/4                             d(-2logL) =  109.37  REJECTED
```

The population change moved `alpha` from 0.407 to 0.405 and left the conclusion
intact: `2/5` is still the only surviving simple rational, `1/2` still rejected.
`sigma_hat` moved 0.0543 -> 0.0438, which is just the ~1.24 amplitude inflation
coming out.

## Factorisation check — src/factorisation_check.py

`f_q = delta_q / E` across three decade windows in t, same population in each:

```
E falls 0.1199 -> 0.0620 -> 0.0317  (factor 3.8)
```

| q | w1 | w2 | w3 | mean | sd |
|---|---|---|---|---|---|
| 5 | -1.082 | -1.085 | -0.935 | -1.034 | 0.086 |
| 7 | +0.666 | +0.458 | +0.469 | +0.531 | 0.117 |
| 11 | +0.897 | +1.167 | +0.955 | +1.006 | 0.142 |
| 13 | +0.997 | +1.116 | +0.924 | +1.012 | 0.097 |
| 23 | -0.558 | -0.673 | -0.497 | -0.576 | 0.090 |
| 43 | -0.443 | -0.591 | -0.419 | -0.484 | 0.093 |

Spread is 8–22% of the mean. `sd` is the sample standard deviation (ddof=1); the
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
| 142 (whole sweep) | 0.358 | 0.444 |

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
