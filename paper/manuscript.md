# Vanishing of eta^11 modulo 13 along square classes: a computation to 10^10

*Draft v3, corrected 2026-08-18. Supersedes all earlier versions, which are in the git
history. Numbers are on a single consistent population (squarefree kernels coprime to
13); earlier drafts mixed two populations in one place, which inflated the Legendre
correlations by a factor 1.198. Every figure below is cross-checked against
`results/EXPECTED.md`; the two must never differ.*

---

## Abstract

Let `p(n)` be the partition function. For `n ≡ 6 (mod 13)` put `m = (24n-1)/13`, so
that `p(n) ≡ 11·a(m) (mod 13)`, where `a(m)` is the `m`-th coefficient of
`∏(1-q^k)^11`. We compute `a(m) mod 13` for all `m ≤ 10⁹` (41,666,667 coefficients)
and measure the density with which those coefficients vanish along squarefree kernels.

Writing `m = t s²` with `t` squarefree, the statistics are organised by `t`, as the
Hecke structure requires — though `a(t) ≠ 0` does not prevent `a(t s²) = 0`, so a
density must be measured on a population fixed independently of the computation's range
(§6.1). The density of `t` with `a(t) ≡ 0 (mod 13)` exceeds `1/13`
but converges toward it: normalised so that `1` is the random baseline, the ratio
falls from `1.3361` to `1.0290` across four and a half decades, and the cumulative excess
count obeys `N(X) = 0.0296·X^0.7241`, or `0.0296·X^0.7242` on extending the computation to
`m ≤ 10¹⁰`, with the exponent `3/4` rejected.

Within this excess, the density depends on the Legendre symbol `(t|q)` — for
essentially every small prime `q`, not only for `q = 13`. The dependence factorises:
`δ_q(t) = E(t)·f_q` to within 8–29% across three decades in `t` over which `E` falls
by a factor 3.8. The `q`-dependence is `f_q ∝ q^(-0.405)` with 95% interval
`[0.355, 0.457]` over 723 primes to `q = 5483`; `1/4`, `1/3`, `2/3` and `3/4` are
rejected at every range measured, and `1/2` is too, though by a margin that narrows as the
data improves. We name no rational: the estimate is not settled.

**The two exponents are not on the same footing, and that contrast is this paper's
clearest statement.** Across a decade of additional data `β` moves by `0.0001` and `3/4`
is rejected more strongly. Over the same extension `α` moves by `0.025`, from `0.405` to
`0.430`, with no narrowing of its 95% interval — width `0.102` against `0.100` — because
that width is set by how many auxiliary primes exist and not by the range of the
computation, so no amount of `m` will settle it. One of the two exponents is a
measurement; the other is not yet.

Existing theory reaches none of this: it gives lower bounds only, the best for
non-vanishing being `≫ √X/log log X` (Bellaïche–Green–Soundararajan) against an
observed `(12/13)X`. We report measurements and a factorisation; the mechanism behind
`β` and behind `α` is open in both cases.

---

## 1. What this paper does and does not claim

This is a computational paper containing no new theorems. Its contributions:

1. an explicit, coefficient-level verification of the Folsom–Kent–Ono eigenform
   congruence at `ℓ = 13`, with eigenvalues, weight exponent and character identified
   numerically, and the boundary `ℓ ≤ 31` measured as a sharp dichotomy (§3, §4);
2. the density of `13 | p(n)` along squarefree kernels out to `m = 10⁹`, with the
   rate at which it approaches `1/13` (§6.1);
3. the factorisation of the Legendre-class dependence and the measurement of `f_q`
   over 723 primes (§6.2, §6.3).

The third is the substance. It addresses a question on which nothing is proved:
whether the coefficients of a half-integral weight form equidistribute modulo `ℓ`,
and if so how fast.

An earlier version of this work claimed a bias specific to `ℓ = 13`. That claim is
withdrawn: extending the measurement to 723 auxiliary primes shows `q = 13` is high
but not an outlier within the family (§6.3). *(The earlier, narrower measurement it
replaced is not reproducible from the archived data and is not relied on here.)*

---

## 2. Reduction

For `n ≡ 6 (mod 13)` write `m = (24n-1)/13`, so `m ≡ 11 (mod 24)`. Let

```
Σ a(m) q^m  =  ∏_{k≥1} (1 - q^k)^11                          [OEIS A010819]
```

whose exponents are exactly the integers `≡ 11 (mod 24)`. Then

```
p(n) ≡ 11 · a(m)   (mod 13).                                              (R)
```

This is [FKO, §6], where it appears as `P₁₃(1; z) ≡ 11·η(z)^11 (mod 13)`. Verified
directly for all `n ≤ 40000`.

Since `(1-x)^13 ≡ 1 - x^13 (mod 13)`, one also has
`Σ a(m)q^m ≡ ∏(1-q^{13k})·∏(1-q^k)^{-2} (mod 13)`. We do not use this, but it makes
the elementary nature of the reduction visible.

The analogous single-eta-power reduction fails for every other prime: we tested
`p(ℓi+δ_ℓ) ≡ c·[qⁱ]∏(1-q^k)^{ℓ-2} (mod ℓ)` for all nine primes
`ℓ = 17, 19, 23, 29, 31, 37, 41, 43, 47` over 120,000 coefficients each, and in every
case it fails already at `i = 1`, while holding identically at `ℓ = 13`. Control
computations therefore run on `p(n)` directly.

Because `m ≡ 11 (mod 24)`, every `m` is odd and coprime to 3, so square factors come
only from primes `≥ 5`.

### Population convention

Throughout, statistics are computed on **squarefree kernels `t` coprime to 13**. The
stratum `13 | t` behaves quite differently — `E = +0.1193` against `+0.0317` for the
bulk in the top decade, a factor 3.8 — and mixing the two makes the ratio `δ_q/E`
ill-defined, since the factorisation of §6.2 presupposes a population with uniform
`E`. The stratum `13 | t` is reported separately where relevant and never folded in.
For `q = 13` the class `(t|13) = 0` is structurally excluded in any case, so
excluding it for all `q` is what makes the 723 primes comparable.

---

## 3. The Hecke structure, explicitly

[FKO, Theorem 1.3] states that `P_ℓ(b; 24z) (mod ℓ^m)` is an eigenform of the
half-integral weight Hecke operators for `5 ≤ ℓ ≤ 31`, with `T(c²)` acting by
[FKO, (1.4)]. At `ℓ = 13`, `b = m = 1` this specialises to a congruence between
coefficients that can be checked directly. We record the explicit constants because
they do not appear in the literature and are needed below.

### Observation 3.1

*For squarefree `t ≡ 11 (mod 24)` and every prime `p ≠ 13`,*

```
a(t p²)      ≡  ( λ_p + ε(p)·(t|p)·p⁴ ) · a(t)                (mod 13)
a(t p^{2j})  ≡  λ_p · a(t p^{2j-2})  -  p⁹ · a(t p^{2j-4})    (mod 13),   j ≥ 2
```

*where `ε(p) = +1` if `p ≡ 2 (mod 3)` and `-1` if `p ≡ 1 (mod 3)`, and*

| p | 5 | 7 | 11 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 | 47 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| λ_p | 10 | 8 | 5 | 1 | 8 | 8 | 4 | 4 | 5 | 9 | 12 | 6 |

**Verification.** The first relation was tested on every squarefree `t ≡ 11 (mod 24)`
with `t < 2·10⁵` and `t p²` in range — 7,600 cases for each of the twelve primes,
91200 congruences in all, the count being uniform in `p` because `t < 2·10⁵` binds
before `t p² ≤ 10⁹` does. The second was tested over the same kernels, 18833
congruences. No exceptions in either. The character rule was
checked separately on all 50 primes `5 ≤ p ≤ 250`, `p ≠ 13`, without exception; within
each class of `p mod 3` the residue `p mod 4` varies freely (`p = 37, 61, 73, 97, 109,
157, 181, 193` are all `≡ 1 mod 4` and all give `ε = -1`), so the conductor is 3, not 12.

The exponents are not fitted. `η^11` has weight `11/2`, so `λ = 5` in [FKO, (1.4)],
giving `p^{λ-1} = p⁴` and `p^{2λ-1} = p⁹` exactly as observed; and `ε(p)` is
`χ₁₂(p)·((-1)^λ|p)` up to sign, since `χ₁₂(p) = (3|p) = (p|3)(-1|p)`.

`λ_p` does not factor through `p mod 13`: the pairs `(5,31)` and `(17,43)` are
congruent mod 13 with different `λ_p` (10 against 4, and 1 against 12), so the system
is not Eisenstein. The converse also occurs — `11 ≡ 37 (mod 13)` with `λ₁₁ = λ₃₇ = 5` —
so `p mod 13` neither determines `λ_p` nor is determined by it. We did not
identify the Shimura lift; it is not an eta quotient on the divisors of 24 (exhaustive
search over weight-10 quotients with two or three divisors), and not of level 1
(`Δ`, `Δ·E₄` … `Δ·E₁₄`, all twists `p^j`, `0 ≤ j ≤ 11`, all at chance level).
*(Both searches predate the current repository and are not reproduced by any script in
it; they are recorded as unverified. The level-288 search below is reproducible.)*

Since `η(24z)^11` lies on `Γ₀(576)` and the Shimura lift maps `S_{k+1/2}(4N)` to
`S_{2k}(2N)`, the lift has level dividing 288. **Every such level has now been
searched, and none matches.** Scoring each weight-10 newform by the correlation
between its normalised eigenvalues `a_q/q^{9/2}` and the measured `f_q` over the 60
primes `5 ≤ q ≤ 293`, the maximum |correlation| is `0.294` at level 288 (45 form–
embedding pairs) and `0.377` over all 142 pairs at all 17 levels dividing 288.

These are not near-misses; they are what noise produces. At `n = 60` a single
correlation has standard error `0.130`, and the null distribution of the maximum over
`N` candidates has median `0.313` for `N = 45` and `0.359` for `N = 142`. The level-288
maximum falls *below* the null median for its own candidate count. So the lift is
either not a newform of level dividing 288, or `f_q` is not a normalised eigenvalue
sequence — and the mechanism behind `f_q` remains open.

### Corollaries

**3.2 (propagation).** If `a(t) ≡ 0` then `a(t s²) ≡ 0` for all `s` coprime to 6.
In terms of `p(n)`: if `13 | p(n₀)` with `24n₀-1 = 13t`, `t` squarefree, then
`13 | p(n)` for every `n` with `24n-1 = 13ts²`. Smallest instance `t = 155`, `n₀ = 84`:

```
n = 84, 2099, 4114, 10159, 14189, 24264, 30309, 44414, 52474, …
```

verified for 847 members of this class, the largest at `n = 541,239,159`. Over all
kernels `t < 10⁵` with `a(t) ≡ 0`, all 426 square classes vanish entirely, without
exception. These generators are recorded as OEIS A399067.

**3.3 (kernel dependence).** `a(ts²)` is `a(t)` times a factor depending only on `t`,
`s` and the eigenvalues, so divisibility statistics on a square class are determined
by its kernel. Verified numerically in §5.5.

---

## 4. The boundary at ℓ = 31

Observation 3.1 can be tested directly on `p(n)`. For prime `ℓ` set
`δ_ℓ ≡ 24^{-1} (mod ℓ)` and `A_ℓ(m) := p((ℓm+1)/24) mod ℓ`; the prediction is that
`A_ℓ(tp²)/A_ℓ(t) (mod ℓ)` is constant on each class of `(t|p)`. If the series is not
congruent to a single eigenform the ratio should be equidistributed, with the largest
class holding about `1/ℓ`.

Computed for `n ≤ 3·10⁷` over `p ∈ {5,7,11,17,19}` (excluding `p = ℓ`), worst case
over `p` and over the classes `(t|p) = ±1`, at kernel bound `t < 4·10⁵`:

| ℓ | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 | 47 |
|---|---|---|---|---|---|---|---|---|---|---|
| purity | 100% | 100% | 100% | 100% | 100% | 100% | 3.06% | 2.69% | 2.62% | 2.43% |
| `1/ℓ` | | | | | | | 2.70% | 2.44% | 2.33% | 2.13% |

**The purity percentages for `ℓ ≥ 37` are finite-sample estimates and must be quoted
with their bound.** "Largest class share" is biased upward at small sample size and
decays toward `1/ℓ` as the bound grows; at `t < 5·10⁴, 1.5·10⁵, 4·10⁵, 6·10⁵` the
`ℓ = 37` figure reads `3.53, 3.27, 3.06, 3.02` percent. Only the `100%` column and the
dichotomy itself are bound-independent, and they are what this table asserts. (An
earlier draft quoted `3.19 / 2.78 / 2.81 / 2.55` with no bound named; those values are
withdrawn — they correspond to no bound this script reproduces.)

The dichotomy is total, and its location is exactly where [FKO] have
`⌊(ℓ-1)/12⌋ - ⌊(ℓ²-1)/24ℓ⌋ = 1`, i.e. the relevant space is one-dimensional. At
`ℓ = 37` the count is `3 - 1 = 2` and the property fails.

The law is uniform where it holds. Fitting the exponent independently for each `ℓ`
over `5 ≤ p < 60` — 14 auxiliary primes — gives `e = k - 1` with `k = (ℓ-3)/2` in
every one of the six cases:

| ℓ | 13 | 17 | 19 | 23 | 29 | 31 |
|---|---|---|---|---|---|---|
| `k = (ℓ-3)/2` | 5 | 7 | 8 | 10 | 13 | 14 |
| fitted exponent | 4 | 6 | 7 | 9 | 12 | 13 |

and the character, extracted with the same fitted exponent, is uniformly
`ε_ℓ(p) = -(p|3)·(-1|p)^{k+1}` — conductor 3 for `ℓ = 13, 17, 29` (`k` odd) and 12 for
`ℓ = 19, 23, 31` (`k` even). Both hold at 6/6 over `5 ≤ p < 60`, and again over
`5 ≤ p < 80` (19 primes). Unlike the purity row, neither table moves with the kernel
bound.

Both are reproduced by `src/pgen.c` and `src/controls.py`; see `results/EXPECTED.md`.

The weight behaves as if the governing form were `η^{ℓ-2}` in every case, even though
the reduction (R) itself holds only at `ℓ = 13`.

---

## 5. The computation

### 5.1 Method

`a(m) mod 13` was computed for all `m ≤ 10⁹` — 41,666,667 coefficients — by eleven
successive multiplications of the pentagonal series, in C, coefficients stored as
`int8`, with a cache-blocked accumulation loop.

Independent checks:

1. against a numpy implementation using different blocking, on 7,692,308 overlapping
   coefficients: exact agreement;
2. against exact big-integer `p(n)` (218 decimal digits at the top) for all
   `n ≤ 40000`, via (R): exact agreement;
3. reference values `p(10) = 42`, `p(54) = 386155`, `p(100) = 190569292`;
4. re-run from scratch by an independent implementation: all eight checks reproduced,
   every count identical, and the 723 `δ_q` and `Z` values agree bit-for-bit.

### 5.2 The vanishing is genuinely 13-adic

Recomputing the coefficients modulo `2³¹ - 1` gives **no exact zeros** for any
`m ≤ 92,000,003` (3,833,334 coefficients). The observed vanishing is not an artefact of
identically zero coefficients.

### 5.3 Only the zero class is anomalous

Distribution of `a(t) mod 13` over squarefree `t`, normalised so that `1` is the
random baseline:

| value | `(t\|13) = +1` | `(t\|13) = -1` | `13 \| t` |
|---|---|---|---|
| **0** | **1.052** | **1.018** | **1.132** |
| 1–12 | 0.994 – 0.997 | 0.996 – 1.000 | 0.981 – 0.996 |

17,640,743 and 17,640,699 samples in the first two columns, 2,714,016 in the third —
all squarefree `m ≤ 10⁹`. The largest deviation of a non-zero class from 1 is 0.0040,
so the twelve non-zero classes are uniform to within 1%: the
phenomenon is an excess of vanishing, not a skew in the distribution of values.

### 5.4 Independent verification of the pipeline

Agreement between two implementations that share an algorithm proves nothing. The
population-mixing error corrected in draft v3 was emitted by the analysis script,
recorded in the expected output, and quoted here — consistent in all three places and
wrong in all three. Three checks were therefore run that do not share an algorithm with
the pipeline (`src/invariants.py`, `src/independent_check.py`, and brute force):

1. **Brute force at `m ≤ 10⁴`.** `∏(1-q^k)^11` built by literal polynomial
   multiplication in exact integers — no pentagonal number theorem, no reduction, no
   blocking — reproduces all 417 coefficients. The squarefree-kernel sieve, whose
   modular-inverse progression offset is where an off-by-one would hide, agrees with
   direct trial-division factorisation on all 41,667 values to `m ≤ 10⁶` and on a random
   sample of 20,000 spanning the full range.
2. **Reimplementation by different algorithms.** Kernels by trial division rather than
   modular-inverse progressions; Legendre symbols by Euler's criterion rather than a
   table of squares; every rate as an exact integer ratio rather than a float mean; and
   `α` from binned deconvolved moments rather than a profiled likelihood, giving `0.420`
   against the maximum-likelihood `0.405`.
3. **Invariants.** Twelve assertions that must reconcile: the windows partition their
   span, `N(X) = zeros − n/13` at every `X`, the two Legendre classes reconcile with the
   pooled rate at weights `n₁, n₂`, `δ_q` negates under swapping the classes, the
   thirteen residue classes partition the population, vanishing propagates up every
   square class, and `δ_q` and `E` are measured on one population. That last is the
   assertion the v3 error would have failed, by `8.2·10⁻³` against a tolerance of
   `10⁻¹²`.

All three passed on every figure the scripts produce — and so did the `E(t)` window
definition of §6.1, which was still not the quantity the section claimed, and took a
second computation at ten times the range to expose. The seventeen corrections below are
of six kinds, three of which were invisible to every mechanism in place when they were
introduced. Appendix A lists them, and Appendix B states what these checks do not
establish — the three above were written by the same author as the pipeline, so they
differ in algorithm and not in origin.

### 5.5 The statistics depend only on the kernel

Contrast computed over all `m` with a given kernel, against the same computed over
squarefree `m` only (where `t = m`, so `m` is far smaller):

| kernel window | all `m` | squarefree `m` only |
|---|---:|---:|
| `[10⁵, 10⁶)` | +0.2289 | +0.2178 |
| `[10⁶, 10⁷)` | +0.1271 | +0.1199 |
| `[10⁷, 10⁸)` | +0.0655 | +0.0620 |
| `[10⁸, 10⁹)` | +0.0317 | +0.0317 |

Agreement to four decimals in the largest window, where the squarefree sample is
largest; the smaller windows differ by 5–6% because the all-`m` column there is
dominated by kernels drawn from a wider range of `m`. This is Corollary 3.3 measured.

The all-`m` column is the one that grows with `MMAX`: at `MMAX = 10¹⁰` its `[10⁵,10⁶)`
entry rises well above the value here, while the squarefree column does not move at all.
Both belong in this comparison — the difference between them *is* the measurement — but
only the squarefree column is a range-free quantity, which is why §6.1 uses it. See §6.1
for the size of the drift.

---

## 6. Results

### 6.1 The excess vanishes

Normalised density `B` of `a(t) ≡ 0 (mod 13)` on **squarefree** kernels coprime to 13,
so that `t = m` and the population of a window is fixed once the window is:

| `t` window | `B` | n |
|---|---:|---:|
| `[10⁴·⁵, 10⁵)` | 1.3361 | 2,413 |
| `[10⁵, 10⁵·⁵)` | 1.2466 | 7,623 |
| `[10⁵·⁵, 10⁶)` | 1.2087 | 24,125 |
| `[10⁶, 10⁶·⁵)` | 1.1547 | 76,287 |
| `[10⁶·⁵, 10⁷)` | 1.1089 | 241,244 |
| `[10⁷, 10⁷·⁵)` | 1.0748 | 762,869 |
| `[10⁷·⁵, 10⁸)` | 1.0579 | 2,412,443 |
| `[10⁸, 10⁸·⁵)` | 1.0399 | 7,628,828 |
| `[10⁸·⁵, 10⁹)` | **1.0290** | 24,124,497 |

Monotone, with no sign of levelling off. The `[10⁴, 10⁴·⁵)` window is omitted: 760
kernels give a standard error of 0.15 on `B`, which is not a measurement.

**These values do not depend on the range of the computation, and the earlier
all-`m` version did.** Taking every `m` whose kernel falls in the window, rather than
the squarefree `m`, admits more square multiples `t·s²` as `MMAX` grows — and since
`a(t) ≠ 0` does not prevent `a(t s²) = 0`, admitting more multiples raises `B`. The
`[10⁵, 10⁵·⁵)` window read `B = 1.2592` at `MMAX = 10⁹` and `1.6269` at `10¹⁰`, with `n`
growing from 169,680 to 533,967. On squarefree kernels every window shared by the two
computations is identical in both `n` and its zero count. The all-`m` quantity is a
different thing, and a range-dependent one; it belongs only in §5.5, where the
comparison between the two columns is the point.

`E(t)` still has no single exponent: fitting all nine windows gives
`E(t) = 6.773·t^(-0.2685)`, the six with `t ≥ 10⁶` give `9.680·t^(-0.2886)`, and the
choice of windows moves the exponent by 0.020. But the restatement removes most of the
apparent pathology. Under the all-`m` definition the two smallest windows looked flat
(`E = 0.7614`, `0.7460`) and the fits diverged by 0.028 with an rms log-residual of
0.167 over all windows; on squarefree kernels the small windows continue the decay
(`0.3361`, `0.2466`) and the rms falls to 0.073. **The flatness was an artefact of the
definition, in the same windows where the definitions differ most.**

The complementary relation now holds. If `E(t) ~ t^{-γ}` exactly then summing gives
`N(X) ~ X^{1-γ}`, so `γ = 1 - β = 0.2759` — which lies *between* the two fitted
exponents, `0.2685` and `0.2886`. Under the all-`m` definition it fell outside both
(`0.2927` and `0.3202`), and an earlier draft read that as evidence against a power law.
Extending to `MMAX = 10¹⁰` narrows the two fits further, to `-0.2633` over twelve windows
and `-0.2805` over the eight with `t ≥ 10⁶`, a divergence of 0.017. The two windows that
range adds give `E = +0.0213` and `+0.0163`, continuing the decay at the same rate: the
top four half-decades read `+0.0399, +0.0290, +0.0213, +0.0163`, successive ratios 0.73,
0.73 and 0.77 against 0.724 for a `t^(-0.28)` law.

We still quote `β` rather than `γ`, because `β` is a fit to a cumulative count and does
not depend on a windowing choice at all.

The stable quantity is the cumulative excess count
`N(X) := #{t ≤ X squarefree, 13 ∤ t : a(t) ≡ 0} - #{t ≤ X}/13`:

| X | `#{t ≤ X}` | `#{a(t) ≡ 0}` | N(X) |
|---|---:|---:|---:|
| 10⁷ | 352,805 | 30,700 | 3,561 |
| 10⁸ | 3,528,117 | 290,094 | 18,700 |
| 10⁹ | 35,281,442 | 2,809,979 | 96,022 |

with `N(X) = 0.0296·X^0.7241`, local slopes in `0.70–0.77`. `β = 3/4` is rejected: the
rms log-residual is 0.063 against 0.020 for the free fit.

**`β` holds across a further decade.** Extending to `MMAX = 10¹⁰` adds two rows —
`X = 3.16·10⁹` with 111,569,767 kernels and 8,803,158 vanishing, giving `N = 220,868`,
and `X = 10¹⁰` with 352,814,785 kernels and 27,663,159 vanishing, giving `N = 523,560` —
and the fit becomes `0.0296·X^0.7242` — the exponent moves by 0.0001 and the prefactor not
at all. Every value for `X ≤ 10⁹` is unchanged. `β = 3/4` is rejected more strongly, rms
0.079 against 0.019 for the free fit. This is the paper's most stable measurement: nine
half-decades of a cumulative count, no windowing choice, and no sensitivity to the range
of the computation.

### 6.2 The Legendre dependence factorises

With `δ_q := 13·[P(a(t)≡0 | (t|q)=+1) - P(a(t)≡0 | (t|q)=-1)]`, the ratio `δ_q/E`
is constant across three decades in `t` over which `E` falls by a factor 3.8:

| q | `t ∈ [10⁶,10⁷)` | `[10⁷,10⁸)` | `[10⁸,10⁹)` | mean | sd |
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

`E` itself falls `+0.1199 → +0.0620 → +0.0317` across the three windows. The spread of
`δ_q/E` is 8–29% of the mean; `sd` is the sample standard deviation (`ddof = 1`) of the
three window values. The `q = 13` row is unchanged from the mixed-population draft, as
it must be: the class `(t|13) = 0` was always excluded there.

So

```
P( a(t) ≡ 0 mod 13 | (t|q) = ε )  =  (1/13)·[ 1 + E(t)·(1 + f_q·ε/2) ]
```

with `E` depending only on `t` and `f_q := δ_q/E` only on `q`. There is one excess
phenomenon, partitioned between the Legendre classes in a `q`-dependent but
`t`-independent ratio.

### 6.3 `q = 13` is high but not an outlier

On the top decade (`E = +0.0317`, `n = 31,753,325`):

| q | δ_q | Z | f_q |
|---|---:|---:|---:|
| 5 | -0.0296 | -21.65 | -0.935 |
| 7 | +0.0148 | +11.13 | +0.469 |
| 11 | +0.0302 | +23.21 | +0.955 |
| **13** | **+0.0292** | **+23.44** | **+0.924** |
| 17 | +0.0113 | +8.77 | +0.356 |
| 23 | -0.0157 | -12.34 | -0.497 |
| 43 | -0.0133 | -10.51 | -0.419 |

`q = 11` gives a slightly larger split than `q = 13`; `q = 5` is comparable with the
opposite sign. Over 723 primes to `q = 5483`, 123 have `|Z| ≥ 3` (chance would give 2).

`q = 13` is high in the family but not an outlier. By raw `|f_q|` it ranks **3rd of
723** — but that ranking is confounded by the decay of §6.4, since `f_q` falls like
`q^{-0.405}` and the small primes therefore dominate any raw ordering. Dividing the
trend out, `|f_q|·q^{0.405} = 2.61` for `q = 13` against a median of 1.21 and a maximum
of 5.32 (attained at `q = 2971`), which places it at the **88th percentile** — above
the middle of the distribution, and well inside it. What is withdrawn is the claim that
`ℓ = 13` plays a distinguished role in its own Legendre statistics; what is not claimed
is that it is average.

### 6.4 The decay of `f_q` is not determined

Fitting `δ_q = c_q·q^{-α}` with `c_q ~ N(0,σ²)` and the binomial measurement noise,
constant in `q`, modelled explicitly. Four measurements:

| population | primes | `α` | 95% interval | width | SE | `σ` |
|---|---:|---:|---|---:|---:|---:|
| `m ≤ 10⁹`, `t ∈ [10⁸,10⁹)` | 723 | 0.405 | [0.355, 0.457] | 0.102 | 0.00125 | 0.0438 |
| `m ≤ 10⁹`, `q < 20000` | 2260 | 0.427 | [0.390, 0.465] | 0.075 | 0.00125 | 0.0508 |
| `m ≤ 10¹⁰`, `t ∈ [10⁹,10¹⁰)` | 723 | 0.430 | [0.382, 0.482] | 0.100 | 0.00039 | 0.0225 |

The 68% intervals are `[0.380, 0.430]`, `[0.410, 0.445]` and `[0.407, 0.455]`, and the
counts of primes at `|Z| ≥ 3` are 123 of 723, 151 of 2260, and 205 of 723 — against 2.0,
6.1 and 2.0 expected by chance. The third population is a different one: at `MMAX = 10¹⁰`
the top decade is `t ∈ [10⁹, 10¹⁰)`, where `E = +0.0175` on `n = 317,533,343` kernels,
with the excluded `13|m` stratum at `E = +0.0578` on 24,425,646.

A fourth line, from a different estimator entirely — binned deconvolved moments rather
than a profiled likelihood, on the first population — gives `0.420`.

**The estimate is not settled.** It moves in one direction under two independent
levers: extending the auxiliary primes by a factor 3 moves it `0.405 → 0.427`, and
extending `m` by a decade moves it `0.405 → 0.430`. These change different things — the
second measures `f_q` on a different range of `t` altogether — and they agree. The
earlier value was low.

**The interval does not narrow, and this is structural.** Extending `m` by a decade cut
the measurement error by a factor 3.2, from `0.00125` to `0.00039`, and moved the 95%
width from `0.102` to `0.100`. That is not a shortfall of data. The model treats `c_q` as
random with spread `σ`, so `α` is estimated from the scatter of one draw per prime, and
its precision is set by how many primes there are rather than by how well each is
measured. Two checks confirm it: holding the primes fixed and shrinking the assumed SE
tenfold moves the width only to `0.082`, while holding the SE and varying the prime count
gives widths `0.202`, `0.148`, `0.100` at 180, 361 and 723 primes — a clean `1/√n`.

The consequence is worth stating as a result rather than as a limitation. **`α` is not
determined by this estimator at any range of `m`**, because `m` supplies no additional
primes. The only lever that narrows the interval is the number of auxiliary primes, and
that is bounded in practice by the cost of `δ_q` at large `q`: reaching a width of `0.02`
would need of order 10⁵ primes.

We therefore name no rational. Against the candidates:

| candidate | `m ≤ 10⁹`, 723 | `m ≤ 10⁹`, 2260 | `m ≤ 10¹⁰`, 723 |
|---|---:|---:|---:|
| 1/4 | 45.88 rejected | 121.65 rejected | 69.29 rejected |
| 1/3 | 8.31 rejected | 28.12 rejected | 17.36 rejected |
| 3/8 | 1.30 consistent | 7.90 rejected | 5.30 rejected |
| 2/5 | 0.02 consistent | 2.00 consistent | 1.55 consistent |
| 5/12 | 0.23 consistent | 0.27 consistent | 0.32 consistent |
| 1/2 | 11.64 rejected | 12.52 rejected | 6.74 rejected |
| 2/3 | 69.60 rejected | 102.23 rejected | 63.20 rejected |
| 3/4 | 109.37 rejected | 164.12 rejected | 105.56 rejected |

`3/8` is now rejected at both extensions and `5/12` is closest at both — but `α̂ = 0.430`
is nearer `3/7 = 0.4286` than `5/12 = 0.4167`, and `3/7` is not in the table. That is the
signature of a moving estimate acquiring whichever rational happens to lie next to it. An
earlier draft called `2/5` "the only simple rational not rejected"; one range later the
same reasoning would name `5/12`, and one range further it may name something else. We
do not repeat it.

**The `1/2` rejection is weakening.** It is rejected at every range measured — `11.64`,
`12.52`, `6.74` — but the trend is the number to report, not the endpoint: the better the
data, the weaker the rejection, because the estimate is moving toward `1/2` faster than
the interval is tightening. `1/2` is what the obvious heuristic predicts, if
`a(t)² ~ L(1/2, g, χ_t)` by Waldspurger and the Euler factor at `q` contributes an
amplitude `λ̃_q·q^{-1/2}` with `λ̃_q ∈ [-2,2]`. We take the heuristic to be **unsettled
rather than refuted**. An earlier draft stated flatly that the picture "is not supported";
on the present evidence that is stronger than the data warrants.

Nothing else in this paper depends on `α`. The factorisation of §6.2 is a statement about
`δ_q/E` being constant in `t`, measured directly. The Shimura search of §3 correlates
against the measured `f_q` and is invariant under any rescaling of them, so its negative
holds whatever the exponent turns out to be. And `β` is unaffected — see below.

*(Caution: an earlier analysis of 60 primes gave `α ≈ 1/2` and a distribution consistent
with Sato–Tate. Restricting the present data to 60 primes reproduces it — `α = 0.495`
with `1/2` consistent at `d(-2logL) = 0.00`, against `0.405` and `1/2` rejected at 723 —
so it was a small-sample artefact, and it is the same shape as the drift documented here.
The `√q` growth of the noise in `δ_q·√q` is the trap; it must be modelled, not eyeballed.)*

### 6.5 Reading

`E(t)` decreases monotonically over the four and a half decades the squarefree windows
cover, and over five and a half at `MMAX = 10¹⁰`; `f_q` decreases in `q`. The data is
consistent with

```
density{ t squarefree : a(t) ≡ 0 (mod 13) }  =  1/13,
```

with the excess and its Legendre partition both being finite-range correlations. This
is a reading of the data, not a theorem; the decay does not exclude a small positive
limit, though it now constrains one — see §9.

---

## 7. Relation to known results

The literature is uniformly existential or a lower bound.

- **Ono–Skinner** (Ann. of Math. 147, 1998): infinitely many squarefree `m` with
  `v_ℓ(a(m)) = 0`; Corollary 4 of the arXiv version (math/9611225) gives the same with
  prescribed `χ_D(p_i) = ε_i` at any finite set of primes. Existence, not density; the
  `p_i` are auxiliary. We have not checked the corollary numbering against the Annals
  text, and it can differ between preprint and published version.
- **Bruinier** (Duke Math. J. 98, 1999): Theorem 1 is the closest structural result —
  if *all* coefficients in one Legendre class `(n/p) = ε` vanish mod `ℓ`, then `λ_p`
  satisfies an explicit congruence. Legendre classes are thus known to organise mod-`ℓ`
  vanishing. But the statement is all-or-nothing, and the hypotheses require
  `(ℓ, Np(p-1)) = 1`, hence `p ≠ ℓ`.
- **Bruinier–Ono** (JNT 99, 2003), **Choi** (2007): "Property A" —
  `#{n ≤ X : a(n) ≡ r} ≫ X` for `r ≡ 0` and `≫ √X/log X` for `r ≢ 0`. Positive density
  of vanishing, no rate, no square-class refinement.
- **Ahlgren–Boylan** (Amer. J. Math. 129, 2007): if two suitable `D` give
  non-vanishing algebraic `L`-values mod `ℓ`, infinitely many do.
- **Bellaïche–Green–Soundararajan** (Res. Math. Sci. 5, 2018):
  `#{n ≤ X : p(n) ≢ 0 (mod ℓ)} ≫ √X/log log X`. The authors note `≫ X` is *expected*
  for `Σp(n)q^{24n-1}` but not known.

The last sets the scale: the best available lower bound for non-vanishing is `√X`;
the observed count is `(12/13)X`. The gap is a square root.

Our §6.2–6.4 measures a *partial* Legendre correlation at every prime, where the
literature treats only the degenerate case of complete vanishing on one class. That
the effect appears at essentially all `q` is consistent with Bruinier's framework
being the right one; what is not available anywhere is a density or a rate.

---

## 8. Questions

1. Is `density{t squarefree : a(t) ≡ 0 (mod 13)} = 1/13`? Even positive density for
   the complement is open. What the data now says is that `E(t)` reaches `+0.0163` in
   the top half-decade at `MMAX = 10¹⁰` and that successive half-decades multiply it by
   0.73, 0.73, 0.77 — no flattening anywhere, so any limit above about `0.016` is
   already excluded. An earlier draft rested part of this question on the two smallest
   windows looking flat; that flatness was an artefact of the window definition (§6.1)
   and is not evidence for a plateau.
2. What is `f_q`? The `q^{-1/2}` prediction is rejected, but no single rational is
   picked out: `3/8`, `2/5` and `5/12` all survive at 723 primes, and the estimate drifts
   from `0.405` to `0.427` on extending to 2260. This is the central mystery, and the
   drift is the part to resolve first.
3. What governs the sign and size of `f_q` across `q`? Ruled out: dependence on
   `λ_q mod 13` alone; proportionality to a Sato–Tate-distributed `λ̃_q/√q`;
   `λ_q ± ε(q)q⁴ ≡ 0` as a predictor of sign.
4. What is the Shimura lift of `η^11` modulo 13? **All levels dividing 288 are now
   excluded** (§3), so it is not a weight-10 newform of level dividing 288 — or `f_q`
   is not a normalised eigenvalue sequence at all. Which of the two is open.
5. Are the `f_q` multiplicatively independent across primes? Untested.

---

## 9. Limitations

1. **Range.** Kernels reach `t ≤ 10⁹` here and `10¹⁰` in the extended run, where
   `E(t) = +0.0163` in the top half-decade with no sign of flattening. Limits above
   about `0.016` are therefore excluded; distinguishing `E(t) → 0` from a limit an order
   of magnitude below that would need several more decades.
2. **Controls are asymmetric.** `ℓ = 13` runs to `10⁹` in `m`; controls to `3·10⁷`
   in `n`. Claims about `ℓ`-specificity are correspondingly weak — and in the one
   place where we extended the range (§6.3), the `ℓ`-specific claim died.
3. **Observation 3.1 is verified, not derived.** It should follow from [FKO, Thm 1.3]
   and [FKO, (1.4)]; we have not written the derivation.
4. **`E(t)`'s exponent depends on which windows are fitted**, by about 0.020 — less
   than an earlier draft claimed, which was measuring a range-dependent quantity (§6.1).
   `N(X)` has no such freedom, so `β` is what we quote.
5. **The excess and its Legendre partition are one phenomenon, not two** (§6.2), so
   any mechanism must explain both together. We have none.

---

## Corrections to draft v3

Applied 2026-08-18. Items 1–5 were identified by the author, 6–10 by cross-checking
every figure against `results/EXPECTED.md`.

1. **§6.1, `E(t)` exponents.** Quoted as `-0.33` and `-0.28`; the fitted values are
   `16.790·t^(-0.3201)` over all ten windows and `10.419·t^(-0.2924)` over the six
   with `t ≥ 10⁶`. The rounded `-0.28` invited the reading that the complementary
   relation `γ = 1 - β` reproduces it. It does not: `1 - β = 0.2759` lies outside both
   fits, and that gap was stated as evidence against a single power law rather than
   suppressed. **Item 17 supersedes this reading**: the gap was an artefact of the
   window definition, and on squarefree kernels `1 - β` lies between the two fits.
2. **§6.2, factorisation table.** The spread was given as 5–15%; with `ddof = 1` it is
   8–29%, and the `ddof` is now stated. The larger error was the table itself — all six
   rows were still mixed-population (`-1.303 / -1.293 / -1.152` at `q = 5`), despite
   this draft's header claiming a single consistent population. Replaced throughout.
   The `q = 13` row is unchanged, as it must be.
3. **§6.3, the standing of `q = 13`.** "Sits mid-pack by `|f_q|`" was wrong: it ranks
   3rd of 723 raw, and 88th percentile after dividing out the `q^{-0.405}` trend. The
   claim is now "high but not an outlier", and both rankings are given, since the raw
   one is confounded by the decay measured in §6.4.
4. **Data section.** Now names `data/delta_q_consistent.json`, its population, and the
   superseded `data/delta_q_mixed_legacy.json`.
5. **§6.1, the `N(X)` table.** Read `3,478 / 18,566 / 96,022`; the correct values are
   `3,561 / 18,700 / 96,022`. Only `X = 10⁹` had been measured on the consistent
   population — the other two were obtained by rescaling the mixed-population values
   by `96022/123563`, which is why exactly one of the three was right. The table now
   carries the underlying counts so each row can be checked as `zeros - n/13`.
6. **§5.3, distribution table.** `1.089 / 1.026 / 1.211` on "3.26 million samples"
   reproduces on no population we can identify. Over all squarefree `m ≤ 10⁹`:
   `1.052 / 1.018 / 1.132`, on 17,640,743 / 17,640,699 / 2,714,016 samples. The
   qualitative claim — non-zero classes uniform to within 1% — survives.
7. **§5.5, kernel-dependence table.** Wrong in both columns; the `[10⁸,10⁹)` entry read
   `+0.0292`, which is `δ₁₃`, not `E`. The correct row is `+0.0317 / +0.0317`, so the
   four-decimal agreement still holds, on different numbers.
8. **§6.3, `f_q` column.** Every entry was one out in the third decimal, from dividing
   `δ_q` by the rounded `E = 0.0317` rather than by the full-precision value.
9. **§3, verification counts.** "11401 cases for each `p ≤ 23`, decreasing to 3175 at
   `p = 47`" cannot be right: `t < 2·10⁵` binds before `t p² ≤ 10⁹` does, so the count
   is uniform in `p` at 7,600, totalling 91200 congruences, with 18833 for the
   three-term relation.
10. **§5.1, number of independent checks.** Seven, now eight — a check that a
    regenerated `results/delta_q.json` reproduces `data/delta_q_consistent.json` was
    added to `test_verify.py`.

11. **§4, the purity row.** Quoted `3.19 / 2.78 / 2.81 / 2.55` with no kernel bound
    named. Section 4 had no script behind it until now; with `src/pgen.c` and
    `src/controls.py` the figures at `t < 4·10⁵` are `3.06 / 2.69 / 2.62 / 2.43`, and
    the earlier values match no bound the script reproduces. The quantity is a
    largest-class share that decays toward `1/ℓ` as the bound grows, so it is
    meaningless without one — the same failure as quoting an `E(t)` exponent without
    saying which windows were fitted. The dichotomy the table asserts is unaffected.

Items 12–16 come from the audit of 2026-08-19, which classified all 101 numeric claims
by whether a script in the repo regenerates them (`paper/claims_audit.tsv`).

**This log is itself checked.** An entry's justification can be withdrawn by a later
entry without the entry noticing — item 1 asserted that the gap between `1 - β` and the
`E(t)` fits was evidence against a power law, and item 17 showed the gap was an artefact
of the window definition, leaving item 1 standing as a live conclusion. So every entry
must state what it superseded, and whenever two entries touch the same section the
earlier must either reference the later or be recorded as reviewed with the reason it
still stands. `paper/check_against_expected.py` enforces this and verifies the numbering
is contiguous; adding an entry forces a re-read of every earlier entry on its section.

12. **§3, congruent primes with different `λ_p`.** The pairs cited were `(5,31)`,
    `(11,37)`, `(17,43)`. But `λ₁₁ = λ₃₇ = 5`: the middle pair is congruent mod 13 with
    *equal* `λ_p`, so it was a counterexample to the sentence citing it. Now cites the
    two valid pairs and notes the third as the converse phenomenon. An independent fit
    of `(λ_p, ε(p))` for all 50 primes `5 ≤ p ≤ 250` confirms the twelve tabulated
    `λ_p`, the character rule at 50/50, and constancy of the ratio on both classes at
    50/50.
13. **§2, range of the direct verification.** Read `n < 30000`; `test_verify.py` checks
    `n ≤ 40000`. The paper understated its own check.
14. **§6.3, the trend-divided median.** Read `1.22`; the value is `1.2146`, so `1.21`.
15. **§6.4, the standing of `2/5`.** "The only simple rational not rejected" was false:
    `3/8` (1.30) and `5/12` (0.23) are equally consistent, and were simply absent from
    the candidate list. With a 95% interval of width 0.102 no candidate inside it could
    have been rejected.
16. **Header, size of the inflation.** Read "about 24%". The factor is
    `E_mixed/E_consistent = 0.037917/0.031656 = 1.198`. The `1.24` came from evaluating
    `(12/13) + (1/13)(E₁₃/E_bulk)` with the wrong stratum weights; the true count
    fractions are 0.9286 and 0.0714.

17. **§6.1, the definition of the E(t) windows.** The table took every `m` whose kernel
    fell in the window rather than the squarefree `m`, which makes it depend on the range
    of the computation: the `[10⁵, 10⁵·⁵)` window reads `B = 1.2592` at `MMAX = 10⁹` and
    `1.6269` at `10¹⁰`. Restated on squarefree kernels, where every window shared by the
    two computations is identical in `n` and in its zero count. This was not a
    transcription error and no invariant caught it — the figure was correctly emitted by
    the script, and only a second computation at a different range revealed that the
    quantity was not the one the section claimed to measure. Three consequences: the
    apparent flatness of the two smallest windows was an artefact, the two window fits
    now diverge by 0.020 rather than 0.028 with the rms log-residual falling from 0.167
    to 0.073, and `1 - β = 0.2759` lies between the two fitted exponents instead of
    outside both. An earlier draft read that last discrepancy as evidence against a power
    law; it was evidence against the definition.

Also corrected without changing a figure: §2's list of primes for which the single-eta
reduction fails omitted `ℓ = 19` and `ℓ = 43`; all nine now appear, each verified to
fail at `i = 1` over 120,000 coefficients. §8's fourth question still described level
288 as untested.

Two claims were *verified* rather than corrected. §5.1's numpy cross-check on
7,692,308 coefficients and §5.2's absence of exact zeros to `m = 92,000,003` both
exceeded what the routine test covers (200,000 coefficients each); both were re-run in
full and both hold. They are now `src/extended_verify.py`.

`python3 paper/check_against_expected.py` asserts that every figure in this paper also
appears in `results/EXPECTED.md`, and that this log carries no withdrawn conclusion.

These seventeen items are of six kinds, listed in Appendix A. The worked examples are
the entries above.

---

## Appendix A. The six failure modes

Worked examples are the corrections log.

1. **Transcription.** Computed once, copied into two documents.
2. **Untested claim.** Asserted beyond any test's range.
3. **Literature at second hand.** Taken from an abstract; no recomputation reaches it.
4. **Unscripted figure inside the checker.** Asserted by the checker, computed by nothing.
5. **Correct computation of the wrong quantity.** Faithfully emitted, but not what the
   section claimed; found only at a second range.
6. **Withdrawn justification left standing.** Overturned by a later entry, still live.

Modes 3, 5 and 6 were invisible to every mechanism in place when they were introduced.
When adding a mechanism, ask what it cannot see.

---

## Appendix B. Extent of AI involvement, and the limits of the verification

The computations, the analysis code, the invariants and much of this text were produced
with an AI assistant (Claude, Anthropic) working interactively with the author across
sessions. The author set the questions, supplied several components and independent
measurements, checked the literature against the original papers, and made the
editorial decisions. Errors that remain are the author's.

**The independence claimed in §5.4 is of algorithms, not of authors.** The pipeline,
`src/independent_check.py`, the invariants and the small-range brute force were written
in the same project by the same model family at different times. Agreement among them
rules out implementation slips — an off-by-one in a sieve, a wrong modular inverse, a
float mean where an integer ratio belongs — because those do not survive being rewritten
a different way. It does not rule out a shared misconception about *which quantity to
compute*, because every implementation inherits it. Agreement between implementations of
common origin is therefore weaker evidence than agreement between implementations
designed independently, and should be read as the former.

This is not a hypothetical reservation. The `E(t)` window definition of §6.1 was emitted
correctly by the pipeline, passed all twelve invariants, reproduced under
reimplementation by different algorithms, and matched a brute-force computation in exact
integers — and was still not the quantity the section claimed, because all four shared
one assumption about which `m` belonged in a window. It fell only to a computation at
ten times the range. That is correction 17, and failure mode 5 of Appendix A; modes 3
and 5 are precisely the classes that same-origin implementations cannot reach.

Three checks in this paper do not have that limitation, having involved a second party:

- the author's independent measurement of `δ_q` over 723 primes, which reproduces
  `data/delta_q_consistent.json` bit-for-bit and is asserted by `test_verify.py`;
- the author's independent fit at `q < 20000`, which replicates the `|Z| ≥ 3` count
  exactly and six of the eight candidate statistics to the digit (§6.4);
- the literature of §7, read by the author against the original papers after two of
  three summaries taken at second hand proved wrong.

A reader should weight §5.4 as evidence against implementation error, which it is, and
should not read it as evidence that the quantities are the intended ones. For that, the
guarantees are the second-party checks above and the agreement of `β` across a decade of
range — not the number of mechanisms that agree.

---

## Data and code

Generation, verification and analysis code, with the expected output of every table
above in `results/EXPECTED.md`: <https://github.com/Jun4M/eta11-mod-13>, archived at
[doi:10.5281/zenodo.22016761](https://doi.org/10.5281/zenodo.22016761) — the concept
DOI, which resolves to the most recent archived version.

The coefficient files are regenerable rather than distributed, being too large to
version:

- `a(m) mod 13` for `m ≤ 10⁹`, `m ≡ 11 (mod 24)`: 41,666,667 bytes, about seven minutes
  on eight cores (`./eta11 1000000000 11`).
- the same for `m ≤ 10¹⁰`: 416,666,667 bytes, about fifty minutes, used for the `β` and
  `E(t)` extensions of §6.1 and the `α` measurement of §6.4
  (`./eta11 10000000000 2`, repeated until it reports COMPLETE).
- `p(n) mod ℓ` for `n ≤ 3·10⁷`, `ℓ ∈ {13,17,19,23,29,31,37,41,43,47}`, about one minute
  (`./pgen 30000000`), for the controls of §4.

Measured values are in the repository:

- `δ_q` and `Z` for 723 primes on the consistent population
  (`data/delta_q_consistent.json`: squarefree kernels with `m ≥ 10⁸` and `13 ∤ m`);
  `f_q = δ_q/E` is derived from these with `E = +0.0317`. The superseded
  mixed-population measurement is retained as `data/delta_q_mixed_legacy.json` so that
  earlier numbers can be traced; it should not be used.
- the same for 2260 primes to `q < 20000` (`results/delta_q_q20000_from1e9.json`) and
  for the `m ≤ 10¹⁰` population (`results/delta_q_1e10.json`), both used in §6.4.
- Generators of the square classes of Corollary 3.2: OEIS A399067.

## References

- S. Ahlgren, M. Boylan, *Central critical values of modular L-functions and
  coefficients of half-integral weight modular forms modulo ℓ*, Amer. J. Math. **129**
  (2007), 429–454.
- A. O. L. Atkin, J. N. O'Brien, *Some properties of p(n) and c(n) modulo powers of
  13*, Trans. Amer. Math. Soc. **126** (1967), 442–459.
- J. Bellaïche, B. Green, K. Soundararajan, *Nonzero coefficients of half-integral
  weight modular forms mod ℓ*, Res. Math. Sci. **5** (2018), Paper 6.
- J. H. Bruinier, *Nonvanishing modulo ℓ of Fourier coefficients of half-integral
  weight modular forms*, Duke Math. J. **98** (1999), 595–611.
- J. H. Bruinier, K. Ono, *Coefficients of half-integral weight modular forms*,
  J. Number Theory **99** (2003), 164–179; corrigendum **104** (2004), 378–379.
- D. Choi, *Distribution of integral Fourier coefficients of a modular form of half
  integral weight modulo primes*, arXiv:0704.0012.
- A. Folsom, Z. A. Kent, K. Ono, *ℓ-adic properties of the partition function*,
  Adv. Math. **229** (2012), 1586–1609. [FKO]
- K. Ono, C. Skinner, *Fourier coefficients of half-integral weight modular forms
  modulo ℓ*, Ann. of Math. **147** (1998), 453–470.
- G. Shimura, *On modular forms of half integral weight*, Ann. of Math. **97** (1973),
  440–481.
- J.-L. Waldspurger, *Sur les coefficients de Fourier des formes modulaires de
  poids demi-entier*, J. Math. Pures Appl. **60** (1981), 375–484.
- OEIS A000041, A010819, A071750, A399067.
