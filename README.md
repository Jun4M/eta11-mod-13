[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22016761.svg)](https://doi.org/10.5281/zenodo.22016761)

# p13 — vanishing of `eta^11` coefficients modulo 13

Computational study of how often `13 | p(n)` along the Ramanujan branch, and how
that frequency correlates with quadratic residue classes.

Everything here is **measurement**, not proof. Two laws have been measured; the
mechanism behind either is unknown.

---

## What is being computed

For `n ≡ 6 (mod 13)` put `m = (24n − 1)/13`. Then

```
p(n) ≡ 11 · a(m)   (mod 13),     Σ a(m) q^m = ∏_{k≥1} (1 − q^k)^11
```

(this is Folsom–Kent–Ono §6; the exponents `m` are exactly the integers `≡ 11 mod 24`).
So all 13-adic statistics of `p(n)` on that branch are statistics of a single
half-integral weight form. Write `m = t·s²` with `t` squarefree — the **kernel**.
By the Hecke structure, everything depends only on `t`.

## The two measured laws

With `E(t) := 13·P(a(t) ≡ 0) − 1` over squarefree kernels, and
`δ_q := 13·[P(a(t)≡0 | (t|q)=+1) − P(a(t)≡0 | (t|q)=−1)]`, both measured over
squarefree `t` with `13 ∤ t`:

**Factorisation.** Verified across three windows spanning `t ∈ [10⁶, 10⁹)`, over
which `E` itself falls by a factor 3.8, with 8–22% spread:

```
P( a(t) ≡ 0 mod 13 | (t|q) = ε )  =  (1/13) · [ 1 + E(t) · (1 + f_q·ε/2) ]
```

`E` depends only on `t`, `f_q := δ_q/E` only on `q`.

The `13 ∤ t` restriction is not cosmetic. The `13 | t` stratum has `E = +0.1193`
against `+0.0317` for the bulk on the top decade, and `(t|13) = 0` is structurally
excluded for `q = 13`, so `δ_q` and `E` must be measured on the same `13 ∤ t`
population or `f_q` is not a well-defined ratio.

**Decay rates** (measured to `m ≤ 10⁹`, i.e. 41,666,667 coefficients):

| quantity | fit | ruled out |
|---|---|---|
| `E(t)` | no single power law — see below; cumulative excess count `N(X) = 0.0296·X^0.7241` | `3/4` for the count exponent |
| `f_q` | `q^(−0.405)`, 95% CI `[0.355, 0.457]`, 723 primes to `q = 5483` | `1/2`, `1/3`, `1/4`, `2/3`, `3/4` — **`2/5` fits exactly** |

**`E(t)` is not well described by a single power law.** Fitting all ten half-decade
windows gives `E(t) = 16.79·t^(−0.3201)`; fitting only the six with `t ≥ 10⁶` gives
`10.42·t^(−0.2924)`. The two smallest windows are essentially flat (`E = 0.7614`,
`0.7460`) and only the large-`t` windows decay, so any single exponent is a statement
about which windows were chosen. **The stable quantity to quote is the cumulative
count exponent `β = 0.72`** (`3/4` rejected: rms log-residual 0.063 against 0.020 for
the free fit) — it is a fit to a cumulative count rather than to windowed
differences. Unlike `α`, `E(t)` has no well-defined confidence interval.

`q = 13` is not an outlier, but it is not mid-pack either: by raw `|f_q|` it ranks
3rd of 723, and once the `q^(−0.405)` trend is divided out it sits at the 88th
percentile (`|f_q|·q^0.405 = 2.61` against a median of 1.22 and a maximum of 5.32).
123 of 723 primes have `|Z| ≥ 3` (chance would give 2).

Both effects appear to → 0, so the data is consistent with
`density{t squarefree : a(t) ≡ 0 (mod 13)} = 1/13`, i.e. equidistribution, with a
measured rate. Nothing in the literature gives a density here: the best known
lower bound for non-vanishing is `≫ √X/log log X` (Bellaïche–Green–Soundararajan
2018) against an observed `(12/13)X`.

---

## Open questions, in order of value

1. **What is `f_q`?** The `q^(−1/2)` prediction from a naive Euler-factor /
   Sato–Tate argument is rejected at 11.6 units of `−2log L`. `2/5` fits, with no
   explanation. This is the central mystery.
2. ~~**Identify the Shimura lift.**~~ **Closed, negative.** All 17 levels dividing 288
   have now been searched — 142 form–embedding pairs of weight-10 newforms. Maximum
   |correlation| with `f_q` is 0.294 at level 288 and 0.377 overall, against a null
   median of 0.313 and 0.358 for those candidate counts. The lift is not a newform of
   level dividing 288, or `f_q` is not a normalised eigenvalue sequence.
   See `src/shimura288.gp` and `results/EXPECTED.md`.
3. **Push to `m ≤ 10¹⁰`.** Would tighten both exponents. Sample size for `f_q`
   scales as `n·E²  ∝ T^0.44`, so a decade in `m` buys a factor 1.7 in precision.
4. **Is `E(t) → 0`?** Monotone decrease over five decades with no flattening, but a
   small positive limit is not excluded.

---

## Layout

```
src/eta11.c          coefficient generator (OpenMP, checkpointed)
src/test_verify.py   independent verification — run this first
src/analyze.py       produces E(t), the excess count, and delta_q
src/fit_alpha.py     maximum-likelihood fit of the delta_q exponent
src/factorisation_check.py   f_q across three windows in t
src/extended_verify.py       wide-range numpy and 2^31-1 checks (~12 min, one-off)
src/pgen.c                   p(n) mod l for the section-4 controls
src/controls.py              the three tables of section 4
src/invariants.py            structural assertions that must reconcile
src/independent_check.py     reimplementation by different algorithms
src/paper_figures.py         every paper figure no other script emits
src/shimura288.gp    PARI/GP search for the Shimura lift
paper/manuscript.md          the write-up
paper/check_against_expected.py  asserts the paper and EXPECTED.md agree
data/delta_q_consistent.json    delta_q, Z for 723 primes — the reference measurement
data/delta_q_mixed_legacy.json  superseded; provenance only, do not use
results/delta_q.json         delta_q, Z, f_q regenerated by analyze.py
results/fq_for_pari.txt      input for the PARI script
```

`data/delta_q_consistent.json` is the reference measurement: 723 primes, squarefree
kernels with `m ≥ 10⁸` and `13 ∤ m`. `data/delta_q_mixed_legacy.json` was measured on a
population that included `13 | m` while its normaliser did not; it is retained only so
that earlier numbers can be traced, and should not be used.

## Running

```bash
gcc -O3 -march=native -fopenmp -o eta11 src/eta11.c

# 10^9 : ~0.1 GB, a few minutes on 8 cores. Rerun until it prints COMPLETE.
./eta11 1000000000 11

python3 src/test_verify.py a13.bin      # must print "all checks passed"
python3 src/analyze.py a13.bin
python3 src/fit_alpha.py results/delta_q.json
python3 src/factorisation_check.py a13.bin
python3 src/invariants.py a13.bin           # must print "all invariants hold"
                                            # 11 relation assertions; part of the
                                            # routine path, not run on suspicion
python3 src/independent_check.py a13.bin
python3 src/paper_figures.py a13.bin
python3 paper/check_against_expected.py     # paper and EXPECTED.md must agree
```

For the section-4 controls:

```bash
gcc -O3 -march=native -funroll-loops -o pgen src/pgen.c
./pgen 30000000            # 61 s, 240 MB, writes res_<l>.bin
python3 src/controls.py    # 1.0 s
```

`controls.py` prints the kernel bound it used. The purity percentages for `ℓ ≥ 37`
move with that bound — never quote one without it.

On macOS `gcc` is the Apple clang shim and does not accept `-fopenmp`; use
`clang -O3 -march=native -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include
-L$(brew --prefix libomp)/lib -lomp`, with a libomp at least as new as clang's major
version.

For `m ≤ 10¹⁰`: 0.83 GB resident, roughly 11 × (cost of one pass). Use the
checkpoint — `./eta11 10000000000 2` repeatedly — so a crash costs at most one pass.

For the Shimura lift:

```bash
gp -s 4000000000 -q src/shimura288.gp
```

## Verification

`src/test_verify.py` is deliberately independent of the generator:

- exact big-integer `p(n)` for `n ≤ 40000` (218 digits at the top), checked against
  the reduction `p(n) ≡ 11 a(m)`;
- `eta^11 mod 13` recomputed in numpy by a different code path;
- the same series mod `2³¹ − 1` to confirm there are **no exact zeros**, so the
  vanishing is genuinely 13-adic;
- the Hecke law `a(t p²) ≡ (λ_p + ε(p)(t|p)p⁴)·a(t)` and the three-term recursion
  `a(t p⁴) ≡ λ_p a(t p²) − p⁹ a(t)`, ~110,000 congruences;
- the square class of `t = 155` vanishing entirely (847 coefficients);
- and, once `analyze.py` has run, that `results/delta_q.json` reproduces the reference
  `data/delta_q_consistent.json` on all 723 primes (skipped if either is absent).

Expected output is in `results/EXPECTED.md`.

## Licence

Licensed under the MIT License, in `LICENSE`. This covers the code, the data files, and
the manuscript in `paper/`. The licence text is identical to the one in `erdos126`.

## Background

- Folsom, Kent, Ono, *ℓ-adic properties of the partition function*, Adv. Math. 229 (2012).
  Theorem 1.3 (eigenform congruence for `5 ≤ ℓ ≤ 31`) and eq. (1.4) (the Hecke
  operator) are what the multiplicative law verifies.
- Bruinier, Duke Math. J. 98 (1999); Bruinier–Ono, JNT 99 (2003); Choi,
  arXiv:0704.0012 — Legendre classes are known to organise mod-ℓ vanishing, but all
  results are existential or lower bounds. Bruinier's Theorem 1 requires `p ≠ ℓ`.
- Bellaïche, Green, Soundararajan, Res. Math. Sci. 5 (2018) — the `√X` bound.

## Status

No new theorems. The Hecke-law verification reproduces known results explicitly.
The two decay laws and the factorisation are, as far as we can tell, not in the
literature — but they are measurements, and the mechanism is open.
