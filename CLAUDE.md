# CLAUDE.md

Instructions for Claude Code working in this repository.

## What this is

A computational number theory project measuring how often the Fourier coefficients
of `eta^11` vanish modulo 13, and how that frequency correlates with quadratic
residue classes. Read `README.md` first — it states the two measured laws and the
open questions.

## Ground rules

**Correctness before anything else.** Every claim in this repo is a numerical
measurement that someone may cite. A wrong number is worse than no number.

- Never report a result from a single computation. Cross-check with an independent
  method before believing it. `src/test_verify.py` is the model: it recomputes
  things by different code paths and with exact arithmetic.
- Run `python3 src/test_verify.py a13.bin` after regenerating `a13.bin`, always.
  If it does not print `all checks passed`, stop and diagnose.
- When fitting exponents, model the measurement noise explicitly. The binomial
  standard error of `delta_q` is constant in `q` (about 0.0012 at `m ≤ 10⁹`), but
  after multiplying by `sqrt(q)` it grows — this artefact has already produced two
  wrong conclusions in this project. `src/fit_alpha.py` handles it correctly; copy
  that approach.
- Report confidence intervals, not point estimates. `alpha = 0.405` alone is
  misleading; `[0.355, 0.457]` is the result. And where there is no well-defined
  interval — `E(t)`, whose exponent depends on which windows are fitted — say so
  rather than quoting a single exponent.

**The paper and `results/EXPECTED.md` must never differ.** Every figure in
`paper/manuscript.md` also appears in `EXPECTED.md`, and `python3
paper/check_against_expected.py` asserts it — 53 figures, string-exact. Run it after
touching either file. When a number changes, change it in both, and add it to the
`FIGURES` list if it is new. The recurring failure here is a number corrected in one
document and left stale in the other.

**Do not overclaim.** The history of this project is a sequence of claims that
looked clean at small sample size and dissolved at larger sample size:
`q^(-1/2)` scaling, Sato–Tate distribution, and a supposed special role for
`q = 13` all died this way. If a pattern appears, the next step is to test it on
more data, not to write it up.

**Check the literature before claiming novelty.** A large part of this project
turned out to be a rediscovery of Folsom–Kent–Ono Theorem 1.3. When something
looks like a new structural result, assume it is known until checked.

## Environment

- C compiler with OpenMP for `src/eta11.c`.
- Python 3 with numpy.
- PARI/GP ≥ 2.11 for `src/shimura288.gp` (`apt install pari-gp` or `brew install pari`).
  GP reads script files **line by line** — every statement must be on one line, and
  `#` cannot be used as the length operator in a script (use `length()`).

## Priority tasks

### 1. Level 288 Shimura lift (highest value)

`src/shimura288.gp` tests whether a weight-10 newform has normalised eigenvalues
`a_q/q^(9/2)` correlating with the measured `f_q`. Levels up to 144 are done and all
correlate at noise level. Level 288 needs a large stack:

```bash
gp -s 4000000000 -q src/shimura288.gp
```

`mfinit([288,10],0)` is the expensive call. If it exhausts memory, try
`default(parisizemax, ...)`, or compute the newspace in pieces with `mfsplit`.

A correlation above 0.9 identifies the lift and probably settles the mechanism.
Below about 0.5 is noise. **Report the actual maximum, whatever it is** — a clean
negative result is valuable here.

### 2. Extend to `m ≤ 10^10`

```bash
./eta11 10000000000 2     # repeat until COMPLETE; checkpoints in state.bin
```

0.83 GB resident. Then rerun `analyze.py` and `fit_alpha.py`. This should narrow
the `alpha` interval enough to test `2/5` against neighbours, and add one more
decade to `E(t)`.

Guard: `state.bin` is tied to a specific `MMAX`. Delete it when changing `MMAX`.

### 3. Mechanism hunting

Open-ended. Things already tried and ruled out, do not repeat:

- `f_q` as a function of `lambda_q mod 13` — no (primes with equal `lambda_q` have
  different `f_q`).
- `f_q ∝ lambda~_q / sqrt(q)` with Sato–Tate `lambda~` — the exponent `1/2` is
  rejected.
- weight-10 eta quotients on divisors of 24 — exhaustive search, no match.
- level-1 forms (`Delta`, `Delta·E_4` … `Delta·E_14`) with all twists `p^j` — noise level.
- `lambda_q ± eps(q) q^4 ≡ 0 (mod 13)` as a predictor of the sign of `f_q` — no.

Untried and worth a look: whether the `f_q` are multiplicatively independent across
primes (test joint conditioning on `(t|q1)` and `(t|q2)`); whether `E(t)` and `f_q`
have a common source in the mod-13 Galois representation; whether the excess
concentrates on `t` with particular factorisation shapes (a partial signal was seen:
`t` with no prime factor below 31623 have roughly half the excess).

## Conventions

- Coefficient files: `a13.bin`, one `signed char` per index `i`, exponent `m = 24i + 11`.
- Results go in `results/` as JSON, so scripts can be chained.
- `data/delta_q_consistent.json` is the reference measurement: 723 primes, squarefree
  kernels with `m >= 1e8` and `13 nmid m`. `test_verify.py` checks that a regenerated
  `results/delta_q.json` reproduces it. Do not overwrite either.
- `data/delta_q_mixed_legacy.json` (formerly `delta_q_measured.json`) was measured on a
  population that included `13 | m` while its normaliser did not. Retained only so that
  earlier numbers can be traced — do not use it.
- New runs write to `results/`.
