# CLAUDE.md

These are the working rules I impose on any assistant used in this repository.
Responsibility for the contents rests with me.

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
- Report confidence intervals, not point estimates — and where the estimate is not
  settled, say that instead of quoting either. `alpha` is the live example: it reads
  0.405, 0.427 and 0.430 on three populations and its interval does not narrow with
  `m`, so §6.4 names no value and no rational. The older phrasing of this rule, kept
  below for its shape, understated that. `alpha = 0.405` alone is
  misleading; `[0.355, 0.457]` is the result. And where there is no well-defined
  interval — `E(t)`, whose exponent depends on which windows are fitted — say so
  rather than quoting a single exponent.

**The paper and `results/EXPECTED.md` must never differ.** Every figure in
`paper/manuscript.md` also appears in `EXPECTED.md`, and `python3
paper/check_against_expected.py` asserts it, string-exact, for every figure in its
`FIGURES` list. Do not quote the count here: it drifted from 53 to 119 while this line
said 53, which is the failure this very rule exists to prevent. Run the checker after
touching either file. When a number changes, change it in both, and add it to the
`FIGURES` list if it is new. The recurring failure here is a number corrected in one
document and left stale in the other.

**A figure may only enter `results/EXPECTED.md` if a script in `src/` emitted it in the
same run.** Not "could in principle be computed" — emitted, by a script, in the run
being recorded. Transcribing a number into EXPECTED.md and into the paper makes
`check_against_expected.py` assert a fiction: the checker compares two documents and
cannot catch a figure that is wrong in both. That is exactly how the population-mixing
bug survived, and the 2026-08-19 audit found seven more figures in the same state,
introduced while fixing it.

If no script emits a figure, the figure does not go in EXPECTED.md. If that means it has
to leave the paper until a script exists, **it leaves.** `paper/claims_audit.tsv` records,
for every numeric claim, which script produces it; a new figure needs a row there.

**The invariants are the checking layer, not the checker.** `src/invariants.py` is part
of the routine path, not something run on suspicion — it is the only layer that can fail
on a number that is consistently wrong everywhere. Add an invariant for every *relation*
the paper asserts, not for every number it quotes: the assertion that `δ_q` and `E` share
a population would have caught the eleven corrections, where no amount of agreement
between documents could. Routine path, in order:

```bash
python3 src/test_verify.py a13.bin         # must print "all checks passed"
python3 src/analyze.py a13.bin
python3 src/invariants.py a13.bin          # must print "all invariants hold"
python3 src/fit_alpha.py results/delta_q.json
python3 src/factorisation_check.py a13.bin
python3 src/paper_figures.py a13.bin
python3 paper/check_against_expected.py    # last, and the weakest of the checks
```

Run `src/independent_check.py` when a script changes: it reimplements each quantity by a
different algorithm (kernels by factorisation, Legendre by Euler's criterion, `α` by
binned moments), which a refactor of the same algorithm cannot substitute for.

**`analyze.py` and `invariants.py` are memory-lean by necessity.** They keep only the
uint32 square part and two bool arrays full-length and form `m`, `t` and the Legendre
class in 2^24 chunks. The earlier int64 version needed 21 GB at `MMAX = 10^10`. Peak is
about 2.2 GB at `10^9` and 6 GB at `10^10`. There is one implementation of the analysis;
keep it that way — two is how the population bug lasted as long as it did.

**A quantity claimed to depend only on t must be checked at two ranges.** The E(t)
window table was correctly emitted by `analyze.py`, recorded in EXPECTED.md and quoted
in the paper, and was still measuring the wrong thing: it took every `m` whose kernel
fell in the window, so it drifted with MMAX (the `[1e5,1e5.5)` window read 1.2592 at
`10^9` and 1.6269 at `10^10`). No invariant, no reimplementation and no brute force at
small range could catch that — only a second computation at a different range did.
Invariant 12 now enforces the property that makes the quantity range-free: every `m`
counted in a window must lie in that window, checked against what `analyze.py` wrote.

**The corrections log is a document that can go stale, and it is checked.** An entry's
justification can be withdrawn by a later entry without the entry noticing: item 1 stood
as a live conclusion for a day after item 17 overturned its reasoning. Every entry states
what it superseded; whenever two entries touch the same section, the earlier must either
reference the later by number or be recorded in `REVIEWED` in
`paper/check_against_expected.py` with the reason it still stands. Adding an entry
therefore forces a re-read of every earlier entry on its section, and the re-read has to
be written down. Do not clear a pair by adding it to `REVIEWED` without actually
re-reading the earlier entry.

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

### 1. ~~Level 288 Shimura lift~~ — DONE 2026-08-19, clean negative

All 17 levels dividing 288 searched, 142 form–embedding pairs. Max |correlation| 0.294
at level 288, 0.377 overall, both at or below the null median for their candidate
counts. `mfinit([288,10],0)` fits in a 4 GB stack (`gp -s 4000000000`); PARI/GP 2.17.4
via `brew install pari`. Do not re-open this without a new idea about what `f_q` is —
re-running the same correlation will not change the answer.

When comparing against `f_q`, note that `results/fq_for_pari.txt` must be regenerated
from `results/delta_q.json`; the version shipped before 2026-08-19 predated the
population fix and was not a rescaling of the corrected values.

### 2. ~~Extend to `m ≤ 10^10`~~ — DONE 2026-08-19

`./eta11 10000000000 2`, repeated until COMPLETE; 0.83 GB resident, about fifty minutes.
`analyze.py` needs 2.6 GB at this range and 7.1 GB peak through the delta_q loop, which
took 7035 s. Guard: `state*.bin` is tied to a specific `MMAX` — delete it when changing
`MMAX`, and note it is gitignored by glob because a rename once put 397 MB into history.

**The stated expectation was wrong, and that is the useful part.** This task said the
extra decade "should narrow the `alpha` interval enough to test `2/5` against
neighbours". It did not: the 95% width went 0.102 -> 0.100 while the measurement error
fell by a factor 3.2. `alpha`'s precision is set by the number of auxiliary primes, not
by `m` — measured directly, widths 0.202, 0.148, 0.100, 0.075 at 180, 361, 723, 2260
primes, against 0.100 -> 0.082 for a tenfold cut in the SE. So no further extension of
`m` will settle `alpha`, and §6.4 is written around that.

What the decade did buy: `beta` confirmed at 0.7242 against 0.7241, moving by 0.0001,
with `3/4` rejected more strongly; two more `E(t)` windows; and the discovery that the
old `E(t)` window definition was MMAX-dependent, which only a second range could show.

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
