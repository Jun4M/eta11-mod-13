#!/usr/bin/env python3
"""
Enforces the rule that paper/manuscript.md and results/EXPECTED.md never differ.

Every figure below must appear verbatim in BOTH files. This is a string check, not a
recomputation -- EXPECTED.md is itself checked against the code by test_verify.py and
by regenerating results/. What this catches is the failure mode this project keeps
hitting: a number corrected in one document and left stale in the other.

Run:  python3 paper/check_against_expected.py
"""
import sys, os

PAPER = "paper/manuscript.md"
EXPECTED = "results/EXPECTED.md"

FIGURES = [
    ("coefficient count",        "41,666,667"),
    ("E first window",           "1.3361"),
    ("E last window",            "1.0290"),
    ("all-m MMAX dependence",    "1.6269"),
    ("superseded all-m fit",     "16.790"),
    ("1e10 E fit, all windows",  "-0.2633"),
    ("1e10 E fit, t >= 1e6",     "-0.2805"),
    ("E at top half-decade, 1e10","+0.0163"),
    ("E at 1e9-1e9.5 window",    "+0.0213"),
    ("excluded limit for E(t)",  "0.016"),
    ("E(t) fit, all windows",    "6.773"),
    ("E(t) fit, all windows",    "-0.2685"),
    ("E(t) fit, t >= 1e6",       "9.680"),
    ("E(t) fit, t >= 1e6",       "-0.2886"),
    ("N(X) prefactor",           "0.0296"),
    ("N(X) exponent",            "0.7241"),
    ("N(X) exponent at 1e10",    "0.7242"),
    ("N(3.16e9)",                "220,868"),
    ("N(1e10)",                  "523,560"),
    ("N(1e10) population",       "352,814,785"),
    ("N(3.16e9) population",     "111,569,767"),
    ("N(1e10) zeros",            "27,663,159"),
    ("N(3.16e9) zeros",          "8,803,158"),
    ("rms at 3/4, 1e10",         "0.079"),
    ("N(1e7)",                   "3,561"),
    ("N(1e8)",                   "18,700"),
    ("N(1e9)",                   "96,022"),
    ("N(X) local slopes",        "0.70"),
    ("N(X) rms free fit",        "0.020"),
    ("N(X) rms beta=3/4",        "0.063"),
    ("top-decade E",             "+0.0317"),
    ("top-decade n",             "31,753,325"),
    ("13|m stratum E",           "+0.1193"),
    ("factorisation E w1",       "0.1199"),
    ("factorisation E w2",       "0.0620"),
    ("factorisation mean q=5",   "-1.034"),
    ("factorisation mean q=7",   "+0.531"),
    ("factorisation mean q=11",  "+1.006"),
    ("factorisation mean q=13",  "+1.012"),
    ("factorisation mean q=23",  "-0.576"),
    ("factorisation mean q=43",  "-0.484"),
    ("factorisation spread",     "8–29%"),
    ("f_q q=5",                  "-0.935"),
    ("f_q q=7",                  "+0.469"),
    ("f_q q=11",                 "+0.955"),
    ("f_q q=13",                 "+0.924"),
    ("f_q q=17",                 "+0.356"),
    ("f_q q=23",                 "-0.497"),
    ("f_q q=43",                 "-0.419"),
    ("primes with |Z|>=3",       "123"),
    ("number of primes",         "723"),
    ("largest prime",            "5483"),
    ("alpha",                    "0.405"),
    ("sigma",                    "0.0438"),
    ("measurement SE",           "0.00125"),
    ("68% interval",             "[0.380, 0.430]"),
    ("95% interval",             "[0.355, 0.457]"),
    ("d(-2logL) 1/4",            "45.88"),
    ("d(-2logL) 1/3",            "8.31"),
    ("d(-2logL) 2/5",            "0.02"),
    ("d(-2logL) 1/2",            "11.64"),
    ("d(-2logL) 2/3",            "69.60"),
    ("d(-2logL) 3/4",            "109.37"),
    ("Hecke congruences",        "91200"),
    ("three-term congruences",   "18833"),
    ("t=155 class size",         "847"),
    ("numpy cross-check",        "7,692,308"),
    ("no exact zeros bound",     "92,000,003"),
    ("N(1e7) population",        "352,805"),
    ("N(1e8) population",        "3,528,117"),
    ("N(1e9) population",        "35,281,442"),
    ("N(1e7) zeros",             "30,700"),
    ("N(1e8) zeros",             "290,094"),
    ("N(1e9) zeros",             "2,809,979"),
    ("purity l=37",              "3.06"),
    ("purity l=41",              "2.69"),
    ("purity l=43",              "2.62"),
    ("purity l=47",              "2.43"),
    ("1/l at l=37",              "2.70"),
    ("1/l at l=47",              "2.13"),
    ("purity drift, TMAX 5e4",   "3.53"),
    ("purity drift, TMAX 6e5",   "3.02"),
    ("withdrawn purity values",  "3.19"),
    ("Shimura max, level 288",   "0.294"),
    ("Shimura max, all levels",  "0.377"),
    ("null median, 45 cands",    "0.313"),
    ("null median, 142 cands",   "0.359"),
    ("correlation SE at n=60",   "0.130"),
    ("delta_q inflation factor",  "1.198"),
    ("trend-divided median",     "1.21"),
    ("q=13 trend-divided value", "2.61"),
    ("q=13 percentile",          "88th"),
    ("max trend-divided",        "5.32"),
    ("argmax prime",             "2971"),
    ("kernels t<1e5 vanishing",  "426"),
    ("nonzero class deviation",  "0.0040"),
    ("d(-2logL) 3/8",            "1.30"),
    ("d(-2logL) 5/12",           "0.23"),
    ("q<20000 alpha",             "0.427"),
    ("q<20000 primes",           "2260"),
    ("q<20000 |Z|>=3 count",     "151"),
    ("1e10 alpha",               "0.430"),
    ("1e10 95% interval",        "[0.382, 0.482]"),
    ("1e10 SE",                  "0.00039"),
    ("1e10 sigma",               "0.0225"),
    ("1e10 |Z|>=3 count",        "205"),
    ("1e10 top-decade E",        "+0.0175"),
    ("1e10 13|m stratum E",      "+0.0578"),
    ("1e10 d 1/2",               "6.74"),
    ("1e10 d 3/8",               "5.30"),
    ("1e10 d 2/5",               "1.55"),
    ("1e10 d 5/12",              "0.32"),
    ("width at 723 primes",      "0.102"),
    ("width at 2260 primes",     "0.075"),
    ("width, SE x 0.1",          "0.082"),
    ("width at 180 primes",      "0.202"),
    ("width at 361 primes",      "0.148"),
    ("moment estimator",         "0.420"),
    ("form-embedding pairs",     "142"),
]

# ---------------------------------------------------------------- corrections log
# The log is not append-only-and-therefore-safe. An entry's justification can be
# withdrawn by a later entry without the entry noticing: item 1 asserted that the
# gap between 1 - beta and the E(t) fits was evidence against a power law, and
# item 17 showed the gap was an artefact of the window definition. Item 1 stood as
# a live conclusion for a day afterwards because nothing was pointed at the log.
#
# Rule enforced below: whenever two entries touch the same section, the earlier one
# must either reference the later by number, or the pair must be recorded here with
# the reason the earlier entry still stands. A new entry therefore forces a re-read
# of every earlier entry on its section, and the re-read has to be written down.
REVIEWED = {
    (1, 5):   "item 5 corrects the N(X) values and so beta; item 1's '1 - beta = 0.2759' "
              "already uses the corrected beta = 0.7241. Item 1 stands, dependent on item 5.",
    (5, 17):  "item 17 changes the E(t) window population only. N(X) was always measured "
              "on squarefree m with 13 nmid m, so item 5 is untouched.",
    (3, 8):   "item 8 shifts f_q by one in the third decimal, uniformly; ranks and "
              "percentiles are invariant under that. Item 3's 3rd-of-723 and 88th "
              "percentile stand, computed from the post-item-8 values.",
    (3, 14):  "item 14 corrects the trend-divided median, which item 3 does not quote. "
              "Item 3 stands.",
    (8, 14):  "item 14's median is computed from the f_q that item 8 corrects, so 14 "
              "depends on 8 rather than overturning it.",
    (9, 12):  "item 12 fits lambda_p and eps independently over the same kernel bound and "
              "confirms the twelve tabulated lambda_p; it does not touch the congruence "
              "counts of item 9.",
}

def check_corrections(paper):
    import re
    from collections import defaultdict
    try:
        start = paper.index("## Corrections to draft v3")
    except ValueError:
        print("[corrections] log section not found -- skipping"); return []
    # end at the next top-level heading, so that numbered lists in later sections
    # (Appendix A, for one) are not parsed as correction entries
    nxt = paper.find("\n## ", start + 1)
    log = paper[start:nxt if nxt != -1 else len(paper)]
    ent = re.findall(r'^\s*(\d+)\.\s+\*\*(.+?)\.\*\*', log, re.M)
    bad = []
    nums = [int(n) for n, _ in ent]
    if nums != list(range(1, len(nums)+1)):
        bad.append(f"entry numbering is not contiguous from 1: {nums}")
    sec = defaultdict(list)
    for num, head in ent:
        # section key is the text before the first comma, minus a leading section sign:
        # "§6.1, ..." -> "6.1"; "Header, ..." -> "Header"; "Data section" -> "Data section"
        key = head.split(',')[0].strip().lstrip('§').strip()
        if not key:
            bad.append(f"item {num} does not declare a section: {head[:40]!r}"); continue
        sec[key].append(int(num))
    for key, items in sec.items():
        items = sorted(items)
        for i in range(len(items)):
            for j in range(i+1, len(items)):
                a, b = items[i], items[j]
                body = re.search(rf'^\s*{a}\.\s+\*\*.*?(?=^\s*\d+\.\s+\*\*|\Z)',
                                 log, re.M | re.S)
                linked = bool(body and re.search(rf'\bitem {b}\b', body.group(0), re.I))
                if not linked and (a, b) not in REVIEWED:
                    bad.append(f"items {a} and {b} both touch section {key}: re-read item {a} "
                               f"and either reference item {b} in it, or record the pair in "
                               f"REVIEWED with the reason item {a} still stands")
    return bad

def main():
    for f in (PAPER, EXPECTED):
        if not os.path.exists(f):
            print(f"FAIL: {f} not found"); return 1
    paper = open(PAPER).read()
    expected = open(EXPECTED).read()
    bad = []
    for label, fig in FIGURES:
        in_p, in_e = fig in paper, fig in expected
        if not (in_p and in_e):
            where = "paper" if in_e else ("EXPECTED" if in_p else "both")
            bad.append((label, fig, where))
    for label, fig, where in bad:
        print(f"[MISSING from {where:8}] {label}: {fig}")
    cbad = check_corrections(paper)
    for msg in cbad:
        print(f"[corrections] {msg}")
    print()
    if bad or cbad:
        if bad:
            print(f"{len(bad)} of {len(FIGURES)} figures do not agree between "
                  f"{PAPER} and {EXPECTED}")
        if cbad:
            print(f"{len(cbad)} unreviewed overlap(s) in the corrections log")
        return 1
    print(f"all {len(FIGURES)} figures agree between {PAPER} and {EXPECTED}")
    print(f"corrections log: numbering contiguous, "
          f"{len(REVIEWED)} same-section overlaps reviewed, all supersessions recorded")
    return 0

sys.exit(main())
