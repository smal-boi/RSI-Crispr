# Handoff — the human arm, after reproducing ρ 0.472

**For:** whoever picks up the human boundary (item 4's original owner).
**Date:** 16 October 2026.
**Code:** `src/sgrna/human_baseline.py`, `src/sgrna/kingdom_transfer.py`.
**Results:** `results/human_baseline.csv`, `human_baseline_knobs.csv`,
`human_additivity.csv`, `kingdom_transfer.csv`, `kingdom_transfer_boot.csv`.

---

## 1. The number reproduces, and it is not about features

`HUMAN_BASELINE_SETTINGS.md` gave 400 selected columns, 1,000 estimators, every
other XGBoost knob at library defaults. Run through the current harness on the
human iRF matrix (17,421 × 6,233), 3 seeds × 5 folds, same selector:

| arm | columns | settings | ρ |
|---|---:|---|---:|
| `harness_default` | 300 | champion | 0.4036 ± 0.0128 |
| `jacky_cols_only` | 400 | champion | 0.4074 ± 0.0148 |
| **`jacky_400_1000`** | **400** | **library defaults, 1,000 trees** | **0.5076 ± 0.0149** |

**ρ 0.472 reproduces and is exceeded.** Nothing about it is historical, and the
caveat in your note is not needed.

The third arm is there because the quoted configuration changes two things at
once. It settles which: **the column count is worth +0.0038 (noise); the
hyperparameters are worth +0.1002.** Had we not split them, the natural reading
would have been "human needs more features", which is the opposite of true.

## 2. Which knob, and why it matters

Each knob moved alone from the 400-column reference at 0.4074:

| change from the champion | ρ | Δ |
|---|---:|---:|
| `max_depth` 4 → 6 | 0.4659 | **+0.0585** |
| `learning_rate` 0.05 → 0.3 | 0.4490 | +0.0416 |
| `n_estimators` 400 → 1,000 | 0.4435 | +0.0361 |
| no row/column subsampling | 0.3999 | −0.0075 |
| all three | 0.5076 | +0.1002 |

Sub-additive, as capacity knobs should be. The biggest single term is **depth**,
which is the interaction-order knob — so the suggestion is that human efficiency
has interaction structure *E. coli*'s lacks. Depth also buys plain capacity, so
that needed a direct test. It has one, and it agrees:

| ridge against boosting, same selected columns | *E. coli* | human |
|---|---:|---:|
| boosting | 0.608 | 0.5076 |
| ridge | 0.577 | 0.3232 |
| **ratio** | **95.6%** | **63.7%** |

**A straight line gets within 0.03 ρ of boosting in *E. coli* and falls 0.18
short in human.** Your "`n_estimators` kept improving, unlike bacteria" was the
same observation from the other end.

So the champion is not arbitrarily mis-tuned — it is tuned for a nearly-linear
signal, and the human signal is not nearly linear. That also **bounds** the
worry: the champion travels fine across the bacterial sets, where the signal
shape is the same, and fails across kingdoms, where it changes.

Caveat on the ratio: the two percentages are not on an identical basis
(*E. coli*'s on 427 decorrelated columns, human's on 400 selected ones), so treat
it as a factor-of-magnitude comparison, not a precise number. **Putting both on
a matched basis is the one cheap thing left here** (§4).

## 3. It cost us a claim — please know this before you quote Part 3

The script behind `results/human_transfer.csv` was never in `src/`, so it was
rebuilt as `kingdom_transfer.py` (6,216 columns shared by name, no change of
representation) and run at both settings, with intervals from 2,000 test-row
resamples instead of from three seeds:

| direction | champion | human-fitted | human-fitted, seeds pooled | 95% CI |
|---|---:|---:|---:|---|
| *E. coli* → human | −0.0462 | −0.0367 | −0.0393 | [−0.0542, −0.0242] |
| human → *E. coli* | −0.0155 | −0.0055 | **−0.0060** | **[−0.0227, +0.0111]** |

The first two columns are means of three single-seed fits; the third averages the
three seeds' predictions before scoring, which is what the interval is
bootstrapped on. The two conventions agree to 0.003 and are kept separate rather
than mixed.

The champion column reproduces the quoted −0.048 and −0.017. The human-fitted
column moves the reverse direction **to zero** — and it is the direction trained
on human, the only one whose hyperparameters were wrong for its training data.

**"Negative in both directions, intervals excluding zero" is withdrawn.** It is
replaced by "no transfer in either direction, measurably negative in one". The
boundary survives; the symmetry does not. Part 3, Finding 5's table, Part 5's
withdrawal list and the README are all updated — the claim to quote now is
*nothing transfers*, not *transfer is negative*.

## 4. What is actually left, in the order I would do it

1. **The matched-basis additivity ratio** (~1 h). Build a decorrelated human
   basis the way Appendix B built *E. coli*'s 427, then ridge against boosting on
   both. This is the only soft number above, and it is the one carrying the
   interesting claim. If it holds, the non-additivity contrast is a result in its
   own right rather than a footnote explaining our hyperparameters.
2. **Re-tune the human arm properly** (~2 h). 0.5076 is a point on a curve you
   already said has no plateau — nobody has searched depth, shrinkage and tree
   count jointly on human. Worth knowing the actual human ceiling before anyone
   writes "the human matrix is less predictable", because at the moment that
   sentence is indistinguishable from "we fitted human badly".
3. **Flank features on human** (half a day, speculative). Family A is the whole
   contribution on the bacterial side and has never been built for human
   coordinates. If human efficiency is interaction-driven, flanking context could
   behave quite differently there. This is a new question, not a loose end — it
   does not block the *E. coli* paper.

**Not worth doing:** re-running the bacterial external arms at per-dataset
hyperparameters. Each is a *difference between arms on one dataset*, so the
mis-specification largely cancels; the absolute figures are lower bounds, and the
report now says so.

## 5. One thing to check if you disagree with any of the above

Everything here is three seeds of five folds on fixed splits, and the knob
decomposition is a one-factor-at-a-time design against a single reference arm —
so it attributes differences, it does not search. If you think the reference arm
is the wrong one (e.g. that 400 columns at *tuned* settings would beat 400 at
champion settings by more than +0.0038), that is a direct disagreement with the
table in §1 and is one run away from being settled:

```bash
python -m sgrna.human_baseline --decompose --arms plus_depth plus_lr
```
