# Results — claim to evidence

This file is an **index**, not a narrative. For every claim the progress report
makes, it names the file that backs it and the command that produced it, so any
number in the report can be checked without reading the report.

The narrative lives in [`../docs/PROGRESS_REPORT.md`](../docs/PROGRESS_REPORT.md);
the chronology, including claims later withdrawn, lives in
[`../docs/FINDINGS_LOG.md`](../docs/FINDINGS_LOG.md).

**What is indexed:** results a claim rests on. That includes nulls when the null
*is* the claim — the dCas9 screen, the periodicity family, the shift control —
because a ruled-out mechanism is a result.

**What is not indexed:** exploratory sweeps, superseded runs, tangential
questions, per-fold dumps (`folds_*.csv`), per-run settings (`*_settings.json`)
and prediction dumps (`preds_*.csv`). Those files are still here and still
readable; they are just not load-bearing. The long-form commentary written up to
10 October 2026 is kept verbatim in
[`RESULTS_archive_2026-10-10.md`](RESULTS_archive_2026-10-10.md).

**Which cohort a number belongs to.** The project pivoted on 13 October 2026 from
Noshay's published 13,880-row matrix and `cut.score` to crisprHAL 2's curated
33,567-guide release and its ALDEx2 score. Both frames appear below and are
marked **[P]** published or **[C]** curated. They are not interchangeable.

---

## The one large effect

| claim | number | frame | evidence | command |
|---|---|---|---|---|
| Flanking DNA composition is the largest effect measured here | +0.164 ρ | [C] | `ablation_curated.csv` | `python -m sgrna.curated_ablation --all` |
| …and +0.080 on the published matrix, so the effect is not an artefact of the pivot | +0.080 ρ | [P] | `ablation_main.csv` | `python -m sgrna.run_ablation` |
| It survives holding out whole genes, losing only the grouping tax | 6% tax | [C] | `ablation_curated_gene.csv`, `grouped_retention.csv` | `python -m sgrna.curated_ablation --all --grouping gene`; `python -m sgrna.rederive --only grouped_retention` |
| It is spatially local: a texture-preserving circular shift destroys it | 8% at 1 kb, 0 at 10 kb | [P] | `shift_null.csv` | `python -m sgrna.shiftnull` |
| It is a smooth gradient, not a motif: k-mer order buys nothing | tetramers −0.0019 | [C] | `curated_ladder.csv`, `kmer_order_noselect.csv` | `python -m sgrna.curated_ladder`; `python -m sgrna.rederive --only kmer_order_noselect` |
| The window widths are not where the gain lives | +0.1634 vs +0.1638 | [C] | `window_robustness.csv` | `python -m sgrna.window_robustness` |
| It is not target GC in disguise: conditioning raises it | +0.044 → +0.082 | [P] | `gc_strata_corr.csv`, `gc_strata_ablation.csv` | `python -m sgrna.gcstrata` |
| It is not gene ranking: it holds inside genes | +0.061 within-gene | [P] | `ablation_grouped.csv` | `python -m sgrna.run_ablation --grouping gene` |
| It is not read depth or mappability | +0.078 → +0.051 | [P] | `ablation_abund_split.csv` (`abund_local_vs_mappability.csv` †) | `python -m sgrna.run_ablation --tag abund_split` |

## The representation

| claim | number | frame | evidence | command |
|---|---|---|---|---|
| 99% of the 6,232 published columns regenerate from a 20-mer | 5,853 of 5,887 `V` columns decoded | [P] | `portability_audit.csv` | `python -m sgrna.expand --build` |
| 412 columns computable from any guide match all 6,232 | 0.00003 ρ | [P] | `portability.csv`, `portability_364.csv` | `python -m sgrna.portability`; `python -m sgrna.rederive --only portability_364` |
| …but no reduced set is non-inferior on the curated cohort | −0.0118 / −0.0142 | [C] | `additivity_reduced.csv` | `python -m sgrna.rederive --only additivity_reduced` |
| A deterministic function of the guide adds no information, but does add performance | +0.011 to +0.019 | [C] | `no_information.csv`, `base_enrichment.csv` | `python -m sgrna.no_information`, `python -m sgrna.base_enrichment` |
| Noshay's quantum columns beat every k-mer encoding of the same bases | +0.0114 | [C] | `curated_ladder.csv` | `python -m sgrna.curated_ladder` |
| Of the 63 hand-named columns, the computable 46 carry the gain | +0.0065 vs −0.0007 | [C] | `handnamed_split.csv` | `python -m sgrna.handnamed_split` |
| `d_mechanics`'s gain is flanking DNA, not re-description | +0.097 of +0.099 | [C] | `no_information.csv`, `d_mechanics_reconstructibility.csv` | `python -m sgrna.no_information`; `python -m sgrna.rederive --only d_mechanics_reconstructibility` |
| Model class barely matters, with or without the one real effect | 0.1076 ρ over 8 classes | [C] | `curated_models.csv` | `python -m sgrna.curated_models` |
| Summing the family table would overstate it 3.5× | +0.2030 against +0.7046 | [C] | `curated_joint.csv` | `python -m sgrna.curated_joint` |
| The 300-column budget and the selector are both free | ±0.0005 | [C] | `selector_control.csv`, `selector_comparison.csv` | `python -m sgrna.rederive --only selector_control` |

## Mechanism — what it is, and what it is not

| claim | number | frame | evidence | command |
|---|---|---|---|---|
| **The effect is specific to cleavage, not Cas9 binding** | −0.009 vs +0.157 | Cui | `cui_screen.csv` | `python -m sgrna.cui_screen` |
| Orientation to the nearest TSS matters | Δ 0.051 [0.020, 0.082] | [P] | `asymmetry_strand_audit.csv` | `python -m sgrna.asymmetry`; the TSS control is Jacky's |
| Not helical phasing: GC skew and periodicity are null | −0.001 alone, +0.0004 on A | [P] | `ablation_jper.csv`, `ablation_jper_joint.csv` | `python -m sgrna.run_ablation --families j_periodicity --tag jper` |
| Not flank shape: shape is a reparameterisation of composition | +0.0007 on top of A | [P] | `ablation_shape_joint.csv`, `cross_redundancy_*.csv` | `python -m sgrna.run_ablation --families i_shape --tag shape_joint` |
| Not supercoiling: the no-antibody control predicts as well | control ≈ ChIP | [P] | `dissection_supercoiling.csv` | `python -m sgrna.rederive --only dissection_supercoiling` |
| Not essentiality (Guo confirmed); expression is not null | +0.00004 vs +0.0213 | [C] | `txn_interactions.csv` | `python -m sgrna.txn_interactions` |
| The family's own interaction rationale fails | +0.0004 over parents | [C] | `txn_interactions.csv` | `python -m sgrna.txn_interactions` |
| The label's ceiling, from crisprHAL's own replicate counts | 0.946 [0.944, 0.948] | [C] | `curated_ceiling.json` | `python -m sgrna.curated_ceiling` |

## Against the state of the art

| claim | number | frame | evidence | command |
|---|---|---|---|---|
| Parity on crisprHAL 2's own eSpCas9 split, their weights | Δρ +0.0021 [−0.0031, +0.0072] | their split | `final_eval_espcas9.csv`, `final_eval_metric_bootstrap.csv` | `python -m sgrna.final_eval --arm compare`; `python -m sgrna.rederive --only final_eval_metric_bootstrap` |
| Parity on their SpCas9 split | +0.005 | their split | `their_split_slicer.csv`, `data_ablation.csv` | `python -m sgrna.rederive --only their_split_slicer` |
| Ahead on our own held-out folds, paired on guides | +0.0115 / +0.0172 / +0.0043 | [C] | `heldout_compare.csv` | `python -m sgrna.heldout_compare` |
| A search of their own finds nothing better than their published config | 0.6987, 4/5 inner folds | [C] | `matched_search_folds.csv`, `crisprhal_matched_folds.csv` | `colab/crisprhal_tune_matched.py` (Colab; the CSVs were transcribed from its log †) |
| The conservative margin, against their best known configuration | +0.0161 | [C] | `crisprhal_tuned_folds.csv`, `tuning.csv` | `python -m sgrna.tune` |
| **Holding the architecture fixed, parity is bought with rows** | label +0.0063, which rows +0.0083, how many +0.0302 | their split | `data_ablation.csv` | `python -m sgrna.data_ablation` |
| …and the feature set still dominates every data effect | +0.1554 | their split | `data_ablation.csv` | `python -m sgrna.data_ablation` |
| Under equal input DNA the ranking depends on the window | theirs +0.023 at 378 nt, ours +0.065 at ±1 kb | [C] | `crisprhal_wide_window.csv` (`seqnet_matched.csv` †) | `colab/crisprhal_wide_window.py`, `python -m sgrna.seqnet` |
| On an assay neither model trained on, theirs transfers better | −0.0258 [−0.0289, −0.0227] | Cui | `cui_final_exam.csv` | `python -m sgrna.cui_final_exam` |
| Data efficiency: the lead is a small-data lead | +0.0133 at 26.9k, +0.0021 at 47.6k | both | `sample_complexity.csv`, `learning_curve.csv` | `python -m sgrna.sample_complexity` |
| Cost: ~150× in CPU-seconds, 92× wall-clock, comparable memory | 65.8 s vs ~10,000 s per fold | — | `cost_cpu_seconds.csv`, `cost_flops.csv`, `cost_panel.csv`, `memory_comparison.csv` | `python -m sgrna.cost` |
| Model size: 6,348 fitted leaf values against 790,533 parameters | 125× | — | `complexity.json` | `python -m sgrna.complexity` |
| **SLICER against Noshay's own model, on Noshay's data** | iRF 0.4785 against SLICER 0.6060 | [P] | `noshay_replication.csv`, `ablation_main.csv`, `shift_null.csv` | `python -m sgrna.noshay_replicate` |

## Where it stops

| claim | number | frame | evidence | command |
|---|---|---|---|---|
| Transfer across nuclease, and across nuclease + organism | 88% / 91% retained | external | `transfer_model.csv`, `transfer.csv` | `python -m sgrna.transfer` |
| No transfer across kingdoms, in either direction | −0.039 [−0.054, −0.024]; −0.006 [−0.023, +0.011] | both | `kingdom_transfer.csv`, `kingdom_transfer_boot.csv` | `python -m sgrna.kingdom_transfer --bootstrap` |
| The human baseline reproduces at the teammate's settings | 0.5076 against 0.472 quoted | human | `human_baseline.csv`, `human_baseline_knobs.csv` | `python -m sgrna.human_baseline` |
| The human signal is far less additive than *E. coli*'s | ridge/boosting 63.7% vs 95.6% | human | `human_additivity.csv` | `python -m sgrna.human_baseline --additivity` |
| Adaptive overfitting is bounded, not measured | +0.0005 ± 0.016, so ≲0.02 | [C] | `sealed_holdout.csv` | `python -m sgrna.sealed_holdout --seal-seeds 3` |
| The 3.3× inflation is the row population, not the label | label +0.004; population ≈ +0.070 | both | `why_inflation.csv`, `ablation_label_published.csv`, `ablation_label_curated.csv` | `python -m sgrna.why_inflation` |
| **and the mechanism is a steeper slope outside Guo's 1,222-gene panel** | β +5.4 inside, +17.0 outside | [C] | `curation_variance.csv` | `python -m sgrna.curation_variance` |
| Gene identity explains more of the curated label than the published one | 0.844 → 0.713 within-gene share | both | `cohort_groups.csv` | `python -m sgrna.cohort_groups` |

---

## † Files whose producer was lost, and what re-deriving them found

Thirteen files cited above were written by code no longer in `src/` or `colab/`.
`python -m sgrna.rederive` now reproduces ten of them from committed code and
**compares rather than overwrites** — re-derivation is a check, and the stored
file stays authoritative. Verdicts in `results/rederivation_check.csv`, outputs
in `results/rederived/`:

| file | verdict |
|---|---|
| `grouped_retention.csv` | **exact** |
| `selector_control.csv` | **exact** |
| `final_eval_metric_bootstrap.csv` | **exact** point estimates; CIs differ by bootstrap noise only |
| `portability_364.csv` | agrees, 0.2% |
| `kmer_order_noselect.csv` | agrees, 0.2% |
| `dissection_supercoiling.csv` | agrees on the claim (control ≥ signal), 1.2%; the stored arm names and `bin100k` grouping could not be reconstructed |
| `their_split_slicer.csv` | agrees, ρ 0.6997 against 0.6985; all metrics within 0.004 except a top-10 statistic |
| `d_mechanics_reconstructibility.csv` | agrees in distribution, median 0.832 against 0.845; the file is keyless so rows cannot be aligned |
| `additivity_reduced.csv` | **two of four arms exact, two wrong** — see below |
| `abund_local_vs_mappability.csv` | **not reproducible**: the definition of its `uniq` columns is lost. Exact-20-mer uniqueness is near-constant in *E. coli* and gives correlations an order of magnitude weaker, so the stored quantity was something else |
| `seqnet_matched.csv` | re-derivable but not re-run here: `python -m sgrna.seqnet --flank {100,250,500,1000} --out results/seqnet_matched.csv` |
| `crisprhal_matched_folds.csv`, `matched_search_folds.csv` | GPU artefacts transcribed from a Colab log. `colab/crisprhal_tune_matched.py` now writes them, but rerunning needs a T4 |

**The exercise paid for itself.** `additivity_reduced.csv`'s two reduced arms
re-derive to within 0.0007 (and its XGBoost column bit-identically), but its two
full-matrix arms do not: the stored "ridge on the raw 6,232 collapses to 0.423"
holds only at `alpha=1`, and a penalty search gives **0.4945**. That claim is
withdrawn in the report. A number nobody can regenerate is a number nobody can
argue with — which is exactly what happened to `human_transfer.csv`.

## Reproducing any row

Every command above assumes `pip install -e .` and the data layout in
[`../README.md`](../README.md). Seeds are fixed (41–45) and every table here was
produced by the committed code in [`../src/sgrna`](../src/sgrna); the two
GPU-dependent rows run on free Colab from [`../colab`](../colab). The harness
check is `ablation_main.csv`'s baseline row, which must reproduce R² 0.2937 /
ρ 0.5278 to four decimal places before any other number is trusted.
