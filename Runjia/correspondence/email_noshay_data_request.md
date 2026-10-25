# Draft email — data request to Noshay et al. (NAR 2023)

**To:** jacobsonda@ornl.gov (Daniel A. Jacobson, corresponding author)
**Cc:** d.kainer@uq.edu.au (David Kainer, co-corresponding)
**Subject:** Two data questions about Noshay et al. 2023 (NAR 51:10147) — QCT sgRNA matrix

---

Dear Professor Jacobson and Dr Kainer,

I'm a student working on sgRNA cutting-efficiency prediction in *E. coli*, and
your 2023 NAR paper is the starting point for our project — we build on the
quantum chemical tensor feature matrix from Supplementary Table 2. Two things
we can't resolve from the published files, and I'd be grateful for your
guidance on either.

**1. The feature matrix for the full 40,468 sgRNAs.** Table 1 reports the full
matrix model on 40,468 sgRNAs (32,374 train), but Supplementary Table 2 covers
13,880. We reimplemented the iRF model from your methods section — 1,000 trees
per forest, ten iterations, five-fold cross-validation on all 6,232 columns — and
on the 13,880 released rows it gives **Pearson 0.4881 ± 0.0119 and R²
0.2345 ± 0.0102**, against your 0.5019 and 0.2491. The shortfall is about what
we would expect from training on a third of the rows, so we read it as agreement
rather than a discrepancy. If the 40,468-row matrix is shareable it would let us
check that directly, and would roughly triple the rows available to us. If it
isn't, a one-line confirmation of which guides the Table 1 models used would
still help us describe the comparison accurately.

**2. Replicate-level counts behind the cutting scores.** We've been trying to
estimate how much of the label is measurement noise, because that sets a ceiling
on what any model can achieve. The standard approach needs two independent
measurements of the same guides. The only pairs available to us are the
WT-SpCas9 and eSpCas9 screens, and those are not replicates — we can predict
the *disagreement* between them at ρ 0.55, so a large part of it is systematic
enzyme behaviour rather than noise, and the ceiling estimate inherits that
confound. Per-replicate or per-timepoint read counts for the Guo et al. screens
you worked from would settle it directly. If you still hold the intermediate
counts from your processing, or can point us to where the replicate-level data
lives, that would be very helpful. We appreciate this may simply not exist in a
shareable form.

**Something that may be useful to you in return.** To apply your feature set to
guides outside the published table, we reverse-engineered it: the `V####`
columns are position-specific k-mer indicators (5,853 of 5,887 identified, no
ambiguities), and the quantum tensor columns recover exactly as k-mer lookup
tables with no conflicts. Together that regenerates 6,169 of the 6,232 columns
from any 20-nt protospacer, verified at 100% agreement over 1.17 million
rebuilt values. The remaining ~63 hand-named columns we could not reproduce
exactly, which is the one place we'd welcome a correction. We're happy to share
the decoding code, and if it's useful we would be glad for it to be linked from
your repository — it makes the resource usable on any screen rather than only
on the guides already processed.

Thank you for making the matrix public in the first place; the project would not
exist without it.

With thanks,

Runjia Wang
RSI09 — sgRNA efficiency project
wangjinwei2012@gmail.com

---

## Notes before sending — not part of the email

- **Check the 40,468 figure** against Table 1 in the PDF before sending; it is
  quoted here from `lit_review/gkad736.pdf` (page with "Table 1. iRF model
  summary and metrics"), where the E. coli full-matrix row reads
  `40 468 [32 374 train]`, R² 0.24906667, Pearson 0.5019173.
- **The reproduction claim is now measured** (`results/noshay_replication.csv`,
  5 folds, 2026-10-04): Pearson 0.4881 ± 0.0119, R² 0.2345 ± 0.0102, and
  Spearman 0.4785 — the last of which their paper does not report, so it is not
  mentioned in the email but is worth offering if they reply.
- One thing deliberately **not** in the email: that swapping iRF for gradient
  boosting on their own matrix and folds is worth +0.049 Spearman. It is a fair
  result and it is in our report, but leading a data request with "your model is
  beatable" would be a poor opening. Raise it later if a conversation starts.
- The replicate ask is really about **Guo et al. 2016/2018**, not about this
  paper's own work, so it is framed as "the screens you worked from" and as a
  pointer request rather than a demand.
- The offer at the end is genuine and costs nothing: it is the single thing most
  likely to make a reply worth their time.
- If no reply in ~2 weeks, the fallback is re-running the counting pipeline from
  SRA (PRJNA450978 and others), which the report already notes as the expensive
  route.
