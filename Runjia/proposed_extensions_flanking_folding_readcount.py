# ============================================================
# PROPOSED EXTENSIONS: flanking-window QCT, real RNA folding, label filtering
#
# Implements three reviewer suggestions for the E. coli sgRNA cut-score model:
#   A. Extend the descriptor window past the 20-mer protospacer into the
#      ~11 nt downstream of the PAM ("target locus" instead of just
#      "protospacer").
#   B. Real ViennaRNA folding features (spacer MFE, spacer:scaffold duplex
#      stability) instead of sequence-composition proxies, with the
#      reported sigmoidal thresholds at -5 and -15 kcal/mol.
#   C. Filter the training label by control read count before fitting.
#
# STATUS / WHAT IS AND ISN'T WIRED UP:
#   ecoli_feature_matrix.csv (the only project dataset in this repo) has no
#   raw sequence column and no control-read-count column -- it is a
#   pre-aggregated supplementary table (see its own row-1 title:
#   "Supplemental Table 2: E.coli iRF feature matrix"). So:
#     - Section B (ViennaRNA folding) needs only a sequence string, so it
#       is fully implemented and runnable right now (see the __main__ demo
#       at the bottom).
#     - Section A splits into two parts. A1 (flank_composition_features:
#       AT fraction, GC skew, melting temp) is fully computable from raw
#       sequence alone and is now implemented and runnable, same as
#       Section B — see genomic_coordinate_lookup.py for how to get real
#       flank sequences to feed it. A2 (compute_flanking_qct_features,
#       the per-position quantum descriptors) still needs (1) the QCT
#       descriptor function itself -- referenced by the reviewer as "your
#       QCT pipeline", which lives with Joshua/Jacky's offline DFT-based
#       computation and is not reusable code in this repo -- and (2) full
#       target-locus sequences, now obtainable via
#       genomic_coordinate_lookup.py once you have a local MG1655
#       GenBank file. A2 is accepted as a parameter so it slots in once
#       available; calling it without one raises a clear error rather
#       than silently doing nothing.
#     - Section C (read-count filtering) needs a control-read-count
#       column that does not exist in the current CSV. It raises a clear
#       error naming the missing column rather than silently no-op'ing.
#
#   The Spearman deltas quoted in the request (+0.143 for flanking context,
#   +0.029-0.036 for read-count filtering) are as reported by the reviewer's
#   "curated re-analysis" -- they have NOT been reproduced against this
#   project's own data, since the raw sequence/read-count files needed to
#   do so aren't in this repo yet (per team discussion, they exist in
#   datasets cited by papers in the literature-review folder, not yet
#   pulled in locally). Treat them as a claim to validate, not a
#   confirmed result, until Sections A and C are run against real inputs.
# ============================================================

# %% CELL 0 — Install missing packages
import sys
import subprocess


def install_if_missing(import_name, pip_name=None):
    pip_name = pip_name or import_name
    try:
        __import__(import_name)
    except ImportError:
        print(f"Installing {pip_name}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pip_name])


for _mod, _pip in [("numpy", None), ("pandas", None), ("RNA", "ViennaRNA")]:
    install_if_missing(_mod, _pip)

# %% CELL 1 — Imports
import numpy as np
import pandas as pd
import RNA  # ViennaRNA Python bindings


# ============================================================
# SECTION A — Extend the descriptor window past the 20-mer
# ============================================================
# A1. Composition features (AT fraction, GC skew, melting temp) — these are
# genuinely computable right now, no DFT/QCT pipeline needed. The literature
# review lists them separately from the quantum descriptors for exactly
# this reason ("AT-richness of the wider flank carries real signal too").
def flank_composition_features(flank_seq: str) -> dict:
    """AT fraction, GC skew, and an approximate melting temperature for a
    flanking-window sequence (upstream or downstream of the protospacer).
    Unlike compute_flanking_qct_features below, this needs only the raw
    sequence — no external QCT descriptor function required.

    GC skew follows the standard (G-C)/(G+C) definition; returns 0.0 for
    a window with no G or C rather than raising, since a skew of 0 is the
    sensible convention when the denominator is zero (not a missing value).

    Melting temperature uses the basic Wallace rule (4*(G+C) + 2*(A+T)),
    which is a coarse approximation valid mainly for short windows (<~50
    nt); for a large ±200 nt window as used in genomic_coordinate_lookup.py,
    prefer a proper nearest-neighbor Tm (e.g. Biopython's
    Bio.SeqUtils.MeltingTemp.Tm_NN) over this if precision matters more
    than dependency-simplicity.
    """
    seq = flank_seq.upper()
    n = len(seq)
    if n == 0:
        raise ValueError("flank_composition_features got an empty sequence.")
    a, t, g, c = (seq.count(b) for b in "ATGC")
    at_fraction = (a + t) / n
    gc_skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
    tm_wallace = 4 * (g + c) + 2 * (a + t)
    return {
        "at_fraction": at_fraction,
        "gc_skew": gc_skew,
        "tm_wallace_c": float(tm_wallace),
    }


def add_flank_composition_features(df: pd.DataFrame, flank_seq_col: str) -> pd.DataFrame:
    """Vectorized wrapper: adds the three flank_composition_features()
    columns, prefixed eng.flank. to match ENGINEERED_PREFIXES."""
    if flank_seq_col not in df.columns:
        raise KeyError(
            f"add_flank_composition_features needs a {flank_seq_col!r} column "
            f"with real flanking sequences; not found (columns: "
            f"{list(df.columns)[:10]}...). Get this from "
            "genomic_coordinate_lookup.py's locus window extraction first."
        )
    feats = df[flank_seq_col].apply(flank_composition_features)
    feats_df = pd.DataFrame(list(feats), index=df.index).add_prefix("eng.flank.")
    return df.join(feats_df)


# A2. Per-position quantum descriptors over the flanking window — this is
# the part that genuinely needs the DFT-based QCT pipeline (Joshua/Jacky's
# offline descriptor computation), not something re-derivable from the
# sequence alone. Kept as a generic looping function with the descriptor
# computation injected, per the original design below.
def compute_flanking_qct_features(locus_sequences, protospacer_len, qct_descriptor_fn,
                                   downstream_nt=11):
    """Run a per-position QCT descriptor function over the ~11 nt
    downstream of the PAM, in addition to the protospacer itself.

    This is deliberately just the *looping* logic -- per the reviewer,
    running an existing per-position descriptor pipeline over more
    positions is mechanical. The actual descriptor computation is
    injected via `qct_descriptor_fn` rather than reimplemented here,
    since that logic already exists (Joshua/Jacky's QCT feature code)
    and duplicating it would risk diverging from the version already
    validated on the protospacer.

    Parameters
    ----------
    locus_sequences : list[str]
        Full target-locus sequences: protospacer + >= downstream_nt of
        flanking sequence past the PAM, one per guide.
    protospacer_len : int
        Length of the protospacer (20 for the 20-mer).
    qct_descriptor_fn : Callable[[str, int], dict]
        Your existing per-position QCT descriptor function. Must accept
        (locus_sequence, position_index) and return a dict of
        {descriptor_name: value} for that position. This is the
        "QCT pipeline" referenced in the request -- plug in the real
        implementation here.
    downstream_nt : int
        How many nt past the protospacer to score. Default 11, matching
        the reviewer's finding that most of the flanking-context signal
        is captured within ~11 nt downstream of the PAM.

    Returns
    -------
    pd.DataFrame, one row per input sequence, columns named
    "flank{offset}.{descriptor_name}" for offset in
    [protospacer_len, protospacer_len + downstream_nt).
    """
    rows = []
    for seq in locus_sequences:
        min_len = protospacer_len + downstream_nt
        if len(seq) < min_len:
            raise ValueError(
                f"Locus sequence shorter ({len(seq)} nt) than protospacer_len "
                f"+ downstream_nt ({min_len} nt): {seq!r}. Flanking features "
                "need the full target-locus sequence, not just the 20-mer."
            )
        row = {}
        for offset in range(protospacer_len, protospacer_len + downstream_nt):
            for name, value in qct_descriptor_fn(seq, offset).items():
                row[f"flank{offset}.{name}"] = value
        rows.append(row)
    return pd.DataFrame(rows)


# ============================================================
# SECTION B — Real ViennaRNA folding features
# ============================================================
def spacer_mfe(spacer_seq):
    """Minimum free energy (kcal/mol) of the spacer folding on itself,
    via ViennaRNA's RNAfold. More negative = more stable secondary
    structure = spacer more likely to be sequestered instead of pairing
    with target DNA."""
    _structure, mfe = RNA.fold(spacer_seq.replace("T", "U"))
    return float(mfe)


def duplex_mfe(spacer_seq, scaffold_seq):
    """Hybridization free energy (kcal/mol) between the spacer and the
    sgRNA scaffold, via ViennaRNA's RNAduplex. More negative = stronger
    (mis-)pairing between spacer and scaffold, competing with correct
    sgRNA folding."""
    duplex = RNA.duplexfold(spacer_seq.replace("T", "U"), scaffold_seq.replace("T", "U"))
    return float(duplex.energy)


def folding_features(spacer_seq, scaffold_seq,
                      spacer_mfe_threshold=-5.0, duplex_mfe_threshold=-15.0):
    """Real ViennaRNA folding features for one guide, replacing sequence-
    composition proxies. Includes both the continuous MFE values and
    binary threshold features at the reviewer-reported sigmoidal
    breakpoints (-5 kcal/mol for spacer self-folding, -15 kcal/mol for
    spacer:scaffold duplex stability) -- a sharp threshold is exactly the
    case where a continuous feature alone under-informs a tree/linear
    model, so both forms are kept."""
    s_mfe = spacer_mfe(spacer_seq)
    d_mfe = duplex_mfe(spacer_seq, scaffold_seq)
    return {
        "spacer_mfe": s_mfe,
        "spacer_mfe_below_neg5": float(s_mfe < spacer_mfe_threshold),
        "duplex_mfe": d_mfe,
        "duplex_mfe_below_neg15": float(d_mfe < duplex_mfe_threshold),
    }


def add_folding_features(df, spacer_col, scaffold_col):
    """Vectorized wrapper: adds the four folding_features() columns to
    every row of df using its spacer/scaffold sequence columns."""
    for col in (spacer_col, scaffold_col):
        if col not in df.columns:
            raise KeyError(
                f"add_folding_features needs a {col!r} column with real sequences; "
                f"not found in the given DataFrame (columns: {list(df.columns)[:10]}...)."
            )
    feats = df.apply(lambda r: folding_features(r[spacer_col], r[scaffold_col]), axis=1)
    return df.join(pd.DataFrame(list(feats), index=df.index))


# ============================================================
# SECTION C — Filter the label by control read count
# ============================================================
def filter_by_control_reads(df, control_read_col, min_reads):
    """Drop rows below a minimum control-read-count threshold before
    fitting. This changes y (label quality), not X (features) -- per the
    reviewer, this can matter more for the CV ceiling than any single
    feature family, since a noisy label caps every model equally
    regardless of how good the features are.

    Raises rather than silently no-op'ing if the column doesn't exist,
    since a silently-skipped filter would make CV numbers look like they
    used clean labels when they didn't."""
    if control_read_col not in df.columns:
        raise KeyError(
            f"filter_by_control_reads needs a {control_read_col!r} column; not "
            f"present in ecoli_feature_matrix.csv or the given DataFrame. This "
            "filter needs the original screen's per-guide control read counts, "
            "not something derivable from the pre-aggregated feature matrix."
        )
    n_before = len(df)
    filtered = df[df[control_read_col] >= min_reads].copy()
    n_after = len(filtered)
    pct_dropped = 100.0 * (n_before - n_after) / n_before if n_before else 0.0
    print(f"Control-read filter (>= {min_reads} reads): kept {n_after}/{n_before} "
          f"rows ({pct_dropped:.1f}% dropped).")
    return filtered


# ============================================================
# DEMO — Section B is fully runnable today; A and C need real inputs
# ============================================================
if __name__ == "__main__":
    print("Demo: Section B (ViennaRNA folding) on example sequences.")
    example_spacer = "GACGCATAAAGATGAGACGC"  # 20-mer, placeholder sequence
    example_scaffold = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC"
    print(folding_features(example_spacer, example_scaffold))

    print("\nDemo: Section A1 (flank composition) on a placeholder flank window.")
    example_flank = "AATATTGCGCGATATTACGCGCGATTATATAGCGCGCGATATTAGCGCGATATAGCGCG"
    print(flank_composition_features(example_flank))

    print(
        "\nSection A2 (per-position QCT descriptors) and Section C (read-count "
        "filtering) are implemented but intentionally not demoed here: they "
        "need a QCT descriptor function (A2) and a real control-read-count "
        "column (C), neither of which exist in ecoli_feature_matrix.csv yet. "
        "Calling them without those inputs raises a clear error rather than a "
        "silent no-op -- see each function's docstring for exactly what to "
        "pass in once the raw screen data is pulled in from the literature-"
        "cited sources. For real flank sequences to feed A1 and A2 both, see "
        "genomic_coordinate_lookup.py."
    )
