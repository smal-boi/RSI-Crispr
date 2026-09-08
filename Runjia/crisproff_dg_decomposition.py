# ============================================================
# CRISPRoff-STYLE ENERGY DECOMPOSITION (literature-review family B)
#
# Implements Section 5, Step 4: dG_H (hybridization), dG_U (guide
# unfolding penalty), dG_O (DNA opening penalty), net dG_B, plus a
# distance-from-optimum term for dG_H's reported sweet spot.
#
# STATUS / IMPORTANT CAVEAT:
#   This is a from-scratch reimplementation of the CRISPRoff/CRISPRspec
#   energy model (Alkan et al. 2022) built directly on ViennaRNA
#   primitives -- the actual CRISPRoff 1.1.1 tool/repo was not fetched
#   (this sandbox has no network access to pull it), so these numbers
#   have NOT been cross-checked against the published tool's output.
#   The decomposition follows the same physical logic the paper
#   describes (net binding = hybridization gain minus the two unfolding/
#   opening penalties), and reuses spacer_mfe() from
#   proposed_extensions_flanking_folding_readcount.py for dG_U so the two
#   scripts stay consistent with each other. Before reporting these as
#   "CRISPRoff features" in the paper, either (a) install the real
#   CRISPRoff package (github.com/RTH-tools/crisproff) on a machine with
#   network access and compare its output against this script's on a
#   handful of guides, or (b) describe them explicitly as "a CRISPRoff-
#   style decomposition" rather than CRISPRoff's own output, since that's
#   what they actually are until cross-checked.
#
#   dG_O (DNA duplex opening penalty) in particular is the roughest
#   approximation here: CRISPRoff computes it from a specific nearest-
#   neighbor DNA:DNA stability table; this script uses ViennaRNA's DNA
#   parameter set as a stand-in. Treat dG_O as directionally right
#   (more GC-rich target = harder to open = larger positive penalty) but
#   not numerically identical to the published tool.
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

# Reported sweet-spot window for dG_H (Alkan et al. 2022), highly efficient
# guides cluster here. Used for the |dG_H - optimum| distance feature.
DG_H_SWEET_SPOT_LOW = -64.5
DG_H_SWEET_SPOT_HIGH = -47.1
DG_H_SWEET_SPOT_MID = (DG_H_SWEET_SPOT_LOW + DG_H_SWEET_SPOT_HIGH) / 2.0


# %% CELL 2 — Component energies
def spacer_unfolding_penalty(spacer_seq: str) -> float:
    """dG_U: the energetic cost of unfolding the spacer's own secondary
    structure before it can hybridize with the target. A spacer with a
    stable self-structure (very negative MFE) pays more to unfold, so
    dG_U is reported here as the POSITIVE cost (-MFE), not the MFE
    itself -- higher dG_U = more of a penalty."""
    _structure, mfe = RNA.fold(spacer_seq.replace("T", "U"))
    return float(-mfe)


def dna_opening_penalty(target_dsdna_seq: str) -> float:
    """dG_O: the energetic cost of locally melting/opening the target
    dsDNA duplex so the spacer can invade and form an R-loop. Approximated
    here via ViennaRNA's DNA nearest-neighbor parameters (RNA.params_load_DNA_Mathews2004)
    applied to the target strand paired with its complement -- i.e. the
    duplex stability of the dsDNA region being opened, reported as a
    positive cost (more stable duplex = larger opening penalty)."""
    RNA.params_load_DNA_Mathews2004()
    fc = RNA.fold_compound(target_dsdna_seq)
    complement = target_dsdna_seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    duplex = RNA.duplexfold(target_dsdna_seq, complement)
    RNA.params_load_RNA_Turner2004()  # restore default RNA params for other calls
    return float(-duplex.energy)


def hybridization_energy(spacer_seq: str, target_strand_seq: str) -> float:
    """dG_H: RNA:DNA hybridization free energy between the spacer and the
    protospacer target strand (the strand the spacer base-pairs with,
    i.e. the strand complementary to the protospacer-as-written). Treated
    as an RNA:RNA duplex via ViennaRNA (T->U substitution) as a stand-in
    for RNA:DNA hybrid parameters, which ViennaRNA does not natively
    support -- the same simplification the existing spacer:scaffold
    duplex_mfe() in proposed_extensions_flanking_folding_readcount.py
    already makes."""
    duplex = RNA.duplexfold(
        spacer_seq.replace("T", "U"), target_strand_seq.replace("T", "U")
    )
    return float(duplex.energy)


# %% CELL 3 — Full decomposition
def crisproff_style_features(spacer_seq: str, target_strand_seq: str) -> dict:
    """CRISPRoff-style energy decomposition for one guide. See module
    docstring for the caveat on numerical fidelity vs. the published tool.

    Returns dG_H, dG_U, dG_O, dG_B (= dG_H - dG_U - dG_O, net binding —
    hybridization gain minus the two penalties), and the sweet-spot
    distance feature so the model doesn't have to discover the
    non-monotonic window itself.
    """
    dg_h = hybridization_energy(spacer_seq, target_strand_seq)
    dg_u = spacer_unfolding_penalty(spacer_seq)
    dg_o = dna_opening_penalty(target_strand_seq)
    dg_b = dg_h - dg_u - dg_o
    dist_from_optimum = (
        0.0 if DG_H_SWEET_SPOT_LOW <= dg_h <= DG_H_SWEET_SPOT_HIGH
        else min(abs(dg_h - DG_H_SWEET_SPOT_LOW), abs(dg_h - DG_H_SWEET_SPOT_HIGH))
    )
    return {
        "dg_h": dg_h,
        "dg_u": dg_u,
        "dg_o": dg_o,
        "dg_b": dg_b,
        "dg_h_dist_from_sweet_spot": dist_from_optimum,
        "dg_h_in_sweet_spot": float(DG_H_SWEET_SPOT_LOW <= dg_h <= DG_H_SWEET_SPOT_HIGH),
    }


def add_crisproff_features(df: pd.DataFrame, spacer_col: str, target_strand_col: str) -> pd.DataFrame:
    """Vectorized wrapper: adds the crisproff_style_features() columns,
    prefixed eng.dg. to match ENGINEERED_PREFIXES in the baseline script."""
    for col in (spacer_col, target_strand_col):
        if col not in df.columns:
            raise KeyError(
                f"add_crisproff_features needs a {col!r} column with real "
                f"sequences; not found (columns: {list(df.columns)[:10]}...)."
            )
    feats = df.apply(
        lambda r: crisproff_style_features(r[spacer_col], r[target_strand_col]), axis=1
    )
    feats_df = pd.DataFrame(list(feats), index=df.index)
    feats_df = feats_df.add_prefix("eng.dg.")
    return df.join(feats_df)


# ============================================================
# DEMO
# ============================================================
if __name__ == "__main__":
    print("Demo: CRISPRoff-style decomposition on example sequences.")
    example_spacer = "GACGCATAAAGATGAGACGC"  # 20-mer, placeholder sequence
    example_target_strand = "GACGCATAAAGATGAGACGC"  # placeholder: normally the
    # genomic target strand from genomic_coordinate_lookup.py's locus window,
    # not literally identical to the spacer -- swap in the real target
    # strand sequence once Step 1's coordinate lookup is wired up.
    result = crisproff_style_features(example_spacer, example_target_strand)
    for k, v in result.items():
        print(f"  {k}: {v}")
    print(
        "\nThis demo uses a placeholder target-strand sequence equal to the "
        "spacer, purely to show the function runs end-to-end -- swap in "
        "genomic_coordinate_lookup.py's real locus window sequence (reverse-"
        "complemented as needed) before trusting these numbers for anything."
    )
