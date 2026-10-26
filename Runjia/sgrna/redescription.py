"""How much is re-description worth, and does it have to mean anything?

Finding 2 says a deterministic function of the 20-mer adds no *information* but
does add *performance* -- about +0.015 to a depth-limited tree, agreed across
four families built from physics: CRISPRoff energies, ViennaRNA folding,
nearest-neighbour thermodynamics, quantum-chemical tables.

Two things follow that the report has not tested.

**First: if the mechanism is "the tree cannot cheaply compute this", the function
need not be physical.** Any deterministic function of the sequence that packages
a higher-order effect should do, including a meaningless one. But performance
cannot rise without bound by bolting on arbitrary functions, so the effect has a
size and a saturation point. `random_tables()` measures both, using functions
built to be structurally identical to the physical ones and semantically empty:
a random weight per (position, base), summed over the guide. That is exactly the
shape of a nearest-neighbour energy sum with the physics replaced by noise.

**Second: the family is the wrong unit.** A family is a folder, not a kind -- 8
of `d_mechanics`'s 30 columns read flanking DNA while 22 are the guide restated,
and the family average hides it. `column_taxonomy()` replaces the folder with a
measurement: regress every one of the 668 family columns on the guide's own
one-hot encoding out of fold, and split by how well the guide predicts it.

    python -m sgrna.redescription --arm random
    python -m sgrna.redescription --arm taxonomy
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

from . import config
from .curated_ablation import FAMILIES, _base, _block

RANDOM_OUT = config.RESULTS / "redescription_random.csv"
TAXONOMY_OUT = config.RESULTS / "redescription_taxonomy.csv"
COLUMNS_OUT = config.RESULTS / "redescription_columns.csv"

BASES = np.array(list("ACGT"))


def _onehot(protos: np.ndarray, order: int = 2) -> np.ndarray:
    """The guide itself, as indicators: 20x4 mono, plus 19x16 di if order >= 2."""
    seqs = np.array([list(str(p)[:20].upper().ljust(20, "A")) for p in protos])
    mono = np.concatenate(
        [(seqs[:, p:p + 1] == BASES[None, :]).astype(np.float32)
         for p in range(20)], axis=1)
    if order < 2:
        return mono
    pairs = np.array([a + b for a in "ACGT" for b in "ACGT"])
    di = np.concatenate(
        [(np.char.add(seqs[:, p], seqs[:, p + 1])[:, None] == pairs[None, :]
          ).astype(np.float32) for p in range(19)], axis=1)
    return np.hstack([mono, di])


def _oof_r2(Y: np.ndarray, A: np.ndarray, n_splits: int = 5) -> np.ndarray:
    """Out-of-fold R^2 of every column of Y predicted by least squares on A.

    One multi-output solve per fold rather than one per column: the design
    matrix is shared, so 668 columns cost no more than one.
    """
    A = np.hstack([A, np.ones((len(A), 1), dtype=np.float32)])
    pred = np.full(Y.shape, np.nan, dtype=np.float64)
    for tr, va in KFold(n_splits, shuffle=True, random_state=41).split(A):
        coef, *_ = np.linalg.lstsq(A[tr], np.nan_to_num(Y[tr]), rcond=None)
        pred[va] = A[va] @ coef
    ss_res = np.nansum((Y - pred) ** 2, axis=0)
    ss_tot = np.nansum((Y - np.nanmean(Y, axis=0)) ** 2, axis=0)
    # A column with no variance has no R^2 to report. Returning NaN here and
    # letting it fall through a chain of comparisons would silently file it
    # under whichever label the last `else` happens to be -- which is how 32
    # constant columns were first labelled "partly determined". They get their
    # own class instead.
    with np.errstate(invalid="ignore", divide="ignore"):
        return 1.0 - ss_res / np.where(ss_tot > 0, ss_tot, np.nan)


def _cv(X, y, n_features=None, seeds=(41, 42, 43), n_splits=5):
    import xgboost as xgb
    from scipy.stats import spearmanr

    from .run_ablation import _impute, champion_params, make_selector

    nf = config.N_FEATURES if n_features is None else n_features
    rhos = []
    for seed in seeds:
        for tr, va in KFold(n_splits, shuffle=True, random_state=seed).split(X):
            A, B = _impute(X[tr], X[va])
            if nf < A.shape[1]:
                sel = make_selector(seed)
                sel.fit(A, y[tr])
                top = np.argsort(sel.feature_importances_)[::-1][:nf]
                A, B = A[:, top], B[:, top]
            m = xgb.XGBRegressor(**champion_params(seed))
            m.fit(A, y[tr])
            rhos.append(float(spearmanr(y[va], m.predict(B)).statistic))
    return float(np.mean(rhos)), float(np.std(rhos, ddof=1))


# --------------------------------------------------------------------------

def random_tables(ks=(4, 16, 64, 256), table_seeds=(0, 1, 2),
                  model_seeds=(41, 42), verbose=True):
    """Meaningless functions of the guide, in the shape of physical ones.

    Three random draws per table count, two model seeds each: the spread across
    *draws* is what matters here -- one set of random tables is itself a sample,
    and a saturation curve read off a single draw would be unreadable.
    """
    X, y, names, protos = _base()
    seqs = np.array([list(str(p)[:20].upper().ljust(20, "A")) for p in protos])
    code = np.searchsorted(BASES, seqs)                  # (n, 20) in 0..3

    base_rho, base_sd = _cv(X, y, seeds=model_seeds)
    rows = [dict(arm="base", n_random=0, table_seed=-1, spearman=base_rho,
                 spearman_sd=base_sd, d_spearman=0.0)]
    if verbose:
        print(f"  base: rho {base_rho:.4f}", flush=True)

    for k in ks:
        for ts in table_seeds:
            rng = np.random.default_rng(1000 * k + ts)
            w = rng.normal(size=(k, 20, 4)).astype(np.float32)
            feats = np.empty((len(code), k), dtype=np.float32)
            for j in range(k):
                feats[:, j] = w[j][np.arange(20)[None, :], code].sum(axis=1)
            rho, sd = _cv(np.hstack([X, feats]), y, seeds=model_seeds)
            rows.append(dict(arm=f"base + {k} random tables", n_random=k,
                             table_seed=ts, spearman=rho, spearman_sd=sd,
                             d_spearman=rho - base_rho))
            if verbose:
                print(f"  {k:4d} random tables (seed {ts}): rho {rho:.4f}  "
                      f"d {rho - base_rho:+.4f}", flush=True)
            pd.DataFrame(rows).to_csv(RANDOM_OUT, index=False)

    out = pd.DataFrame(rows)
    if verbose:
        m = out[out.n_random > 0].groupby("n_random")["d_spearman"].mean()
        print("\n  mean gain by table count:")
        for k, v in m.items():
            print(f"    {k:5d}: {v:+.4f}")
        print(f"  for comparison, the physical families buy +0.011 to +0.019")
        print(f"written -> {RANDOM_OUT}")
    return out


def column_taxonomy(hi=0.99, lo=0.50, verbose=True):
    """Split all 668 family columns by how well the guide alone predicts them."""
    X, y, names, protos = _base()
    A = _onehot(protos, order=2)
    if verbose:
        print(f"  guide one-hot basis: {A.shape[1]} columns "
              f"(20x4 mono + 19x16 di)", flush=True)

    blocks, labels = [], []
    for fam in FAMILIES:
        try:
            v, fn = _block(fam, protos)
        except FileNotFoundError as exc:
            if verbose:
                print(f"  {fam}: skipped ({exc})")
            continue
        blocks.append(v)
        labels += [(fam, n) for n in fn]
    V = np.hstack(blocks)
    r2 = _oof_r2(V, A)

    cols = pd.DataFrame(labels, columns=["family", "column"])
    cols["oof_r2_from_guide"] = r2
    kind = np.where(np.isnan(r2), "constant (no variance)",
            np.where(r2 >= hi, "determined by the guide",
             np.where(r2 <= lo, "outside the guide", "partly determined")))
    cols["kind"] = kind
    cols.to_csv(COLUMNS_OUT, index=False)
    if verbose:
        print("\n  columns by kind:")
        print(cols["kind"].value_counts().to_string())
        print("\n  by family:")
        print(pd.crosstab(cols["family"], cols["kind"]).to_string())

    base_rho, base_sd = _cv(X, y)
    rows = [dict(arm="base", n_cols=0, spearman=base_rho, spearman_sd=base_sd,
                 d_spearman=0.0)]
    if verbose:
        print(f"\n  base: rho {base_rho:.4f}", flush=True)
    for kind in ("determined by the guide", "partly determined",
                 "outside the guide", "constant (no variance)"):
        idx = np.where(cols["kind"].to_numpy() == kind)[0]
        if len(idx) == 0:
            continue
        rho, sd = _cv(np.hstack([X, V[:, idx]]), y)
        rows.append(dict(arm=f"base + {kind}", n_cols=len(idx), spearman=rho,
                         spearman_sd=sd, d_spearman=rho - base_rho))
        if verbose:
            print(f"  + {kind:26s} ({len(idx):3d} cols): rho {rho:.4f}  "
                  f"d {rho - base_rho:+.4f}", flush=True)
        pd.DataFrame(rows).to_csv(TAXONOMY_OUT, index=False)

    rho, sd = _cv(np.hstack([X, V]), y)
    rows.append(dict(arm="base + all 668", n_cols=V.shape[1], spearman=rho,
                     spearman_sd=sd, d_spearman=rho - base_rho))
    if verbose:
        print(f"  + all {V.shape[1]} family columns: rho {rho:.4f}  "
              f"d {rho - base_rho:+.4f}")
    out = pd.DataFrame(rows)
    out.to_csv(TAXONOMY_OUT, index=False)
    if verbose:
        print(f"written -> {TAXONOMY_OUT} and {COLUMNS_OUT}")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--arm", choices=["random", "taxonomy", "families", "both"],
                    default="both")
    a = ap.parse_args(argv)
    if a.arm in ("taxonomy", "both"):
        column_taxonomy()
    if a.arm in ("random", "both"):
        random_tables()
    if a.arm == "families":
        family_split()
    return 0



# --------------------------------------------------------------------------
# Every family split the way d_mechanics was
# --------------------------------------------------------------------------

FAMILY_SPLIT_OUT = config.RESULTS / "redescription_family_split.csv"


def family_split(families=("a_flank", "b_energy", "c_folding", "d_mechanics",
                           "f_methylation", "h_offtarget"),
                 hi=0.99, lo=0.50, verbose=True):
    """Does the d_mechanics lesson generalise? Split each family and measure.

    `d_mechanics` turned out to be 22 columns restating the guide bundled with 8
    reading flanking DNA, and the family-level number hid it. The column probe
    says three other families are mixtures too -- `a_flank` carries 28 columns
    the guide largely explains, `b_energy` is 9 guide columns and 1 outside,
    `h_offtarget` splits 6 against 7 -- so each is run as two arms here rather
    than one. A family whose halves behave alike was a fair bundle; one whose
    halves diverge was not.
    """
    X, y, names, protos = _base()
    A = _onehot(protos, order=2)
    base_rho, base_sd = _cv(X, y)
    rows = [dict(family="(base)", part="base", n_cols=0, spearman=base_rho,
                 spearman_sd=base_sd, d_spearman=0.0)]
    if verbose:
        print(f"  base: rho {base_rho:.4f}", flush=True)

    for fam in families:
        try:
            v, fn = _block(fam, protos)
        except FileNotFoundError as exc:
            if verbose:
                print(f"  {fam}: skipped ({exc})")
            continue
        r2 = _oof_r2(v, A)
        guide = np.where(r2 > lo)[0]          # determined or partly determined
        outside = np.where(r2 <= lo)[0]       # NaN columns are in neither
        const = np.where(np.isnan(r2))[0]
        for part, idx in (("guide-explained", guide), ("outside", outside),
                          ("constant", const),
                          ("whole family", np.arange(v.shape[1]))):
            if len(idx) == 0:
                continue
            rho, sd = _cv(np.hstack([X, v[:, idx]]), y)
            rows.append(dict(family=fam, part=part, n_cols=len(idx),
                             spearman=rho, spearman_sd=sd,
                             d_spearman=rho - base_rho))
            if verbose:
                print(f"  {fam:15s} {part:16s} ({len(idx):3d} cols): "
                      f"d {rho - base_rho:+.4f}", flush=True)
            pd.DataFrame(rows).to_csv(FAMILY_SPLIT_OUT, index=False)

    out = pd.DataFrame(rows)
    if verbose:
        print(f"written -> {FAMILY_SPLIT_OUT}")
    return out

if __name__ == "__main__":
    raise SystemExit(main())
