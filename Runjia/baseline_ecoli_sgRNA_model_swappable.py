# ============================================================
# FIVE-SEED BASELINE MODEL: E. coli CRISPR sgRNA Cut-Score Prediction
# Project: Enhancing sgRNA efficiency in CRISPR-Cas9 genome editing
#
# This is a MODEL-SWAPPABLE variant of "USE THIS!! FINAL_baseline_
# ecoli_sgRNA_model_5seeds.py". Structure, leakage-safe CV, five-seed
# evaluation, and artifact saving are unchanged. The only addition is a
# MODEL_NAME switch (CELL 2) plus a build_champion() factory (CELL 3b)
# so the champion (prediction) model can be swapped without touching the
# CV loop. The feature SELECTOR stays XGBoost regardless of MODEL_NAME —
# see the note in CELL 3b for why.
#
# Report/trust the five-seed CV scores from CELL 7 / output CSVs.
# The final model trained on all data in CELL 8 is for production/inference
# only; it is not a validation score.
#
# Meant to run on a LOCAL Colab runtime (i.e. Colab connected to your own
# machine's kernel): no google.colab import, no Drive mount, plain local
# file paths, and CELL 0 pip-installs whatever's missing on your machine.
# ============================================================

# %% CELL 0 — Install missing packages if needed
import sys
import subprocess


def install_if_missing(import_name, pip_name=None):
    pip_name = pip_name or import_name
    try:
        __import__(import_name)
    except ImportError:
        print(f"Installing {pip_name}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pip_name])


for _mod, _pip in [
    ("numpy", None), ("pandas", None), ("sklearn", "scikit-learn"),
    ("scipy", None), ("xgboost", None), ("joblib", None),
]:
    install_if_missing(_mod, _pip)

# %% CELL 1 — Imports
import os
import json
import hashlib
import warnings
from collections import Counter

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import xgboost as xgb
import joblib
from scipy.stats import spearmanr, pearsonr
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_squared_error

# %% CELL 2 — Config
# Put this script in the same folder as ecoli_feature_matrix.csv, or edit DATA_PATH.
DATA_PATH = "ecoli_feature_matrix.csv"
TARGET_COLUMN = "cut.score"

# Frozen baseline settings.
N_FEATURES = 300
N_ESTIMATORS = 400

# --- Model switch -------------------------------------------------------
# Options: "xgboost" (default/original champion), "random_forest",
# "lightgbm", "ridge", "catboost". Add more cases in build_champion()
# (CELL 3b) as needed. Only installs the extra package for whichever
# model you pick.
MODEL_NAME = "xgboost"

_EXTRA_PACKAGE_BY_MODEL = {
    "lightgbm": "lightgbm",
    "catboost": "catboost",
}
if MODEL_NAME in _EXTRA_PACKAGE_BY_MODEL:
    install_if_missing(_EXTRA_PACKAGE_BY_MODEL[MODEL_NAME])
# --------------------------------------------------------------------

N_SPLITS = 5
SEEDS = [41, 42, 43, 44, 45]
FINAL_PRODUCTION_SEED = 42
SHUFFLE = True

OUTPUT_DIR = f"{MODEL_NAME}_ecoli"
os.makedirs(OUTPUT_DIR, exist_ok=True)

NON_FEATURE_COLUMNS = [
    "sgRNAID", "sequence", "Sequence", "sgRNA", "sgrna",
    "guide", "Guide", "guide_sequence", "GuideSequence",
    "target_sequence", "TargetSequence",
]

# Prefixes for engineered feature columns if they appear in a future CSV.
# This baseline script does not create engineered features; this is only metadata.
ENGINEERED_PREFIXES = ("regional_qct", "eng.")


def is_engineered_feature(name) -> bool:
    return str(name).strip().startswith(ENGINEERED_PREFIXES)


def engineered_family(name) -> str:
    s = str(name).strip()
    if s.startswith("regional_qct"):
        return "regional_qct"
    if s.startswith("eng."):
        return "eng_extra"
    return "original"


# %% CELL 3 — Helpers
def stable_spearman(y_true, y_pred) -> float:
    rho, _ = spearmanr(y_true, y_pred)
    return 0.0 if np.isnan(rho) else float(rho)


def stable_pearson(y_true, y_pred) -> float:
    rho, _ = pearsonr(y_true, y_pred)
    return 0.0 if np.isnan(rho) else float(rho)


def median_impute_train_val(X_train, X_val):
    """Median-impute NaNs (and +/-inf, treated as missing) using stats fit
    on the TRAIN fold only. Tree models (XGBoost/RandomForest/LightGBM)
    tolerate inf silently, but linear models like Ridge will throw
    'array must not contain infs or NaNs' if any slip through, so we
    sanitize inf -> NaN before imputing regardless of MODEL_NAME."""
    X_train = X_train.astype(np.float32, copy=True)
    X_val = X_val.astype(np.float32, copy=True)
    X_train[np.isinf(X_train)] = np.nan
    X_val[np.isinf(X_val)] = np.nan
    medians = np.nanmedian(X_train, axis=0)
    medians = np.where(np.isnan(medians), 0.0, medians)
    tr_rows, tr_cols = np.where(np.isnan(X_train))
    va_rows, va_cols = np.where(np.isnan(X_val))
    X_train[tr_rows, tr_cols] = medians[tr_cols]
    X_val[va_rows, va_cols] = medians[va_cols]
    return X_train, X_val, medians


def median_impute_full(X):
    """Median-impute full data (NaN and +/-inf) for the final production
    model only."""
    X = X.astype(np.float32, copy=True)
    X[np.isinf(X)] = np.nan
    medians = np.nanmedian(X, axis=0)
    medians = np.where(np.isnan(medians), 0.0, medians)
    rows, cols = np.where(np.isnan(X))
    X[rows, cols] = medians[cols]
    return X, medians


def clip_train_val(X_train, X_val, iqr_multiplier=3.0):
    """Winsorize extreme values using an IQR-based bound fit on the TRAIN
    fold only. Guards against outlier/sentinel values that are finite (so
    not caught by NaN/Inf sanitization above) but large enough to overflow
    Ridge's Gram matrix computation (X.T @ X) during Cholesky solve.

    Uses IQR rather than a fixed percentile cutoff: a sentinel/error value
    (e.g. a "missing" flag like -999 baked into a raw quantum-descriptor
    column) can appear in well over 0.5% of rows, in which case a fixed
    0.5/99.5 percentile clip leaves it untouched since it's still inside
    the kept range. Q1/Q3 come from the central 50% of each column, so
    they stay uncontaminated unless outliers make up more than 25% of one
    tail — far more robust to unknown contamination rates.

    Negligible effect on tree models since bounds this wide barely move
    their split points."""
    q1 = np.percentile(X_train, 25, axis=0)
    q3 = np.percentile(X_train, 75, axis=0)
    iqr = q3 - q1
    iqr_safe = np.where(iqr == 0, 1.0, iqr)  # avoid zero-width bounds on constant columns
    lower = q1 - iqr_multiplier * iqr_safe
    upper = q3 + iqr_multiplier * iqr_safe
    X_train_clipped = np.clip(X_train, lower, upper)
    X_val_clipped = np.clip(X_val, lower, upper)
    return X_train_clipped, X_val_clipped


def feature_name_hash(feature_names) -> str:
    payload = json.dumps(list(feature_names), separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def selector_params(seed: int) -> dict:
    """Feature selector always stays XGBoost regardless of MODEL_NAME: it
    needs a fast, reliable feature_importances_ signal, and changing the
    selector would make CV scores across different MODEL_NAME runs harder
    to compare apples-to-apples (they'd be selecting different features,
    not just predicting differently from the same features)."""
    return dict(
        n_estimators=300,
        learning_rate=0.05,
        max_depth=4,
        subsample=0.8,
        colsample_bytree=0.8,
        reg_alpha=0.0,
        reg_lambda=1.0,
        min_child_weight=1,
        random_state=seed,
        n_jobs=-1,
        tree_method="hist",
        importance_type="gain",
    )


def _to_float64(X):
    """Named (not lambda) so the Ridge pipeline stays picklable — CELL 9
    joblib.dump()s champion_model, and joblib/pickle can't serialize a
    lambda or other local/anonymous function."""
    return X.astype(np.float64)


# %% CELL 3b — Champion model factory (this is the swap point)
def build_champion(seed: int, n_estimators: int = N_ESTIMATORS):
    """Builds the prediction model used on the selected top-N_FEATURES
    columns. Swap MODEL_NAME in CELL 2 to change what this returns —
    nothing else in the script needs to change."""
    if MODEL_NAME == "xgboost":
        return xgb.XGBRegressor(
            n_estimators=n_estimators,
            learning_rate=0.05, max_depth=4,
            subsample=0.8, colsample_bytree=0.8,
            reg_alpha=0.0, reg_lambda=1.0, min_child_weight=1,
            random_state=seed, n_jobs=-1, tree_method="hist",
        )
    elif MODEL_NAME == "random_forest":
        from sklearn.ensemble import RandomForestRegressor
        return RandomForestRegressor(
            n_estimators=n_estimators, max_depth=10,
            random_state=seed, n_jobs=-1,
        )
    elif MODEL_NAME == "lightgbm":
        import lightgbm as lgb
        return lgb.LGBMRegressor(
            n_estimators=n_estimators, learning_rate=0.05, max_depth=4,
            subsample=0.8, colsample_bytree=0.8,
            random_state=seed, n_jobs=-1, verbosity=-1,
        )
    elif MODEL_NAME == "catboost":
        from catboost import CatBoostRegressor
        # depth=4 / learning_rate=0.05 mirror the XGBoost champion params
        # for a fair comparison. subsample+bootstrap_type="Bernoulli" is
        # CatBoost's equivalent of XGBoost's row subsampling; rsm is its
        # column-subsampling equivalent (colsample_bytree). CatBoost
        # natively handles NaN, but our data is already imputed upstream
        # so this is just a safety net.
        return CatBoostRegressor(
            iterations=n_estimators, learning_rate=0.05, depth=4,
            bootstrap_type="Bernoulli", subsample=0.8, rsm=0.8,
            random_seed=seed, thread_count=-1, verbose=False,
            allow_writing_files=False,
        )
    elif MODEL_NAME == "ridge":
        from sklearn.linear_model import Ridge
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler, FunctionTransformer
        # Linear baseline — useful to show whether the tree models are
        # actually adding value over a simple fit. n_estimators is unused.
        # Wrapped in StandardScaler: Ridge's regularization penalizes all
        # coefficients on the same scale, so unscaled features (some near
        # 0, others large) would get penalized very unevenly and can also
        # cause numerical instability in the solver even after inf/NaN
        # cleanup. Trees don't need this since splits are scale-invariant.
        #
        # Upcast to float64 first: our data stays float32 through
        # imputation/clipping (see CELL 3). If a selected column has a
        # tiny-but-nonzero variance, StandardScaler divides by that tiny
        # number and can produce values large enough that X.T @ X
        # overflows to inf in float32's ~3.4e38 range during Ridge's
        # Cholesky solve, even though every individual input value is
        # finite. float64's much larger dynamic range (~1.8e308) avoids
        # this; trees don't need this cast since they never form a Gram
        # matrix, so it's confined to this branch only.
        return make_pipeline(
            FunctionTransformer(_to_float64),
            StandardScaler(),
            Ridge(alpha=1.0, random_state=seed),
        )
    else:
        raise ValueError(f"Unknown MODEL_NAME: {MODEL_NAME!r}")


# %% CELL 4 — Load data
print(f"Loading: {DATA_PATH}")
if not os.path.exists(DATA_PATH):
    raise FileNotFoundError(
        f"Could not find {DATA_PATH!r}. Put this script in the same folder as "
        "ecoli_feature_matrix.csv, or edit DATA_PATH."
    )

df = pd.read_csv(DATA_PATH, header=1, low_memory=False)  # row 1 is a title row

required = {"sgRNAID", TARGET_COLUMN}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"CSV is missing columns {missing}.")

df = df.dropna(subset=["sgRNAID", TARGET_COLUMN])
y = pd.to_numeric(df[TARGET_COLUMN], errors="coerce").to_numpy(dtype=np.float64)

drop_cols = [TARGET_COLUMN] + [c for c in NON_FEATURE_COLUMNS if c in df.columns]
X_df = df.drop(columns=drop_cols).apply(pd.to_numeric, errors="coerce")

all_nan_cols = X_df.columns[X_df.isna().all()].tolist()
if all_nan_cols:
    print(f"Dropping {len(all_nan_cols)} all-NaN feature columns.")
    X_df = X_df.drop(columns=all_nan_cols)

feature_names = list(X_df.columns)
engineered_mask = np.array([is_engineered_feature(n) for n in feature_names], dtype=bool)
feature_family = np.array([engineered_family(n) for n in feature_names], dtype=object)
X_raw = X_df.to_numpy(dtype=np.float32)

valid_y = ~np.isnan(y)
if valid_y.sum() < len(y):
    print(f"Dropping {len(y) - valid_y.sum()} rows with missing {TARGET_COLUMN}.")
    X_raw = X_raw[valid_y]
    y = y[valid_y]

if N_FEATURES > X_raw.shape[1]:
    raise ValueError(f"N_FEATURES={N_FEATURES} exceeds available features ({X_raw.shape[1]}).")

feature_hash = feature_name_hash(feature_names)
n_eng_loaded = int(engineered_mask.sum())
n_regional_loaded = int((feature_family == "regional_qct").sum())
n_eng_extra_loaded = int((feature_family == "eng_extra").sum())
print(f"Rows: {X_raw.shape[0]} | Features: {X_raw.shape[1]} | Target: {TARGET_COLUMN}")
print(f"Model: {MODEL_NAME}")
print(
    f"Engineered features LOADED: {n_eng_loaded} "
    f"(regional_qct={n_regional_loaded}, eng.*={n_eng_extra_loaded})"
)
print(f"Feature hash: {feature_hash}")

# %% CELL 5 — Leakage-safe CV for one seed
def run_cv_for_seed(seed: int, n_features: int = N_FEATURES, n_estimators: int = N_ESTIMATORS):
    """Run leakage-safe fold-internal feature selection (always XGBoost)
    and champion-model evaluation (whatever build_champion() returns)."""
    print("\n" + "=" * 70)
    print(f"STARTING SEED {seed} | MODEL={MODEL_NAME}")
    print("=" * 70)

    kf = KFold(n_splits=N_SPLITS, shuffle=SHUFFLE, random_state=seed)
    fold_records = []
    selected_counter = Counter()

    for fold, (train_idx, val_idx) in enumerate(kf.split(X_raw), start=1):
        X_tr_raw, X_val_raw = X_raw[train_idx], X_raw[val_idx]
        y_tr, y_val = y[train_idx], y[val_idx]

        # Fit imputer on training fold only.
        X_tr, X_val, _ = median_impute_train_val(X_tr_raw, X_val_raw)

        # Winsorize outliers/sentinels using training fold only (see CELL 3 note).
        X_tr, X_val = clip_train_val(X_tr, X_val)

        # Fit selector on training fold only (always XGBoost — see CELL 3b note).
        selector = xgb.XGBRegressor(**selector_params(seed))
        selector.fit(X_tr, y_tr)
        top_idx = np.argsort(selector.feature_importances_)[::-1][:n_features]
        selected_names = [feature_names[i] for i in top_idx]
        selected_counter.update(selected_names)

        n_eng = int(engineered_mask[top_idx].sum())
        n_regional = int((feature_family[top_idx] == "regional_qct").sum())
        n_eng_extra = int((feature_family[top_idx] == "eng_extra").sum())

        # Train champion model (MODEL_NAME) on selected features only.
        model = build_champion(seed, n_estimators)
        model.fit(X_tr[:, top_idx], y_tr)
        preds = model.predict(X_val[:, top_idx])

        row = {
            "seed": seed,
            "fold": fold,
            "model": MODEL_NAME,
            "n_train": int(len(train_idx)),
            "n_val": int(len(val_idx)),
            "r2": float(r2_score(y_val, preds)),
            "spearman": stable_spearman(y_val, preds),
            "pearson": stable_pearson(y_val, preds),
            "mse": float(mean_squared_error(y_val, preds)),
            "eng_selected": n_eng,
            "regional_qct_selected": n_regional,
            "eng_extra_selected": n_eng_extra,
            "selected_feature_names": " | ".join(selected_names),
        }
        fold_records.append(row)
        print(
            f"Seed {seed}, fold {fold}: "
            f"R^2={row['r2']:.4f} | Spearman={row['spearman']:.4f} | "
            f"Pearson={row['pearson']:.4f} | MSE={row['mse']:.4f} | "
            f"eng={n_eng} (regional_qct={n_regional}, eng.*={n_eng_extra})"
        )

    return fold_records, selected_counter


# %% CELL 6 — Run all five seeds
all_fold_records = []
overall_selected_counter = Counter()

for seed in SEEDS:
    fold_records, selected_counter = run_cv_for_seed(seed)
    all_fold_records.extend(fold_records)
    overall_selected_counter.update(selected_counter)

fold_df = pd.DataFrame(all_fold_records)
fold_df.to_csv(os.path.join(OUTPUT_DIR, "all_seed_fold_results.csv"), index=False)

summary_rows = []
for seed, g in fold_df.groupby("seed", sort=True):
    summary_rows.append({
        "seed": int(seed),
        "mean_r2": float(g["r2"].mean()),
        "std_r2": float(g["r2"].std(ddof=0)),
        "mean_spearman": float(g["spearman"].mean()),
        "std_spearman": float(g["spearman"].std(ddof=0)),
        "mean_pearson": float(g["pearson"].mean()),
        "std_pearson": float(g["pearson"].std(ddof=0)),
        "mean_mse": float(g["mse"].mean()),
        "std_mse": float(g["mse"].std(ddof=0)),
        "mean_eng_selected": float(g["eng_selected"].mean()),
        "mean_regional_qct_selected": float(g["regional_qct_selected"].mean()),
        "mean_eng_extra_selected": float(g["eng_extra_selected"].mean()),
    })
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUTPUT_DIR, "per_seed_summary.csv"), index=False)

selected_counts_df = pd.DataFrame([
    {
        "feature_name": name,
        "number_of_seed_folds_selected": count,
        "out_of_total_seed_folds": len(SEEDS) * N_SPLITS,
    }
    for name, count in overall_selected_counter.most_common()
])
selected_counts_df.to_csv(os.path.join(OUTPUT_DIR, "selected_feature_counts.csv"), index=False)

# %% CELL 7 — Print five-seed CV report
print("\n" + "=" * 70)
print(f"PER-SEED SUMMARY (model={MODEL_NAME})")
print("=" * 70)
print(summary_df.to_string(index=False))

print("\n" + "=" * 70)
print(f"OVERALL BASELINE SUMMARY ACROSS 5 SEEDS (model={MODEL_NAME})")
print("=" * 70)
print(f"Mean R^2:       {summary_df['mean_r2'].mean():.6f} ± {summary_df['mean_r2'].std(ddof=0):.6f}")
print(f"Mean Spearman:  {summary_df['mean_spearman'].mean():.6f} ± {summary_df['mean_spearman'].std(ddof=0):.6f}")
print(f"Mean Pearson:   {summary_df['mean_pearson'].mean():.6f} ± {summary_df['mean_pearson'].std(ddof=0):.6f}")
print(f"Mean MSE:       {summary_df['mean_mse'].mean():.6f} ± {summary_df['mean_mse'].std(ddof=0):.6f}")

print("\n" + "=" * 70)
print("TOP SELECTED BASELINE FEATURES ACROSS ALL SEEDS")
print("=" * 70)
print(selected_counts_df.head(30).to_string(index=False))

# %% CELL 8 — Train final production model on all data using seed 42
# This is for deployment only. Do not report its in-sample fit as validation performance.
print(f"\nTraining final production selector/model ({MODEL_NAME}) on all data using seed 42...")
X_full, imputer_medians = median_impute_full(X_raw)
X_full, _ = clip_train_val(X_full, X_full)  # winsorize against its own percentiles

final_selector = xgb.XGBRegressor(**selector_params(FINAL_PRODUCTION_SEED))
final_selector.fit(X_full, y)

importance_df = pd.DataFrame({
    "feature_index": np.arange(len(feature_names)),
    "feature_name": feature_names,
    "importance_gain": final_selector.feature_importances_,
}).sort_values("importance_gain", ascending=False)

selected_feature_indices = importance_df.head(N_FEATURES)["feature_index"].to_numpy()
selected_feature_names = [feature_names[i] for i in selected_feature_indices]
X_selected = X_full[:, selected_feature_indices]

champion_model = build_champion(FINAL_PRODUCTION_SEED, N_ESTIMATORS)
champion_model.fit(X_selected, y)
print("Final production model trained.")

# %% CELL 9 — Save artifacts and metadata
joblib.dump(champion_model, os.path.join(OUTPUT_DIR, f"champion_{MODEL_NAME}_model.pkl"))
joblib.dump(final_selector, os.path.join(OUTPUT_DIR, "feature_selector_model.pkl"))
joblib.dump(selected_feature_indices, os.path.join(OUTPUT_DIR, "selected_feature_indices.pkl"))
np.save(os.path.join(OUTPUT_DIR, "imputer_medians.npy"), imputer_medians)

with open(os.path.join(OUTPUT_DIR, "full_feature_names_ordered.json"), "w") as f:
    json.dump(feature_names, f, indent=2)
with open(os.path.join(OUTPUT_DIR, "selected_feature_names_ordered.json"), "w") as f:
    json.dump(selected_feature_names, f, indent=2)

importance_df.to_csv(os.path.join(OUTPUT_DIR, "final_feature_importance_ranking.csv"), index=False)

metadata = {
    "experiment": f"five_seed_baseline_{MODEL_NAME}_top300",
    "project": "Enhancing sgRNA efficiency in CRISPR-Cas9 genome editing",
    "species": "Escherichia coli",
    "model_name": MODEL_NAME,
    "data_path": DATA_PATH,
    "target_column": TARGET_COLUMN,
    "dropped_columns": drop_cols,
    "all_nan_columns_dropped": all_nan_cols,
    "x_shape_full": list(X_raw.shape),
    "x_shape_selected_final_production": list(X_selected.shape),
    "seeds": SEEDS,
    "n_splits": N_SPLITS,
    "n_features": N_FEATURES,
    "n_estimators_selector": 300,
    "n_estimators_final_model": N_ESTIMATORS,
    "final_production_seed": FINAL_PRODUCTION_SEED,
    "feature_order_hash": feature_hash,
    "overall_mean_r2_across_seed_means": float(summary_df["mean_r2"].mean()),
    "overall_std_r2_across_seed_means": float(summary_df["mean_r2"].std(ddof=0)),
    "overall_mean_spearman_across_seed_means": float(summary_df["mean_spearman"].mean()),
    "overall_std_spearman_across_seed_means": float(summary_df["mean_spearman"].std(ddof=0)),
    "overall_mean_pearson_across_seed_means": float(summary_df["mean_pearson"].mean()),
    "overall_std_pearson_across_seed_means": float(summary_df["mean_pearson"].std(ddof=0)),
    "overall_mean_mse_across_seed_means": float(summary_df["mean_mse"].mean()),
    "overall_std_mse_across_seed_means": float(summary_df["mean_mse"].std(ddof=0)),
    "engineered_features_loaded": n_eng_loaded,
    "regional_qct_loaded": n_regional_loaded,
    "eng_extra_loaded": n_eng_extra_loaded,
    "note": (
        "Leakage-safe five-seed evaluation. Feature selector is always "
        "XGBoost, fit inside each CV training fold only; champion model is "
        "MODEL_NAME, also fit per fold on the selected features. The final "
        "production model is trained on all data and is not a validation score."
    ),
}
with open(os.path.join(OUTPUT_DIR, "model_metadata.json"), "w") as f:
    json.dump(metadata, f, indent=2)

print(f"\nSaved five-seed baseline outputs and production artifacts to: {OUTPUT_DIR}/")
print("\nTop 15 final production selected features:")
for i, name in enumerate(selected_feature_names[:15], start=1):
    print(f"{i:2d}. {name}")
