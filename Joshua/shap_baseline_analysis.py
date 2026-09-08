"""Explain saved baseline models on their validation folds; never train a model.

Usage: python shap_baseline_analysis.py [baseline_directory] [--dataset CSV]

Required artifact contract (paths relative to baseline_directory):
  metadata.json: existing baseline metadata plus dataset_sha256 (input CSV bytes)
  feature_names.json: original feature names in model input order
  models/seed_41_fold_1.json: final_model.save_model(...) for each seed/fold
  folds/seed_41_fold_1.json: {
    "seed": 41, "fold": 1, "selected_feature_indices": [...300...],
    "selected_feature_names": [...300 in model input order...],
    "training_medians": [...one per ORIGINAL feature, NaN medians replaced by 0...]
  }
Optional validation_indices in each fold JSON must match the reconstructed split.
These files must be saved during the original baseline run, not guessed afterward.
Dependencies: numpy pandas scipy scikit-learn xgboost shap matplotlib.
"""
import argparse
import hashlib
import json
import re
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

SEEDS = [41, 42, 43, 44, 45]
METRICS = ["mean_abs_shap", "mean_signed_shap", "median_abs_shap", "std_abs_shap"]


def read_json(path):
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def names_hash(names):
    return hashlib.sha256(json.dumps(names, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest()


def preflight(root):
    """Check every saved artifact before loading data or computing explanations."""
    missing = []
    for path in [root / "metadata.json", root / "feature_names.json"]:
        if not path.is_file():
            missing.append(str(path))
    for seed in SEEDS:
        for fold in range(1, 6):
            for folder in ("models", "folds"):
                path = root / folder / f"seed_{seed}_fold_{fold}.json"
                if not path.is_file():
                    missing.append(str(path))
    if missing:
        raise ValueError(
            "Missing required saved baseline artifacts:\n  " + "\n  ".join(missing)
            + "\nThe baseline script must save each FINAL trained model using model.save_model(), "
            "its selected indices/names in model input order, and the original training-fold "
            "medians (after replacing undefined medians with zero). Save original feature ordering "
            "in feature_names.json and the input CSV SHA-256 in metadata.json['dataset_sha256']. "
            "See the artifact schema in this script's opening docstring. Score tables and selected "
            "names alone cannot recover trained models or medians. This analysis will not retrain."
        )
    meta, names = read_json(root / "metadata.json"), read_json(root / "feature_names.json")
    expected = {"seeds": SEEDS, "n_splits": 5, "n_features_selected": 300,
                "n_estimators_selector": 300, "n_estimators_final_model": 400}
    for key, value in expected.items():
        if meta.get(key) != value:
            raise ValueError(f"Baseline metadata {key}: expected {value!r}, got {meta.get(key)!r}.")
    if not meta.get("dataset_sha256"):
        raise ValueError("Missing metadata.dataset_sha256: save the SHA-256 of the original training CSV during the baseline run to verify data and row ordering.")
    if not isinstance(names, list) or len(set(names)) != len(names):
        raise ValueError("feature_names.json must be an ordered list of unique original feature names.")
    if meta.get("feature_hash") != names_hash(names):
        raise ValueError("Saved feature ordering does not match metadata.feature_hash.")
    artifacts = {}
    for seed in SEEDS:
        for fold in range(1, 6):
            path = root / "folds" / f"seed_{seed}_fold_{fold}.json"
            item = read_json(path)
            for key in ("seed", "fold", "selected_feature_indices", "selected_feature_names", "training_medians"):
                if key not in item:
                    raise ValueError(f"{path} is missing {key}.")
            indices = item["selected_feature_indices"]
            if item["seed"] != seed or item["fold"] != fold:
                raise ValueError(f"Seed/fold mismatch in {path}.")
            if len(indices) != 300 or len(set(indices)) != 300 or any(type(i) is not int or not 0 <= i < len(names) for i in indices):
                raise ValueError(f"{path}: expected 300 unique valid selected indices.")
            if item["selected_feature_names"] != [names[i] for i in indices]:
                raise ValueError(f"Selected feature names/order disagree with original indices in {path}.")
            if len(item["training_medians"]) != len(names):
                raise ValueError(f"{path}: training_medians must contain one value per original feature.")
            artifacts[(seed, fold)] = item
    return meta, names, artifacts


def position_of(name):
    # Only an explicit pN prefix is interpreted as position; no guessed k-mer spans.
    match = re.match(r"^p(\d+)(?=[A-Za-z_.]|$)", name, re.I)
    return int(match.group(1)) if match and 1 <= int(match.group(1)) <= 20 else None


def category_of(name):
    lower = name.lower()
    if re.match(r"^p\d+(basepair|monomer|dimer|trimer|tetramer)\.", lower):
        return "QCT / quantum"
    if any(token in lower for token in ("temperature", "structure", "melting", "temp", "mfe", "deltag")):
        return "thermodynamic"
    if re.match(r"^[ACGTU]{2,}sgRNA\.raw$", name):
        return "k-mer"
    if re.match(r"^p\d+[._]?[ACGTU](?:[._]|raw|$)", name):
        return "one-hot / nucleotide identity"
    if position_of(name) is not None or "distance" in lower:
        return "positional"
    return "other"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("baseline_directory", nargs="?", default="five_seed_baseline_xgboost_top300")
    parser.add_argument("--dataset", default="ecoli_feature_matrix.csv")
    parser.add_argument("--output", default="shap_baseline_analysis")
    parser.add_argument("--save-raw", action="store_true", help="Keep compressed validation SHAP arrays for each fold.")
    args = parser.parse_args()
    root, dataset, output = Path(args.baseline_directory), Path(args.dataset), Path(args.output)
    print(f"Checking saved baseline artifacts in {root.resolve()}...", flush=True)
    meta, original_names, artifacts = preflight(root)
    if sha256_file(dataset) != meta["dataset_sha256"]:
        raise ValueError("Dataset SHA-256 mismatch: use the exact CSV used to train the saved baseline models.")
    import numpy as np
    import pandas as pd
    import shap
    import xgboost as xgb
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.model_selection import KFold

    df = pd.read_csv(dataset, header=1, low_memory=False)
    y_valid = pd.to_numeric(df["cut.score"], errors="coerce").notna()
    dropped = meta.get("dropped_columns", ["cut.score", "sgRNAID"])
    if "cut.score" not in dropped or "sgRNAID" not in dropped:
        raise ValueError("Baseline dropped_columns must exclude cut.score and sgRNAID.")
    numeric = df.drop(columns=dropped).apply(pd.to_numeric, errors="coerce")
    numeric = numeric.drop(columns=numeric.columns[numeric.isna().all()])
    numeric = numeric.loc[y_valid].reset_index(drop=True)
    if list(numeric.columns) != original_names or numeric.shape != (meta["n_samples"], meta["n_original_features"]):
        raise ValueError("Reconstructed original feature order or data shape differs from saved baseline.")
    X = numeric.to_numpy(dtype=np.float32)
    del numeric, df
    output.mkdir(parents=True, exist_ok=True)
    rows, summaries = [], []
    with tempfile.TemporaryDirectory(prefix="shap_validation_") as temporary:
        cache = output / "raw_validation_shap" if args.save_raw else Path(temporary)
        cache.mkdir(parents=True, exist_ok=True)
        for seed in SEEDS:
            for fold, (_, val) in enumerate(KFold(5, shuffle=True, random_state=seed).split(X), 1):
                print(f"Explaining saved model: seed {seed}, fold {fold} (validation only)...", flush=True)
                item = artifacts[(seed, fold)]
                if "validation_indices" in item and item["validation_indices"] != val.tolist():
                    raise ValueError(f"Saved validation indices disagree for seed {seed}, fold {fold}.")
                indices = np.asarray(item["selected_feature_indices"])
                medians = np.asarray(item["training_medians"], dtype=np.float32)
                if not np.isfinite(medians).all():
                    raise ValueError(f"Nonfinite saved imputation medians in seed {seed}, fold {fold}.")
                values = X[np.ix_(val, indices)].copy()
                r, c = np.where(np.isnan(values))
                values[r, c] = medians[indices[c]]
                if not np.isfinite(values).all():
                    raise ValueError("Validation features contain infinity; cannot silently alter baseline preprocessing.")
                model = xgb.XGBRegressor()
                model.load_model(root / "models" / f"seed_{seed}_fold_{fold}.json")
                booster = model.get_booster()
                if booster.num_features() != 300 or booster.num_boosted_rounds() != 400:
                    raise ValueError(f"Unexpected saved model dimensions/rounds: seed {seed}, fold {fold}.")
                names = item["selected_feature_names"]
                if booster.feature_names is not None and booster.feature_names != names:
                    raise ValueError("Saved model feature names do not match fold metadata order.")
                explainer = shap.TreeExplainer(model, feature_perturbation="tree_path_dependent", model_output="raw")
                sv = np.asarray(explainer.shap_values(values, check_additivity=True))
                if sv.shape != values.shape:
                    raise ValueError(f"Unexpected SHAP shape {sv.shape}; expected {values.shape}.")
                predictions = booster.predict(xgb.DMatrix(values, feature_names=booster.feature_names), output_margin=True)
                np.testing.assert_allclose(sv.sum(axis=1) + np.asarray(explainer.expected_value).item(), predictions, rtol=1e-4, atol=1e-3)
                absolute = np.abs(sv)
                importance = absolute.mean(axis=0)
                ranks = np.argsort(np.argsort(-importance, kind="stable"), kind="stable") + 1
                for j, name in enumerate(names):
                    rows.append(dict(seed=seed, fold=fold, feature=name, mean_abs_shap=float(importance[j]),
                        mean_signed_shap=float(sv[:, j].mean()), median_abs_shap=float(np.median(absolute[:, j])),
                        std_abs_shap=float(absolute[:, j].std()), selection_rank=j+1,
                        shap_rank=int(ranks[j])))
                best = int(np.argmax(importance))
                summaries.append(dict(seed=seed, fold=fold, number_of_selected_features=300,
                    mean_abs_shap=float(absolute.mean()), top_feature=names[best], top_feature_abs_shap=float(importance[best])))
                np.savez_compressed(cache / f"seed_{seed}_fold_{fold}.npz", names=np.asarray(names),
                                    shap_values=sv, feature_values=values, validation_indices=val)

        importance_frame = pd.DataFrame(rows)
        stable = importance_frame.groupby("feature", sort=False).agg(
            selection_count=("fold", "size"), mean_abs_shap_when_selected=("mean_abs_shap", "mean"),
            std_abs_shap_when_selected=("mean_abs_shap", lambda s: s.std(ddof=0)),
            mean_signed_shap_when_selected=("mean_signed_shap", "mean"),
            top10_count=("shap_rank", lambda s: int((s <= 10).sum())),
            top20_count=("shap_rank", lambda s: int((s <= 20).sum())),
            top50_count=("shap_rank", lambda s: int((s <= 50).sum()))).reindex(original_names)
        for column in ("selection_count", "top10_count", "top20_count", "top50_count"):
            stable[column] = stable[column].fillna(0).astype(int)
        stable["selection_frequency"] = stable.selection_count / 25
        stable.index.name = "feature"
        stable = stable.reset_index()
        ordered = stable[stable.selection_count > 0].sort_values("mean_abs_shap_when_selected", ascending=False)
        top = ordered.head(50).rename(columns={"mean_abs_shap_when_selected": "mean_abs_shap", "std_abs_shap_when_selected": "std_abs_shap", "mean_signed_shap_when_selected": "mean_signed_shap"})
        top.insert(0, "rank", range(1, len(top)+1))
        stable.to_csv(output / "feature_stability.csv", index=False)
        importance_frame.drop(columns="shap_rank").to_csv(output / "fold_feature_importance.csv", index=False)
        top.to_csv(output / "top_50_shap_features.csv", index=False)
        pd.DataFrame(summaries).to_csv(output / "seed_fold_shap_summary.csv", index=False)
        importance_frame["position"] = importance_frame.feature.map(position_of)
        importance_frame["category"] = importance_frame.feature.map(category_of)

        def aggregate(column, labels):
            records = []
            for label in labels:
                group = importance_frame[importance_frame[column] == label]
                records.append({column: label, "total_abs_shap": group.mean_abs_shap.sum()/25,
                    "mean_abs_shap": group.mean_abs_shap.mean() if len(group) else float("nan"),
                    "selected_feature_count": group.feature.nunique(),
                    "selection_frequency": len(group[["seed", "fold"]].drop_duplicates())/25})
            return pd.DataFrame(records)

        positions = aggregate("position", range(1, 21))
        categories = aggregate("category", sorted(importance_frame.category.unique()))
        categories["proportion_total_shap"] = categories.total_abs_shap / categories.total_abs_shap.sum()
        positions.to_csv(output / "position_shap_summary.csv", index=False)
        categories.to_csv(output / "feature_category_shap_summary.csv", index=False)

        def bar(labels, values, filename, xlabel, horizontal=False):
            fig, ax = plt.subplots(figsize=(11, max(5, len(labels)*.22) if horizontal else 5))
            if horizontal:
                ax.barh(list(labels)[::-1], list(values)[::-1]); ax.set_xlabel(xlabel)
            else:
                ax.bar([str(x) for x in labels], values); ax.set_ylabel(xlabel)
                ax.tick_params(axis="x", rotation=30)
            fig.tight_layout(); fig.savefig(output / filename, dpi=180); plt.close(fig)
        bar(positions.position, positions.mean_abs_shap, "shap_by_position.png", "Mean |SHAP| per selected feature-fold")
        bar(categories.category, categories.total_abs_shap, "shap_by_feature_category.png", "Mean per-fold total |SHAP|")
        bar(top.feature, top.mean_abs_shap, "shap_global_summary.png", "Mean |SHAP| when selected", True)

        stable_top = ordered.sort_values(["selection_count", "mean_abs_shap_when_selected"], ascending=False).head(20)
        samples = {name: [[], []] for name in stable_top.feature}
        for path in sorted(cache.glob("*.npz")):
            with np.load(path, allow_pickle=False) as data:
                for j, name in enumerate(data["names"]):
                    if name in samples:
                        samples[name][0].append(data["feature_values"][:, j])
                        samples[name][1].append(data["shap_values"][:, j])
        directions = []
        for rank, (name, (xs, ys)) in enumerate(samples.items(), 1):
            x, s = np.concatenate(xs), np.concatenate(ys)
            correlation = float(np.corrcoef(x, s)[0, 1]) if x.std() > 0 and s.std() > 0 else float("nan")
            directions.append(dict(feature=name, mean_signed_shap=float(s.mean()),
                shap_q05=float(np.quantile(s, .05)), shap_q25=float(np.quantile(s, .25)),
                shap_q50=float(np.quantile(s, .5)), shap_q75=float(np.quantile(s, .75)),
                shap_q95=float(np.quantile(s, .95)), value_shap_correlation=correlation))
            if rank <= 10:
                take = np.random.default_rng(42).choice(len(x), min(10000, len(x)), replace=False)
                fig, ax = plt.subplots(figsize=(8, 5)); ax.scatter(x[take], s[take], s=5, alpha=.2)
                ax.set(xlabel=name, ylabel="Validation SHAP contribution", title="Model association (not a causal effect)")
                fig.tight_layout(); fig.savefig(output / f"dependence_{rank:02d}.png", dpi=180); plt.close(fig)
        pd.DataFrame(directions).to_csv(output / "feature_direction_summary.csv", index=False)

    report = ["="*60, "BASELINE XGBOOST SHAP ANALYSIS", "="*60,
        "Number of seeds: 5\nSeeds: 41, 42, 43, 44, 45\nFolds per seed: 5\nTotal seed-fold models: 25\nFeatures selected per fold: 300",
        "TOP 30 MOST STABLE / IMPORTANT FEATURES", ordered.sort_values(["selection_count", "mean_abs_shap_when_selected"], ascending=False).head(30).to_string(index=False),
        "TOP POSITIONS BY SHAP IMPORTANCE", positions.sort_values("total_abs_shap", ascending=False).to_string(index=False),
        "FEATURE CATEGORY SHAP SUMMARY", categories.to_string(index=False),
        "PAM-PROXIMAL POSITION SUMMARY", positions[positions.position >= 16].to_string(index=False)]
    leaders = stable_top.head(5)
    report.append("Frequently selected contributors: " + "; ".join(f"{r.feature} ({r.selection_count}/25 folds)" for r in leaders.itertuples()) + ".")
    position_total = positions.total_abs_shap.sum()
    if position_total > 0:
        fraction = positions.loc[positions.position >= 18, "total_abs_shap"].sum()/position_total
        report.append(f"Positions p18-p20 account for {fraction:.1%} of contribution assigned to explicit p1-p20 prefixes. "
                      "This describes relative model contribution; it does not establish biological importance.")
    report.append("Category proportions compare QCT with features whose names support sequence/thermodynamic labels; unclassified features remain 'other'. "
                  "Selection and top-50 counts describe repetition across the 25 models. Repeated CV reuses samples, so these are not independent replications. "
                  "Signed SHAP and dependence plots describe model associations, not causation.")
    print("\n\n".join(report))
    (output / "report.txt").write_text("\n\n".join(report), encoding="utf-8")
    analysis_meta = dict(input_baseline_directory=str(root.resolve()), dataset=str(dataset.resolve()), seeds=SEEDS,
        folds=5, number_of_selected_features=300, number_of_models_analyzed=25,
        shap_method="TreeExplainer, tree_path_dependent, raw output; saved model paths as background",
        validation_only=True, timestamp=datetime.now(timezone.utc).isoformat(), feature_hash=names_hash(original_names),
        dataset_sha256=meta["dataset_sha256"], raw_arrays_saved=args.save_raw,
        selection_rank="Saved selector/model input order; shap_rank used internally for top counts",
        aggregation="Feature SHAP means/std over selected folds only. Unselected features have NaN SHAP statistics. std uses ddof=0.",
        position_rule="Explicit leading p1-p20 only; oligomers assigned to encoded start position, without expanding span.",
        group_aggregation="total_abs_shap: sum feature-fold mean absolute SHAP /25; mean_abs_shap: mean over selected feature-folds; selected_feature_count: distinct selected features; selection_frequency: fraction of folds with any member selected.",
        direction_analysis="All validation observations where selected, including repeated observations across seeds; plots subsample at most 10000 points.",
        versions=dict(numpy=np.__version__, pandas=pd.__version__, shap=shap.__version__, xgboost=xgb.__version__))
    (output / "metadata.json").write_text(json.dumps(analysis_meta, indent=2), encoding="utf-8")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError, ImportError) as error:
        print(f"SHAP ANALYSIS STOPPED: {error}", file=sys.stderr)
        sys.exit(1)
