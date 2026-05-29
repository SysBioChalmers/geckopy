"""Compare max vs median aggregation across real BRENDA + ecTestGEM."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import cobra

sys.path.insert(0, "src")
from geckopy import ModelAdapter, make_ec_model, fuzzy_kcat_matching
from geckopy.databases import load_brenda_data, load_phyl_dist

ROOT = Path("/mnt/c/Work/GitHub/geckopy")
BRENDA_DIR = ROOT / "src" / "geckopy" / "data" / "brenda"
EXAMPLE = ROOT / "examples" / "ecTestGEM"

print("=" * 70)
print("Loading bundled BRENDA snapshot")
brenda = load_brenda_data(BRENDA_DIR)
kcat = brenda.kcat
print(f"  rows: {len(kcat)},  unique EC: {kcat['ec_code'].nunique()}")

# ---- Global per-EC stats: how much does median deviate from max? -----
groups = kcat[kcat["kcat"] > 0].groupby("ec_code")["kcat"]
records = []
for ec, vals in groups:
    v = vals.values
    if len(v) < 5:
        continue
    records.append({
        "ec": ec, "n": len(v),
        "max": float(np.max(v)),
        "median": float(np.median(v)),
        "p90": float(np.quantile(v, 0.9)),
    })
df = pd.DataFrame(records)
df["ratio_med_max"] = df["median"] / df["max"]
df["log10_max"] = np.log10(df["max"])
df["log10_med"] = np.log10(df["median"])
df["log_spread"] = df["log10_max"] - df["log10_med"]

print(f"\nGlobal stats over {len(df)} EC codes with >=5 measurements:")
print(f"  ratio median/max: quartiles "
      f"{df['ratio_med_max'].quantile([.25,.5,.75]).round(4).tolist()}")
print(f"  log10(max/median): mean {df['log_spread'].mean():.2f}, "
      f"median {df['log_spread'].median():.2f}, "
      f"p90 {df['log_spread'].quantile(0.9):.2f}")
print(f"  fraction with max >= 10x median: "
      f"{(df['log_spread'] >= 1).mean()*100:.1f}%")
print(f"  fraction with max >= 100x median: "
      f"{(df['log_spread'] >= 2).mean()*100:.1f}%")
print(f"  fraction with max >= 1000x median: "
      f"{(df['log_spread'] >= 3).mean()*100:.1f}%")

# ---- Per-reaction comparison on ecTestGEM -----
print("\n" + "=" * 70)
print("Fuzzy matching ecTestGEM (org='testus testus') against bundled BRENDA")
adapter = ModelAdapter.from_folder(EXAMPLE)
cobra_model = cobra.io.read_sbml_model(str(adapter.params.conv_gem))
ec_model = make_ec_model(cobra_model, adapter)
phyl = load_phyl_dist(adapter.get_phyl_dist_path())

print("  running aggregate='max' ...")
res_max = fuzzy_kcat_matching(ec_model, brenda, phyl, aggregate="max")
print("  running aggregate='median' ...")
res_med = fuzzy_kcat_matching(ec_model, brenda, phyl, aggregate="median")

cmp = res_max.merge(
    res_med, on="rxn_id", how="outer", suffixes=("_max", "_med"),
)
cmp = cmp[["rxn_id", "kcat_max", "kcat_med", "wildcard_level_max"]]
cmp = cmp[cmp["kcat_max"].notna() | cmp["kcat_med"].notna()]
cmp["ratio_med_max"] = cmp["kcat_med"] / cmp["kcat_max"]
print(cmp.to_string(index=False))

# ---- Sample real EC codes for per-reaction-style comparison -----
print("\n" + "=" * 70)
print("Sample: real EC codes (sorted by n_measurements desc, top 20)")
df_sorted = df.sort_values("n", ascending=False).head(20)
sample = df_sorted[["ec", "n", "max", "median", "p90", "ratio_med_max"]]
sample = sample.round({"max": 2, "median": 2, "p90": 2, "ratio_med_max": 4})
print(sample.to_string(index=False))
