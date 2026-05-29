"""Empirical comparison: max vs median over the **raw** BRENDA JSON
(skipping the snapshot's per-(ec, substrate, organism) max-collapse).

Two layers of aggregation get untangled here:

1. *Within-triple* -- for each (ec, substrate, organism) triple with
   multiple raw measurements, how much info does the snapshot's max
   discard? Compared to median, mean.
2. *Per-EC* -- aggregating triples across an EC (the question the
   previous analysis answered, but starting from raw rows instead of
   already-max-collapsed snapshot rows).

Run from the repo root: `python scripts/kcat_aggregation_raw.py`.
"""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from geckopy.databases.brenda.parse import parse_brenda_json

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "src" / "geckopy" / "data" / "brenda" / "_cache" / "brenda_2026_1.json"

# Trypsin (3.4.21.4): BRENDA values include many non-physiological mutants
# and unrelated substrates; treated as unreliable for this comparison.
EXCLUDED_ECS = {"3.4.21.4"}

print(f"Streaming raw BRENDA rows from {RAW.name} ...")
# Group raw measurements by (ec, substrate, organism) triple so we can
# look at the within-triple distribution before any aggregation.
triples: dict[tuple[str, str, str], list[float]] = defaultdict(list)
n_raw = 0
n_skipped_marker = 0
n_skipped_nonpos = 0
for row in parse_brenda_json(RAW):
    if row.kind != "kcat":
        continue
    if row.ec in EXCLUDED_ECS:
        continue
    n_raw += 1
    v = row.value
    if v <= 0:
        n_skipped_nonpos += 1
        continue
    triples[(row.ec, row.substrate, row.organism)].append(v)

print(
    f"  raw kcat rows (trypsin excluded): {n_raw:,}  "
    f"(of which non-positive dropped: {n_skipped_nonpos:,})"
)
print(f"  unique (ec, substrate, organism) triples: {len(triples):,}")

# -- Layer 1: within-triple ------------------------------------------------- #
# How much info does the snapshot's max-per-triple collapse hide?
triple_sizes = np.array([len(v) for v in triples.values()])
n_multi = int((triple_sizes >= 2).sum())
print(
    f"\nTriples with >= 2 raw measurements (within-triple aggregation matters): "
    f"{n_multi:,} ({100 * n_multi / len(triples):.1f}%)"
)

multi_rows = []
for (ec, sub, org), vals in triples.items():
    if len(vals) < 2:
        continue
    v = np.asarray(vals, dtype=float)
    multi_rows.append({
        "ec": ec, "substrate": sub, "organism": org,
        "n": len(v),
        "max": float(np.max(v)),
        "median": float(np.median(v)),
    })
multi = pd.DataFrame(multi_rows)
multi["log_spread"] = np.log10(multi["max"] / multi["median"])

qs = multi["log_spread"].quantile([0.25, 0.5, 0.75]).round(3).tolist()
print(
    f"  within-triple log10(max/median): quartiles {qs}  "
    f"mean {multi['log_spread'].mean():.3f}  "
    f"p90 {multi['log_spread'].quantile(0.9):.3f}"
)
print(
    f"  within-triple fraction max >= 10x median: "
    f"{100 * (multi['log_spread'] >= 1).mean():.1f}%; "
    f">= 100x: {100 * (multi['log_spread'] >= 2).mean():.1f}%"
)


# -- Layer 2: per-EC -------------------------------------------------------- #
# Collapse triples to one max and one median per (ec, sub, org), then
# aggregate across triples per EC. We track both:
#   "max-of-max"  (what the snapshot does today, then per-EC max -> picks the
#                  globally largest measurement)
#   "median-of-medians" (snapshot collapse = median, per-EC = median)
ec_to_max_rows: dict[str, list[float]] = defaultdict(list)
ec_to_median_rows: dict[str, list[float]] = defaultdict(list)
for (ec, _, _), vals in triples.items():
    v = np.asarray(vals, dtype=float)
    ec_to_max_rows[ec].append(float(np.max(v)))
    ec_to_median_rows[ec].append(float(np.median(v)))

records = []
for ec in ec_to_max_rows:
    max_rows = np.asarray(ec_to_max_rows[ec])
    med_rows = np.asarray(ec_to_median_rows[ec])
    if len(max_rows) < 5:
        continue
    records.append({
        "ec": ec,
        "n_triples": len(max_rows),
        # Snapshot=max, per-EC=max  (current behaviour)
        "kcat_maxmax": float(np.max(max_rows)),
        # Snapshot=max, per-EC=median  (today's "median" code path)
        "kcat_maxmed": float(np.median(max_rows)),
        # Snapshot=median, per-EC=median  (the consistent median path)
        "kcat_medmed": float(np.median(med_rows)),
        # Snapshot=median, per-EC=max  (mixed)
        "kcat_medmax": float(np.max(med_rows)),
    })

df = pd.DataFrame(records)
df["spread_maxmax_over_medmed"] = np.log10(df["kcat_maxmax"] / df["kcat_medmed"])
df["spread_maxmed_over_medmed"] = np.log10(df["kcat_maxmed"] / df["kcat_medmed"])
df["spread_maxmax_over_maxmed"] = np.log10(df["kcat_maxmax"] / df["kcat_maxmed"])

print(
    f"\nPer-EC stats over {len(df)} EC codes (>=5 triples, excluding "
    f"{sorted(EXCLUDED_ECS)}):"
)
def _q(col, label):
    qs = df[col].quantile([0.25, 0.5, 0.75]).round(3).tolist()
    print(
        f"  {label:<38}  quartiles {qs}  mean {df[col].mean():.3f}  "
        f"p90 {df[col].quantile(0.9):.3f}"
    )

_q("spread_maxmax_over_medmed", "log10(maxmax / medmed)")
_q("spread_maxmed_over_medmed", "log10(maxmed / medmed) -- today's flag")
_q("spread_maxmax_over_maxmed", "log10(maxmax / maxmed) -- snap-only swap")

for thresh, label in ((1, "10x"), (2, "100x"), (3, "1000x")):
    f_mm = 100 * (df["spread_maxmax_over_medmed"] >= thresh).mean()
    f_mt = 100 * (df["spread_maxmed_over_medmed"] >= thresh).mean()
    print(
        f"  maxmax/medmed >= {label}: {f_mm:5.1f}%   "
        f"maxmed/medmed >= {label}: {f_mt:5.1f}%"
    )

print("\nTop 20 best-measured ECs (by n_triples, excluding trypsin):")
top = df.sort_values("n_triples", ascending=False).head(20).copy()
out = top[
    ["ec", "n_triples", "kcat_maxmax", "kcat_maxmed", "kcat_medmed"]
].round({"kcat_maxmax": 2, "kcat_maxmed": 2, "kcat_medmed": 2})
out["maxmax/medmed"] = (out["kcat_maxmax"] / out["kcat_medmed"]).round(1)
out["maxmed/medmed"] = (out["kcat_maxmed"] / out["kcat_medmed"]).round(2)
print(out.to_string(index=False))

# Esterase 3.1.1.1 focus -- the user flagged this as suspicious.
ester = next((r for r in records if r["ec"] == "3.1.1.1"), None)
if ester is not None:
    print(
        "\nEsterase (3.1.1.1) focus:\n"
        f"  n_triples={ester['n_triples']}, kcat_maxmax={ester['kcat_maxmax']:.0f}, "
        f"kcat_maxmed={ester['kcat_maxmed']:.2f}, kcat_medmed={ester['kcat_medmed']:.2f}"
    )
