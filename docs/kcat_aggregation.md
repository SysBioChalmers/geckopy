# kcat aggregation: real-world impact of `max` vs `median`

Compiled by running the per-EC distribution across the bundled BRENDA snapshot
(`src/geckopy/data/brenda/max_kcat.tsv`, 39 318 rows, 3 567 unique EC codes).
For each EC code with ≥ 5 reported kcat measurements (1 723 codes), we computed
the max and the median.

## Global statistic

| | quartiles | mean | p90 |
|---|---|---|---|
| `median / max` ratio | 0.013 / 0.073 / 0.236 | – | – |
| `log10(max / median)` | – | 1.39 | 2.78 |

So for the typical (median) EC code, `max` is **~14×** the median; on average
**~25×**. The right tail is much heavier:

| fraction of EC codes where `max` ≥ … × `median` | |
|---|---|
| 10× | **57 %** |
| 100× | **22.5 %** |
| 1 000× | **7.9 %** |

I.e. for *more than half* of all enzymes BRENDA has data for, the current
default of taking the largest reported turnover inflates the kcat by an order
of magnitude or more relative to the median; for ~8 %, by *three* orders of
magnitude.

## Top 20 best-measured EC codes

Sorted by number of BRENDA measurements (so the statistics are least noisy).

| EC | n | max [1/s] | median [1/s] | max / median |
|---:|---:|---:|---:|---:|
| 3.5.2.6 | 601 | 7 150 | 20 | 358× |
| 3.1.1.1 | 355 | 431 000 | 29 | **14 862×** |
| 1.1.1.1 | 318 | 3 500 | 3.1 | 1 129× |
| 3.2.1.21 | 289 | 1 214 290 | 14 | **86 735×** |
| 1.10.3.2 | 269 | 158 300 | 42.1 | 3 760× |
| 3.4.21.4 | 265 | 910 000 | 5.1 | **178 431×** |
| 2.3.2.5 | 264 | 220 | 22.95 | 9.6× |
| 3.1.1.73 | 235 | 14 793 | 10 | 1 479× |
| 3.4.21.83 | 213 | 11 300 | 18 | 628× |
| 1.1.99.18 | 203 | 156 | 17 | 9.2× |
| 1.14.18.1 | 181 | 219 000 | 26.7 | 8 202× |
| 4.1.1.39 | 168 | 23.1 | 4.2 | 5.5× |
| 3.4.21.35 | 165 | 1 380 | 0.81 | 1 704× |
| 2.5.1.18 | 161 | 9 600 | 9 | 1 067× |
| 3.4.24.15 | 156 | 179 | 4.75 | 38× |
| 3.2.1.20 | 139 | 1 278 | 26.2 | 49× |
| 3.1.8.1 | 133 | 30 900 | 11.72 | 2 637× |
| 3.1.2.20 | 128 | 190 | 2.15 | 88× |
| 3.4.24.16 | 122 | 11.9 | 1.55 | 7.7× |
| 3.4.11.1 | 112 | 34 500 | 2.94 | **11 735×** |

The extremes — trypsin (3.4.21.4) at 910 000 /s vs median 5 /s, esterase
(3.1.1.1) at 431 000 /s vs median 29 /s — are characteristic of small, fast
hydrolases where engineered mutants and non-physiological substrates land in
BRENDA alongside the natural reaction. Taking `max` picks those outliers.

## Implication

Inside the ec-pipeline an inflated kcat shows up as a *deflated* protein
demand (stoichiometric coefficient on `prot_<id>` is `-MW / (kcat·3600)`), so
**the model systematically thinks each reaction needs less enzyme than it
actually does**. Switching the default to `median` would tighten the protein
budget for most reactions; for the ones with the most extreme spread (heavily
studied hydrolases / kinases / cytochromes) the change is 1–3 orders of
magnitude on the kcat and a commensurate increase in predicted enzyme demand.

## Recommendation

The flag is already in place (`fuzzy_kcat_matching(..., aggregate="median")`,
`fill_kcats_from_isozymes(..., aggregate="median")`); flipping the **default**
to `median` is a one-line change. The trade-off is reproducibility against
MATLAB GECKO results — anyone re-running an existing study would see different
numbers. A reasonable middle ground would be to flip the default for new
projects (so the more defensible choice is the out-of-box one) and document
how to pin `aggregate="max"` for MATLAB-compatible runs.

## Reproduce

The script is in `scripts/kcat_aggregation_compare.py` (also generated this
report); ~30 s to run against the bundled snapshot.
