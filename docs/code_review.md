# geckopy code review — open items

One decision left.

## Default kcat aggregation: flip `max` / `mean` → `median`?

The flags are wired (`fuzzy_kcat_matching(aggregate=…)`,
`fill_kcats_from_isozymes(aggregate=…)`); only the **shipped defaults** are
open. Empirical comparison on the bundled BRENDA snapshot
([`docs/kcat_aggregation.md`](kcat_aggregation.md)):

- 1 723 EC codes with ≥ 5 BRENDA measurements; median ratio `median/max` =
  **0.073**, i.e. `max` is **~14×** the median for the typical reaction.
- **57 %** of EC codes have `max ≥ 10×` the median; **22.5 %** ≥ 100×;
  **7.9 %** ≥ 1 000×.
- Heavily-studied hydrolases / kinases / cytochromes are 10⁴–10⁵× spread
  (trypsin 3.4.21.4: max 910 000, median 5.1; esterase 3.1.1.1: max 431 000,
  median 29).

Choosing `max` picks engineered mutants / non-physiological substrates and
inflates kcats by ≥ 1 order of magnitude for over half of all reactions.
Downstream this *deflates* predicted enzyme demand commensurately.

The decision: flip the default to `median` (more defensible scientifically,
breaks MATLAB-GECKO byte-for-byte reproducibility) or keep `max` (faithful
port, biased). A middle option: flip the default but document
`aggregate="max"` for MATLAB-compatible re-runs.

Also still open: extending the flag to (a) `brenda/parse.py`'s range collapse
(always upper bound today) and (b) the default of `apply_kcat_list`
(`criteria="max"` today, already supports `"median"`).
