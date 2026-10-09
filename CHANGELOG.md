# Changelog

## v1.0.1 (unreleased — tag after merging; Zenodo will assign a new version DOI)

### Added
- **Cross-stage repeat corroboration in Stage C** (`graph_depth_screen.py`): scans Stage A's all-vs-all self-alignment output for internal repeats and the GetOrganelle log for graph-disentangling warnings, flagging either even when the per-node depth table shows nothing anomalous (`--self-alignment-tsv`, `--getorganelle-log`, `--internal-repeat-min-len`, `--internal-repeat-min-pident`). Motivated by the *Potamopyrgus antipodarum* validation (SRR22947759).
- **Pre-flight seed check** (`check_seed.py`, `-S` to skip): warns, before Stage A, when a seed shares no exact k-mer with the reads.
- **Per-lineage summary** (`summarize_lineages.py` → `report/lineage_summary.tsv`) with a recommendation, including detection of Stage B collapsing onto the sibling lineage (*Perumytilus purpuratus* validation, SRR7504397).
- `-V` prints the version; the version is also written to `run.log`.

### Changed (breaking)
- **Animal mitochondrial genomes only.** `-a` now accepts only `animal` (now the default; the previous default is gone); NOVOPlasty is always run with `Type = mito` and GetOrganelle with `-F animal_mt`. `make_novoplasty_config.py --organism-type` likewise.
- `graph_depth_screen.py`: removed the optional external-genome screen together with its two `repeat_candidates.tsv` columns. The remaining columns and their order are unchanged.

### Fixed
- In dual-lineage mode (`-M`), a lineage whose Stage A recovered nothing no longer aborts the whole run under `set -e`; Stage B/C is attempted for each lineage independently and the run ends normally (exit status 2 if any lineage failed).
- A failure of `graph_depth_screen.py` is now reported explicitly instead of being silently ignored/aborting.

### Documentation
- README: seed-choice guidance, pre-flight check, per-lineage summary and exit status.

## v1.0.0
First archived release (Zenodo DOI 10.5281/zenodo.23039343).
