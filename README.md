# repeatMito

**Iterative reseeding and graph-based detection of repeat-mediated mitochondrial genome assembly breaks, from short reads alone — with a DUI dual-lineage mode for animal mitochondrial genomes.**

repeatMito combines two widely used short-read organelle assemblers — [NOVOPlasty](https://github.com/ndierckx/NOVOPlasty) and [GetOrganelle](https://github.com/Kinggerm/GetOrganelle) — in an iterative, cross-validated workflow, with [Bandage](https://github.com/rrwick/Bandage)-based assembly-graph depth screening to flag candidate repeat units. It is scoped to two problems:

1. **Repeat-mediated assembly breaks.** Large repeated sequences that a short-read assembly either collapses or fragments around, producing a *multipartite* mitochondrial genome.
2. **Doubly Uniparental Inheritance (DUI).** Recovering two divergent, co-occurring F-type/M-type mitochondrial lineages from the same heteroplasmic read pool (typically male gonad tissue in bivalves) without letting them cross-contaminate or merge into a chimeric consensus.

repeatMito is scoped to mitochondrial genome assembly. For assemblies that also need to resolve other organelle genome compartments alongside the mitochondrial genome, see the companion pipeline **duOrganelle**; use repeatMito when the problem is a repeat-mediated genome break and/or DUI.

Authors: Fernanda E. Angulo, José J. Núñez (Universidad Austral de Chile).

## Installation

```bash
git clone https://github.com/LESS-UACh/repeatMito.git
cd repeatMito
chmod +x repeatmito.sh
```

repeatMito is a wrapper around the tools listed under [Requirements](#requirements) below — install those first and make sure each is callable on `PATH` (`novoplasty.pl`/`NOVOPlasty4.3.1.pl`, `get_organelle_from_reads.py`, `blastn`, `minimap2`, `samtools`, `python3`, and optionally `Bandage`). No separate build step is required for repeatMito itself.

## What it does

Short reads cannot span large repeats the way long reads can, so a naive short-read assembly of a repeat-mediated multipartite genome tends to either collapse the repeat (hiding alternative conformations) or fragment into multiple disconnected contigs around it. repeatMito does not try to physically phase every alternative conformation. Instead it:

1. **Iteratively reseeds NOVOPlasty** to greedily recover as much non-redundant mitochondrial sequence as possible across several rounds (Stage A).
2. **Feeds everything recovered into GetOrganelle** as a combined multi-seed, letting it do an independent, read-based SPAdes assembly and produce a proper assembly graph (Stage B).
3. **Screens the resulting graph** for nodes with anomalous coverage depth or disconnection from the main graph component — candidate repeat units mediating the multipartite structure (Stage C).

In dual-lineage (DUI) mode, this same Stage A/B/C pipeline runs twice, fully independently, once per lineage — see below.

## Pipeline stages

### Stage A — Iterative NOVOPlasty reseeding
An initial NOVOPlasty assembly is seeded with a single conserved mitochondrial gene. All-vs-all BLASTn detects redundant/overlapping contigs among the round's output; these are reoriented to a common strand and reconciled into a majority-rule consensus (minimap2 + samtools pileup). That consensus (or, in later rounds, whatever new sequence was recovered) is supplied back to NOVOPlasty as an extension seed. This repeats until a round recovers less novel sequence than `MIN_NOVEL_BP` (default 500 bp).

### Stage B — GetOrganelle multi-seed assembly
Every non-redundant segment accumulated across all of Stage A's rounds is consolidated into one multi-sequence FASTA and supplied to GetOrganelle as a combined seed, guiding SPAdes-based read recruitment and de novo assembly. This is a genuinely independent assembly step — GetOrganelle re-assembles from the raw reads, using the accumulated segments only to recruit the right reads, so it is not simply echoing back whatever Stage A produced.

### Stage C — Graph depth screening
Per-node length and read-depth are parsed from the resulting GFA. Nodes disconnected from the main graph component and/or with coverage substantially exceeding the graph median are flagged as repeat-unit candidates, and cross-checked against the accumulated segment set to separate previously-seen sequence from newly-resolved content.

## Background: why the repeat-detection screen exists

Redundancy/consensus clustering in Stage A is designed to merge genuinely overlapping contigs into a single consensus. But a large genomic repeat (for example an inverted repeat) can, under that same logic, be collapsed into an artificially short "core" sequence, because contig-level clustering alone cannot always distinguish true redundancy from a real repeat structure with multiple valid conformations. That risk motivated Stage C's independent graph-depth/disconnection screening: re-assembling with GetOrganelle and checking the resulting graph for node depth or connectivity anomalies catches this kind of over-collapse, which the redundancy/consensus step alone cannot reliably detect. It is why repeatMito's checklist below (GC content, BLAST against known references, coordinate-level gene checks) is worth running on any "novel segment" before trusting it, whatever organism you're assembling.

## Dual-lineage mode (DUI: Doubly Uniparental Inheritance)

Some animal lineages — most studied cases are bivalve molluscs — exhibit Doubly Uniparental Inheritance: individuals (typically males) carry **two distinct mitochondrial genomes**, an F-type (maternally transmitted, present in somatic tissue and female gonad) and an M-type (paternally transmitted, present chiefly in male gonad tissue). F and M mitogenomes are homologous but often substantially diverged from one another, even in otherwise-conserved protein-coding genes. Sequencing male gonad tissue recovers reads from **both** genomes at once.

This is a fundamentally different problem from the repeat-mediated case above. There, the goal is to correctly **merge** redundant/overlapping contigs into one consensus. Here, the goal is the opposite: **keep two genuinely different genomes from being merged into a chimeric consensus.** If F-type and M-type contigs happen to share enough sequence identity in a conserved region to pass the ordinary redundancy-clustering threshold, the majority-consensus step would blend them into a biologically meaningless hybrid.

**Dual-lineage mode does not replace repeat detection — it runs alongside it, independently per lineage.** Each of the two tracks (F and M) is the *same* Stage A/B/C pipeline described above, run in full: Stage A's redundancy clustering and Stage C's graph-depth screening (`flag_repeat_candidate`, based on coverage fold-over-median and disconnection from the main graph component) execute exactly as in single-lineage mode, on each lineage's own assembly graph. So if the F-type or M-type mitogenome itself has repeat-mediated structure (e.g. a duplicated control region, documented in some animal groups), it gets flagged within that lineage's own output, at the same time the cross-lineage screen is keeping F and M from contaminating each other. The two kinds of screening are orthogonal and both always run in dual-lineage mode.

### How dual-lineage mode works

Enable it with `-M <M_type_seed.fasta>` (requires `-a animal`). `-s` becomes the F-type seed.

1. **Two fully independent tracks.** F-type and M-type each get their own complete Stage A iterative-reseeding loop, their own accumulated segment set, and their own Stage B/C run, under `outdir/F/` and `outdir/M/`. Raw contigs and consensus sequences from the two tracks are **never** pooled or all-vs-all clustered together — the structural separation (different seed, different NOVOPlasty run, different clustering pass) is what actually prevents F/M chimeras, not a post-hoc filter.
2. **Cross-lineage divergence screen**, at three layers (raw contig, cluster consensus, graph node): strong, *diverged* homology to the sibling lineage is the **expected**, unflagged outcome (that's just F and M being homologous). Only **near-identity** (`--cross-lineage-max-identity`, default 98%) over a substantial fraction of the sequence (`--cross-lineage-min-coverage-frac`, default 0.5) is flagged — that combination signals index-hopping, a nuclear mitochondrial pseudogene (NUMT), or reads cross-recruited between the two haplotypes during assembly, not genuine independent divergence. The screen uses a more sensitive BLASTn (`word_size 11` vs. `20` elsewhere) so it can actually detect diverged homology, not just near-identical matches.
3. **Every sequence gets a divergence report, flagged or not** (`cross_lineage_divergence_report.tsv` per round, and `pident_vs_sibling_lineage`/`pct_matches_sibling_lineage` columns in the final `repeat_candidates.tsv`), so you can directly check the observed F/M divergence against what's published for the taxon — not just see a pass/fail flag.
4. **Rounds interleave F and M by round number.** Because the two tracks run independently, each track's cross-lineage screen compares against the sibling's most recently *completed* round — so round 1 of whichever track runs first has nothing to screen against yet (nothing accumulated on the sibling side), and from round 2 onward both tracks compare against genuinely current sibling data. This is logged explicitly and is a known, intentional limitation rather than a bug.

### Usage

```bash
./repeatmito.sh -m paired \
  -1 reads_R1.fq.gz -2 reads_R2.fq.gz \
  -a animal \
  -s F_type_reference_or_gene.fasta \
  -M M_type_reference_or_gene.fasta \
  -o output_dir \
  [-t threads] [-r max_rounds] [-x cross_lineage_max_identity] [-y cross_lineage_min_coverage_frac]
```

| Flag | Meaning | Default |
|---|---|---|
| `-s` | F-type seed (published F-type reference, own species or close relative) | — |
| `-M` | M-type seed (published M-type reference). Presence enables dual-lineage mode. | none (single-lineage mode) |
| `-x` | Cross-lineage screen: percent identity to the sibling lineage at/above which a sequence is flagged as suspect contamination rather than real divergence | `98.0` |
| `-y` | Cross-lineage screen: minimum fraction of sequence length covered by near-identical sibling matches to trigger the flag | `0.5` |

All other flags (`-t`, `-r`, `-g`, etc.) behave as in single-lineage mode and apply to both tracks.

### Output layout (dual-lineage mode)

```
output_dir/
├── F/                                   Full Stage A/B/C output tree for the F-type track
│   ├── consensus/
│   │   ├── round_N/
│   │   │   ├── cross_lineage_screen_report_RAW_contigs.tsv   Pre-clustering screen (per raw contig)
│   │   │   ├── flagged_cross_lineage_raw_contigs.fasta        Raw contigs excluded before clustering
│   │   │   ├── cross_lineage_divergence_report.tsv            Post-consensus screen (identity/coverage, ALL sequences)
│   │   │   ├── flagged_cross_lineage_contamination.fasta      Cluster consensuses excluded after building
│   │   │   └── novel_segment.fasta
│   │   └── accumulated_nonredundant.fasta
│   ├── getorganelle/
│   └── bandage/
│       └── repeat_candidates.tsv        Same repeat-candidate columns as single-lineage mode (fold_over_median,
│                                         flag_repeat_candidate, …) PLUS pident_vs_sibling_lineage /
│                                         flag_suspect_cross_lineage_contamination -- both screens run together
├── M/                                   Same structure, M-type track (screened against F's accumulated set)
└── report/run.log                       Combined log for both tracks
```

### Before trusting a lineage assignment

Neither the structural separation nor the identity screen is a substitute for checking the result. Before treating a sequence as confirmed F-type or M-type:

1. **Check the divergence report itself**, not just the flag column. A sequence sitting just under the `-x` threshold is still worth a manual look, especially early on when you don't yet know what typical F/M divergence looks like for this species/gene.
2. **BLAST each assembled lineage against its own published reference and against the sibling lineage's reference.** A genuine F-type sequence should be clearly closer to F references (across related taxa) than to M references, and vice versa — not just "not near-identical to the sibling track's accumulated set" (which only catches contamination picked up *during this pipeline run*, not a sequence that was mis-seeded from the start).
3. **Watch relative depth of coverage between the two tracks.** In DUI, M-type mtDNA is often present at lower copy number than F-type even in male gonad tissue (though this varies by species and can invert). A track that assembles suspiciously cleanly at F-type-like depth when seeded with the M-type reference is worth a second look.
4. Coding genes are usually far more conserved between F and M than control-region/non-coding sequence — expect the cross-lineage identity to vary considerably across a genome, not to sit at one uniform value.

## Requirements

NOVOPlasty ≥ 4.3.5 (distributed under its own custom, non-commercial license — see [Licensing note](#licensing-note-on-dependencies) below), GetOrganelle ≥ 1.7.7 (GPLv3), BLAST+ ≥ 2.10, minimap2, SAMtools, Python 3, and (optionally, for visualization) Bandage — all on `PATH`.

## Usage (single-lineage mode)

```bash
# Paired-end (default)
./repeatmito.sh -m paired \
  -1 reads_R1.fq.gz -2 reads_R2.fq.gz \
  -s seed.fasta \
  -o output_dir \
  [-t threads] [-r max_rounds] [-g genome_range] [-a plant|animal]

# Single-end
./repeatmito.sh -m single \
  -U reads.fq.gz \
  -s seed.fasta \
  -o output_dir \
  [-t threads] [-r max_rounds] [-g genome_range] [-a plant|animal]
```

| Flag | Meaning | Default |
|---|---|---|
| `-m` | `paired` or `single` read layout | `paired` |
| `-1` / `-2` | Forward/reverse reads (paired mode) | — |
| `-U` | Combined/unpaired reads (single mode) | — |
| `-s` | Seed FASTA (a single conserved mitochondrial gene, or a longer sequence to extend directly); F-type seed in dual-lineage mode | — |
| `-g` | NOVOPlasty `Genome Range` (e.g. `300000-800000`). Narrowing this once you have a size expectation (e.g. from related species) can help NOVOPlasty stop extending into non-target sequence. | `12000-2000000` |
| `-a` | `plant` → NOVOPlasty `Type=mito_plant`, GetOrganelle `-F embplant_mt`. `animal` → `Type=mito`, `-F animal_mt`. Kept for repeat-mediated assembly-break cases in either kingdom. | `plant` |
| `-o` | Output directory | — |
| `-t` | Threads | 8 |
| `-r` | Max Stage A reseeding rounds (per lineage, in dual-lineage mode) | 5 |
| `-M` | M-type seed FASTA — enables DUI dual-lineage mode (see above); requires `-a animal` | none |
| `-x` | Dual-lineage mode: cross-lineage max identity (%) threshold | 98.0 |
| `-y` | Dual-lineage mode: cross-lineage min coverage fraction threshold | 0.5 |

`redundancy_and_consensus.py`'s raw-contig pre-clustering cross-lineage screen has its own, stricter default threshold (`--raw-cross-lineage-min-coverage-frac` 0.3) and is not currently exposed as a `repeatmito.sh` flag — edit the script directly if you need to change it for a given run.

## Interpreting the output (single-lineage mode)

```
output_dir/
├── novoplasty_runs/round_N/           Raw NOVOPlasty output for each Stage A round
├── consensus/
│   ├── round_N/
│   │   └── novel_segment.fasta         What this round actually contributed
│   └── accumulated_nonredundant.fasta  Everything Stage A recovered, across all rounds
├── getorganelle/                       Stage B GetOrganelle run (contains the final .gfa)
├── bandage/
│   ├── repeat_candidates.tsv           Per-node length/coverage/repeat-flag table
│   └── graph_depth.png                 Bandage depth-coloured graph image (if Bandage is available)
└── report/run.log                       Full pipeline log
```

**Before trusting any "novel segment" as genuinely mitochondrial**, it is worth independently re-checking it — repeatMito's Stage C screen catches anomalous coverage/connectivity, but that is not the same as confirming organelle identity. A simple, effective checklist:

1. **GC content.** Plant mitochondrial genomes generally run ~44–48% GC; plant chloroplast genomes tend to run cooler (~36–38%). A segment far outside the mitochondrial range for confirmed relatives is worth a second look regardless of what it BLASTs to.
2. **BLASTN against the closest relative(s) with a published, annotated mitochondrial genome.** Real coverage (even fragmented, even at moderate identity) across a substantial fraction of the segment is the strongest positive signal that it's genuinely mitochondrial.
3. **BLASTN against an assembled/reference plastome, for plant runs.** A handful of percent coverage at scattered, short, low-identity hits is normal background noise. A single large block (kb-scale) at ≥95% identity is a real MTPT and worth annotating as such, not discarding — but if it's a large *fraction* of the segment's total length, or if several segments are dominated by it, suspect contamination or a chimeric contig rather than genuine ancient transfer. (repeatMito does not run this check automatically — see "Background" above; duOrganelle does this by default.)
4. **Coordinate-level gene check.** If you have candidate ORFs/tRNAs on a segment, look up what falls at the *same coordinates* on a relative's mitochondrial genome (and, for plant runs, on the plastome). A "candidate mitochondrial ORF" that lands exactly on a well-known plastid gene (e.g. an ATP synthase subunit, a ribosomal protein, *rbcL*, *matK*, …) is a red flag for cross-organelle mislabeling, even at plausible-looking amino-acid identity — cross-organelle paralogs (like mitochondrial *atp1* vs. plastid *atpA*) can produce real, non-trivial homology scores while being the wrong gene entirely.

## Limitations

Like any short-read-only approach, repeatMito cannot physically phase every alternative recombination product the way long-read spanning can. The reproducible recovery of the same repeat unit across independent reassembly attempts, together with the coverage-depth diagnostic, provide internally cross-validated evidence for multipartite structure — but long-read sequencing or PCR-based junction validation is recommended as a complementary step wherever feasible, particularly before formal annotation or GenBank deposition.

repeatMito does not screen for chloroplast/plastid contamination; use duOrganelle for plant datasets where that protection matters.

Dual-lineage mode has the same short-read limitation, plus its own: the cross-lineage screen is a statistical safety net (identity/coverage thresholds), not a guarantee, and it can only compare against what the pipeline itself has accumulated so far in this run — always cross-check final F/M assignments against independent, published references for the taxon (see "Before trusting a lineage assignment" above) before formal annotation or deposition.

## License

repeatMito's own code is released under the [MIT License](LICENSE) — see the `LICENSE` file for the full text.

### Licensing note on dependencies

repeatMito's MIT license covers only the code in this repository; it does not change the license terms of the external tools it calls. In particular:

- **NOVOPlasty** is distributed under its own custom license (not MIT/GPL/BSD) that restricts use to non-commercial purposes and does not permit redistributing derivative works. If you plan to use repeatMito commercially, you will need to resolve NOVOPlasty's terms separately with its authors (or substitute a compatible alternative). See [ndierckx/NOVOPlasty](https://github.com/ndierckx/NOVOPlasty) for the current license text.
- **GetOrganelle** is GPLv3.
- **Bandage** is GPLv3.

repeatMito invokes these tools as separate external programs rather than bundling or linking their source code, so this repository's own MIT licensing is unaffected — but the restrictions above still apply to how *you* use those dependencies.

## Citing repeatMito

If repeatMito is useful in your work, please cite it. A machine-readable citation is provided in [`CITATION.cff`](CITATION.cff) (GitHub renders a "Cite this repository" button from this file automatically), for example:

> Angulo, F. E., & Núñez, J. J. (2026). *repeatMito* (Version 1.0.0) [Computer software]. https://github.com/LESS-UACh/repeatMito

If you archive a tagged release on Zenodo, replace the GitHub URL above with the release-specific Zenodo DOI so the citation points to an immutable version rather than a moving branch.
