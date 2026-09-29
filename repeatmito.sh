#!/bin/bash
# =============================================================================
# repeatMito: Iterative reseeding and graph-based resolution of repeat-mediated
#             multipartite mitochondrial genome structure from short reads alone
#             -- plus a DUI dual-lineage mode for animal mitochondrial genomes
#             with Doubly Uniparental Inheritance (independent F-type / M-type
#             mitochondrial lineages).
#
# repeatMito is scoped to two problems:
#   1. Repeat-mediated assembly breaks: large repeated sequences that a
#      short-read assembly either collapses or fragments around, producing a
#      multipartite mitochondrial genome (first documented in a plant,
#      Fuchsia magellanica; the same failure mode occurs in animal mitogenomes,
#      e.g. Brachidontes).
#   2. Doubly Uniparental Inheritance (DUI): recovering two divergent,
#      co-occurring F-type/M-type mitochondrial lineages from the same
#      heteroplasmic read pool without letting them cross-contaminate.
# Joint mitochondrial + chloroplast/plastid genome assembly is handled by the
# companion pipeline duOrganelle, not by repeatMito.
#
# Combines NOVOPlasty (iterative seed-extension) with GetOrganelle + Bandage
# (assembly-graph depth screening) to detect and characterize repeat-mediated
# multipartite mitochondrial genome structure using Illumina-only data.
#
# Supports both paired-end and single-end read layouts (-m paired|single).
#
# Usage (paired-end, default, single-lineage -- unchanged from before):
#   ./repeatmito.sh -m paired -1 reads_R1.fq.gz -2 reads_R2.fq.gz \
#                    -s seed.fasta -o output_dir [-t threads] [-r max_rounds]
#
# Usage (single-end, single-lineage):
#   ./repeatmito.sh -m single -U reads.fq.gz \
#                    -s seed.fasta -o output_dir [-t threads] [-r max_rounds]
#
# Usage (DUI dual-lineage mode, animal only -- give BOTH -s and -M):
#   ./repeatmito.sh -m paired -1 reads_R1.fq.gz -2 reads_R2.fq.gz -a animal \
#                    -s F_type_seed.fasta -M M_type_seed.fasta -o output_dir \
#                    [-x cross_lineage_max_identity] [-y cross_lineage_min_coverage_frac]
#
#   -s is the F-type (maternally-transmitted) seed, -M is the M-type
#   (paternally-transmitted) seed. Each lineage gets its own fully independent
#   Stage A/B/C run, under output_dir/F and output_dir/M respectively. The two
#   tracks are NEVER pooled or clustered together -- each track's raw contigs
#   and cluster consensuses are instead screened, every round, against the
#   SIBLING lineage's currently-accumulated set: ordinary diverged homology
#   (expected -- F and M are the same organelle type, divergent) is never
#   flagged, but near-identity is, since that signals index-hopping, a NUMT,
#   or read cross-recruitment between haplotypes rather than real divergence.
#   See redundancy_and_consensus.py / graph_depth_screen.py --cross-lineage-*.
#
# Requires on PATH: NOVOPlasty4.3.5.pl (or newer), blastn, makeblastdb,
#                    minimap2, samtools, get_organelle_from_reads.py, Bandage
# =============================================================================
set -euo pipefail

MODE="paired"
THREADS=8
MAX_ROUNDS=5
MIN_NOVEL_BP=500        # stop iterating a lineage if a round yields less new sequence than this
BLAST_IDENTITY=90
BLAST_MINLEN=200

# Internal placeholder label used for single-lineage mode. bash associative
# arrays reject "" as a subscript ("bad array subscript"), so single-lineage
# mode cannot use an empty label as its array key the way lineage_root's
# directory logic does -- this sentinel stands in for that case everywhere
# an array key or a per-lineage log tag is needed. It is never used as a
# real lineage label (those are always "F"/"M", hardcoded below).
SINGLE_LABEL="_single_"

GENOME_RANGE="12000-2000000"
ORGANISM_TYPE="plant"

SEED_M=""                       # M-type seed: presence of -M enables DUI dual-lineage mode
CROSS_MAX_IDENTITY=98.0
CROSS_MIN_COV=0.5

usage() {
  echo "Usage (paired, single-lineage): $0 -m paired -1 R1.fq.gz -2 R2.fq.gz -s seed.fasta -o outdir [-t threads] [-r max_rounds] [-g genome_range] [-a plant|animal]"
  echo "Usage (single, single-lineage): $0 -m single -U reads.fq.gz -s seed.fasta -o outdir [-t threads] [-r max_rounds] [-g genome_range] [-a plant|animal]"
  echo "Usage (DUI dual-lineage, animal only): $0 -m paired -1 R1.fq.gz -2 R2.fq.gz -a animal -s F_seed.fasta -M M_seed.fasta -o outdir [-x cross_lineage_max_identity] [-y cross_lineage_min_coverage_frac]"
  echo "  -g sets NOVOPlasty's 'Genome Range' (e.g. -g 300000-800000). Default: 12000-2000000."
  echo "  -a sets NOVOPlasty's 'Type': plant -> mito_plant (default), animal -> mito."
  echo "  -M <M_type_seed.fasta> enables DUI dual-lineage mode (requires -a animal). -s is then the"
  echo "     F-type seed. Each lineage is assembled independently under outdir/F and outdir/M."
  echo "  -x cross-lineage max identity (%, default 98.0) and -y cross-lineage min coverage fraction"
  echo "     (0-1, default 0.5): thresholds for flagging suspect F/M cross-contamination in dual-lineage mode."
  echo "  Joint mitochondrial + chloroplast/plastid assembly is out of scope here -- see duOrganelle."
  exit 1
}

while getopts "m:1:2:U:s:o:t:r:g:a:M:x:y:h" opt; do
  case $opt in
    m) MODE="$OPTARG" ;;
    1) R1="$OPTARG" ;;
    2) R2="$OPTARG" ;;
    U) SEREADS="$OPTARG" ;;
    s) SEED="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    r) MAX_ROUNDS="$OPTARG" ;;
    g) GENOME_RANGE="$OPTARG" ;;
    a) ORGANISM_TYPE="$OPTARG" ;;
    M) SEED_M="$OPTARG" ;;
    x) CROSS_MAX_IDENTITY="$OPTARG" ;;
    y) CROSS_MIN_COV="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

if [ "$ORGANISM_TYPE" != "plant" ] && [ "$ORGANISM_TYPE" != "animal" ]; then
  echo "ERROR: -a must be 'plant' or 'animal' (got '$ORGANISM_TYPE')"; usage
fi

DUAL_MODE=false
if [ -n "$SEED_M" ]; then
  DUAL_MODE=true
  if [ "$ORGANISM_TYPE" != "animal" ]; then
    echo "ERROR: -M (DUI dual-lineage mode) requires -a animal."; exit 1
  fi
fi

if [ "$MODE" != "paired" ] && [ "$MODE" != "single" ]; then
  echo "ERROR: -m must be 'paired' or 'single'"; usage
fi

if [ "$MODE" = "paired" ]; then
  { [ -z "${R1:-}" ] || [ -z "${R2:-}" ]; } && { echo "ERROR: paired mode requires -1 and -2"; usage; }
else
  [ -z "${SEREADS:-}" ] && { echo "ERROR: single mode requires -U"; usage; }
fi
[ -z "${SEED:-}" ] || [ -z "${OUTDIR:-}" ] && usage
if $DUAL_MODE && [ ! -s "$SEED_M" ]; then
  echo "ERROR: -M seed fasta '$SEED_M' not found or empty."; exit 1
fi

mkdir -p "$OUTDIR/report"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== repeatMito pipeline started: $(date) (mode=$MODE, organism=$ORGANISM_TYPE, dual_lineage=$DUAL_MODE) ===" | tee "$OUTDIR/report/run.log"
if [ "$MODE" = "single" ]; then
  echo "NOTE: single-end mode disables NOVOPlasty's paired-read cross-validation of seed extensions." | tee -a "$OUTDIR/report/run.log"
  echo "      Stage A redundancy detection + majority-vote consensus is relied upon more heavily as a result." | tee -a "$OUTDIR/report/run.log"
fi
if $DUAL_MODE; then
  echo "=== DUI dual-lineage mode: F-type seed=$SEED, M-type seed=$SEED_M ===" | tee -a "$OUTDIR/report/run.log"
  echo "    Cross-lineage screen thresholds: max identity=$CROSS_MAX_IDENTITY%, min coverage frac=$CROSS_MIN_COV" | tee -a "$OUTDIR/report/run.log"
fi

# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

# lineage_root <label> -> output directory root for that lineage.
# label="" (single-lineage mode) resolves to $OUTDIR itself, preserving the
# original flat directory layout exactly when -M is not used.
lineage_root() {
  local label="$1"
  if [ "$label" = "$SINGLE_LABEL" ]; then echo "$OUTDIR"; else echo "$OUTDIR/$label"; fi
}

# find_novoplasty_output <round_dir>
# NOVOPlasty names its output FASTA differently depending on outcome (see
# https://github.com/ndierckx/NOVOPlasty/wiki/Output-files):
#   Circularized_assembly_*  -- fully closed into one contig (best case)
#   Contigs_*                -- assembled but not circularized (multiple contigs)
#   Merged_contigs_*         -- multiple contigs, merge attempted
#   Option_*                 -- alternative possible circular configurations
#   contigs_tmp_*             -- last-resort fallback, only if nothing else was written
# Searched in that priority order; echoes the first match found (path) and
# returns 0, or returns 1 if none exist.
find_novoplasty_output() {
  local round_dir="$1"
  local f
  for pattern in "Circularized_assembly*" "Contigs_*" "Merged_contigs_*" "Option_*" "contigs_tmp_*"; do
    f=$(find "$round_dir" -maxdepth 1 \( -iname "${pattern}.fasta" -o -iname "${pattern}.txt" \) 2>/dev/null | sort | head -1)
    if [ -n "$f" ]; then
      echo "$f"
      return 0
    fi
  done
  return 1
}

# run_stage_a_round <label> <round> <seed_path> <accumulated_path> <other_accumulated_or_empty>
# Runs one Stage A round (NOVOPlasty -> redundancy_and_consensus.py) for one lineage.
# Returns 0 and updates <root>/.last_novel_segment_<label> if the round should continue
# iterating (novel yield >= MIN_NOVEL_BP); returns 1 if this lineage should stop.
run_stage_a_round() {
  local label="$1" round="$2" seed_path="$3" accumulated_path="$4" other_accum="$5"
  local root; root=$(lineage_root "$label")
  local tag=""; [ "$label" != "$SINGLE_LABEL" ] && tag=" [$label]"

  echo "--- Round $round${tag}: NOVOPlasty with seed $seed_path (mode=$MODE) ---" | tee -a "$OUTDIR/report/run.log"
  local round_dir="$root/novoplasty_runs/round_${round}"
  mkdir -p "$round_dir"

  local proj="round${round}"; [ "$label" != "$SINGLE_LABEL" ] && proj="round${round}_${label}"

  if [ "$MODE" = "paired" ]; then
    python3 "$SCRIPT_DIR/make_novoplasty_config.py" \
      --seed "$seed_path" --mode paired --r1 "$R1" --r2 "$R2" \
      --genome-range "$GENOME_RANGE" --organism-type "$ORGANISM_TYPE" \
      --out "$round_dir/config.txt" --project "$proj"
  else
    python3 "$SCRIPT_DIR/make_novoplasty_config.py" \
      --seed "$seed_path" --mode single --reads "$SEREADS" \
      --genome-range "$GENOME_RANGE" --organism-type "$ORGANISM_TYPE" \
      --out "$round_dir/config.txt" --project "$proj"
  fi

  ( cd "$round_dir" && NOVOPlasty4.3.5.pl -c config.txt ) \
    > "$root/round_logs/novoplasty_round${round}.log" 2>&1 || {
      echo "WARNING: NOVOPlasty round $round${tag} failed or produced no contigs; stopping iteration for this lineage." | tee -a "$OUTDIR/report/run.log"
      return 1
    }

  local contigs
  if ! contigs=$(find_novoplasty_output "$round_dir"); then
    echo "No NOVOPlasty output produced in round $round${tag}; stopping iteration for this lineage." | tee -a "$OUTDIR/report/run.log"
    echo "  NOVOPlasty exited without error but wrote none of Circularized_assembly_*/Contigs_*/" | tee -a "$OUTDIR/report/run.log"
    echo "  Merged_contigs_*/Option_*/contigs_tmp_* -- this usually means the seed never found a match" | tee -a "$OUTDIR/report/run.log"
    echo "  in the reads, or extension failed internally. Last 30 lines of" | tee -a "$OUTDIR/report/run.log"
    echo "  $root/round_logs/novoplasty_round${round}.log:" | tee -a "$OUTDIR/report/run.log"
    tail -n 30 "$root/round_logs/novoplasty_round${round}.log" 2>/dev/null | sed 's/^/    | /' | tee -a "$OUTDIR/report/run.log"
    return 1
  fi

  local circularized=false
  case "$(basename "$contigs")" in
    Circularized_assembly*) circularized=true ;;
  esac
  echo "Round $round${tag}: NOVOPlasty output: $(basename "$contigs")$($circularized && echo " -- GENOME CIRCULARIZED (single closed contig)")" | tee -a "$OUTDIR/report/run.log"

  local cons_dir="$root/consensus/round_${round}"
  python3 "$SCRIPT_DIR/redundancy_and_consensus.py" \
    --contigs "$contigs" \
    --accumulated "$accumulated_path" \
    --outdir "$cons_dir" \
    --identity "$BLAST_IDENTITY" --minlen "$BLAST_MINLEN" \
    --threads "$THREADS" \
    ${other_accum:+--cross-lineage-reference "$other_accum" --cross-lineage-max-identity "$CROSS_MAX_IDENTITY" --cross-lineage-min-coverage-frac "$CROSS_MIN_COV"}

  local new_seg="$cons_dir/novel_segment.fasta"
  local new_bp
  new_bp=$(python3 -c "
seq=''
try:
    with open('$new_seg') as f:
        for l in f:
            if not l.startswith('>'):
                seq+=l.strip()
except FileNotFoundError:
    pass
print(len(seq))
")
  echo "Round $round${tag}: $new_bp bp of novel non-redundant sequence recovered." | tee -a "$OUTDIR/report/run.log"
  cat "$new_seg" >> "$accumulated_path" 2>/dev/null || true
  echo "$new_seg" > "$root/.last_novel_segment_${label}"

  if $circularized; then
    echo "Lineage${tag:- (single)} reports a circularized genome as of round $round; stopping reseeding for" | tee -a "$OUTDIR/report/run.log"
    echo "  this lineage (re-extending an already-closed circular genome is not meaningful/safe)." | tee -a "$OUTDIR/report/run.log"
    return 1
  fi

  if [ "$new_bp" -lt "$MIN_NOVEL_BP" ]; then
    echo "Novel yield below threshold ($MIN_NOVEL_BP bp) for lineage${tag:- (single)}; stopping iteration." | tee -a "$OUTDIR/report/run.log"
    return 1
  fi
  return 0
}

# run_stage_bc <label> <accumulated_path> <other_accumulated_or_empty>
run_stage_bc() {
  local label="$1" accumulated="$2" other_accum="$3"
  local root; root=$(lineage_root "$label")
  local tag=""; [ "$label" != "$SINGLE_LABEL" ] && tag=" [$label]"

  if [ ! -s "$accumulated" ]; then
    echo "SKIPPING Stage B/C${tag}: Stage A recovered no sequence for this lineage (accumulated set is empty)." | tee -a "$OUTDIR/report/run.log"
    echo "  Nothing to seed GetOrganelle with. See the round_1 NOVOPlasty log under $root/round_logs/ for why." | tee -a "$OUTDIR/report/run.log"
    return 1
  fi

  local getorg_f
  if [ "$ORGANISM_TYPE" = "plant" ]; then getorg_f="embplant_mt"; else getorg_f="animal_mt"; fi

  echo "--- Stage B${tag}: GetOrganelle multi-seed assembly (mode=$MODE, -F $getorg_f) ---" | tee -a "$OUTDIR/report/run.log"
  # NOTE: "$root/getorganelle" is deliberately NOT pre-created here.
  # get_organelle_from_reads.py refuses to write into an output directory that
  # already exists (even empty) unless --continue/--overwrite is passed.
  if [ "$MODE" = "paired" ]; then
    get_organelle_from_reads.py \
      -1 "$R1" -2 "$R2" \
      -s "$accumulated" \
      -F "$getorg_f" -R 30 -t "$THREADS" \
      -o "$root/getorganelle" --overwrite \
      > "$root/round_logs/getorganelle.log" 2>&1
  else
    get_organelle_from_reads.py \
      -u "$SEREADS" \
      -s "$accumulated" \
      -F "$getorg_f" -R 30 -t "$THREADS" \
      -o "$root/getorganelle" --overwrite \
      > "$root/round_logs/getorganelle.log" 2>&1
  fi

  if [ ! -d "$root/getorganelle" ]; then
    echo "ERROR: GetOrganelle${tag} did not create an output directory at all (exited without producing" | tee -a "$OUTDIR/report/run.log"
    echo "  $root/getorganelle) -- check $root/round_logs/getorganelle.log for the actual error. Last 30 lines:" | tee -a "$OUTDIR/report/run.log"
    tail -n 30 "$root/round_logs/getorganelle.log" 2>/dev/null | sed 's/^/    | /' | tee -a "$OUTDIR/report/run.log"
    return 1
  fi

  # GetOrganelle's output directory can contain several .gfa files -- one per
  # intermediate k-mer round, plus the final "*.selected_graph.gfa" (the
  # organelle-only assembly graph with correct per-node depth annotation).
  # Prefer that one explicitly; only fall back to "any .gfa" if it's missing.
  local gfa; gfa=$(find "$root/getorganelle" -iname "*.selected_graph.gfa" | sort | head -1)
  if [ -z "$gfa" ]; then
    echo "NOTE: no *.selected_graph.gfa found${tag}; falling back to the first *.gfa found -- verify this" | tee -a "$OUTDIR/report/run.log"
    echo "  is really the final assembly graph, not an intermediate k-mer round." | tee -a "$OUTDIR/report/run.log"
    gfa=$(find "$root/getorganelle" -iname "*.gfa" | sort | head -1)
  fi
  if [ -z "$gfa" ]; then
    echo "ERROR: GetOrganelle produced no .gfa assembly graph${tag}. Last 30 lines of $root/round_logs/getorganelle.log:" | tee -a "$OUTDIR/report/run.log"
    tail -n 30 "$root/round_logs/getorganelle.log" 2>/dev/null | sed 's/^/    | /' | tee -a "$OUTDIR/report/run.log"
    return 1
  fi
  echo "Using assembly graph${tag}: $gfa" | tee -a "$OUTDIR/report/run.log"

  echo "--- Stage C${tag}: Graph node depth screening (repeat/cross-lineage detection) ---" | tee -a "$OUTDIR/report/run.log"
  python3 "$SCRIPT_DIR/graph_depth_screen.py" \
    --gfa "$gfa" \
    --accumulated "$accumulated" \
    --outdir "$root/bandage" \
    --threads "$THREADS" \
    ${other_accum:+--cross-lineage-reference "$other_accum" --cross-lineage-max-identity "$CROSS_MAX_IDENTITY" --cross-lineage-min-coverage-frac "$CROSS_MIN_COV"}

  QT_QPA_PLATFORM=offscreen Bandage image "$gfa" "$root/bandage/graph_depth.png" \
    --colour depth --depth --height 1200 2>/dev/null || \
    echo "NOTE: Bandage image rendering skipped (no offscreen Qt platform available)${tag}." | tee -a "$OUTDIR/report/run.log"

  echo "Summary report${tag}: $root/bandage/repeat_candidates.tsv" | tee -a "$OUTDIR/report/run.log"
}

# -----------------------------------------------------------------------
# Set up lineage tracks
# -----------------------------------------------------------------------
if $DUAL_MODE; then
  LABELS=(F M)
else
  LABELS=("$SINGLE_LABEL")
fi

declare -A ACCUM CUR_SEED ACTIVE
for label in "${LABELS[@]}"; do
  root=$(lineage_root "$label")
  mkdir -p "$root"/{round_logs,novoplasty_runs,consensus,bandage,report}
  ACCUM[$label]="$root/consensus/accumulated_nonredundant.fasta"
  > "${ACCUM[$label]}"
  ACTIVE[$label]=1
done

if $DUAL_MODE; then
  CUR_SEED[F]="$SEED"
  CUR_SEED[M]="$SEED_M"
else
  CUR_SEED["$SINGLE_LABEL"]="$SEED"
fi

# -----------------------------------------------------------------------
# STAGE A: Iterative NOVOPlasty reseeding with redundancy-aware consensus,
# one or two lineage tracks in lock-step by round number. In dual-lineage
# mode, each track's cross-lineage screen sees the sibling's state as of its
# most recently COMPLETED round -- so round 1 of whichever track runs first
# has nothing to cross-screen against yet (nothing accumulated on the other
# side); from round 2 onward both tracks screen against genuinely up-to-date
# sibling data. This is logged explicitly below and in each round's output.
# -----------------------------------------------------------------------
for round in $(seq 1 "$MAX_ROUNDS"); do
  any_active=false
  for label in "${LABELS[@]}"; do
    [ "${ACTIVE[$label]}" = "1" ] || continue
    any_active=true
    other=""
    if $DUAL_MODE; then
      sibling="M"; [ "$label" = "F" ] || sibling="F"
      other="${ACCUM[$sibling]}"
    fi
    if run_stage_a_round "$label" "$round" "${CUR_SEED[$label]}" "${ACCUM[$label]}" "$other"; then
      root=$(lineage_root "$label")
      CUR_SEED[$label]="$(cat "$root/.last_novel_segment_${label}")"
    else
      ACTIVE[$label]=0
    fi
  done
  $any_active || break
done

for label in "${LABELS[@]}"; do
  tag=""; [ "$label" != "$SINGLE_LABEL" ] && tag=" [$label]"
  echo "Accumulated non-redundant segment set${tag}: ${ACCUM[$label]}" | tee -a "$OUTDIR/report/run.log"
done

# -----------------------------------------------------------------------
# STAGE B + C: GetOrganelle multi-seed assembly + graph depth screening,
# run independently per lineage track.
# -----------------------------------------------------------------------
for label in "${LABELS[@]}"; do
  other=""
  if $DUAL_MODE; then
    sibling="M"; [ "$label" = "F" ] || sibling="F"
    other="${ACCUM[$sibling]}"
  fi
  run_stage_bc "$label" "${ACCUM[$label]}" "$other"
done

echo "=== repeatMito pipeline finished: $(date) ===" | tee -a "$OUTDIR/report/run.log"
if $DUAL_MODE; then
  echo "F-type summary report: $OUTDIR/F/bandage/repeat_candidates.tsv"
  echo "M-type summary report: $OUTDIR/M/bandage/repeat_candidates.tsv"
  echo "Check both lineages' cross_lineage_divergence_report.tsv (under consensus/round_N/) and"
  echo "the pident_vs_sibling_lineage / flag_suspect_cross_lineage_contamination columns of each"
  echo "repeat_candidates.tsv before treating any sequence as confirmed F- or M-type mitochondrial content."
else
  echo "Summary report: $OUTDIR/bandage/repeat_candidates.tsv"
fi
