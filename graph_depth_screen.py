#!/usr/bin/env python3
"""
graph_depth_screen.py

Parses a GetOrganelle/SPAdes assembly graph (.gfa), computes per-node
length and read-depth (coverage), and flags candidate repeat units:
nodes whose coverage substantially exceeds the graph median AND/OR that
are disconnected from the main (largest) connected component.

Flagged nodes are cross-checked by BLASTn against the accumulated
non-redundant segment set to separate previously known sequence from
newly resolved content.

Output: repeat_candidates.tsv, a per-node summary table.

DUI / dual-lineage mode
------------------------
When --cross-lineage-reference is given (the sibling lineage's currently
accumulated segment set, e.g. the M-type track's accumulated fasta when this
invocation is processing the F-type track, or vice versa, in a Doubly
Uniparental Inheritance dual-seed run), every graph node is ALSO screened
against it, as a final safety net in addition to the Stage A NOVOPlasty-level
and redundancy_and_consensus.py-level cross-lineage checks. Strong, DIVERGED
homology to the sibling lineage is EXPECTED (F and M mitochondrial genomes
are homologous but genuinely divergent), so only NEAR-IDENTITY to the
sibling lineage is flagged, as that is the signature of index-hopping, a
NUMT, or read cross-recruitment between the two haplotypes during assembly --
not real independent divergence.
Cross-stage repeat corroboration
---------------------------------
The per-node depth/connectivity screen above only looks at the FINAL,
already-simplified GetOrganelle graph. GetOrganelle sometimes cannot
disentangle a genuinely repeat-mediated locus into a clean circular path and
falls back to forcing a single linear scaffold through it (logged as
"Incomplete/Complicated graph" / "Disentangling unsuccessful") -- which
erases the elevated-depth signature this script otherwise relies on, so the
per-node table can show nothing anomalous even though a real repeat is
present (confirmed independently in validation against a published external
Potamopyrgus antipodarum dataset, where the true locus was instead visible
in Stage A's own all-vs-all self-alignment and in GetOrganelle's log text).

When --self-alignment-tsv (one or more Stage A redundancy_and_consensus.py
all_vs_all.tsv files) and/or --getorganelle-log are given, this script also
checks those two independent, already-computed sources of evidence and
reports them as pipeline-level warnings (not tied to a single graph node,
since the two stages don't share a coordinate system/contig naming scheme)
so a repeat-mediated locus that Stage C's own node table misses is still
surfaced instead of silently passing as "no repeat found".
"""
import argparse
import glob
import os
import re
import subprocess
import statistics
import sys
from collections import defaultdict


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


def find_internal_repeats(tsv_paths, min_len=50, min_pident=90.0):
    """Scan one or more Stage A all_vs_all.tsv (self-BLAST) files for
    non-trivial internal repeats: hits between a sequence and itself that
    are NOT simply the trivial full-length self-match (qstart==sstart and
    qend==send, i.e. length==qlen==slen). Returns a list of dicts describing
    each such hit (source file, identity, length, coordinates, orientation).
    """
    findings = []
    for tsv in tsv_paths:
        if not tsv or not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
            continue
        with open(tsv) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 14:
                    continue
                qseqid, sseqid, pident, length, mismatch, gapopen, \
                    qstart, qend, sstart, send, evalue, bitscore, qlen, slen = parts
                if qseqid != sseqid:
                    continue  # only interested in a sequence vs. itself
                pident = float(pident)
                length = int(length)
                qstart, qend, sstart, send = int(qstart), int(qend), int(sstart), int(send)
                qlen = int(qlen)
                is_trivial_full_length = (
                    qstart == sstart and qend == send and length >= qlen - 1
                )
                if is_trivial_full_length:
                    continue
                if length < min_len or pident < min_pident:
                    continue
                orientation = "inverted" if sstart > send else "direct"
                findings.append({
                    "source": tsv, "seqid": qseqid, "pident": pident, "length": length,
                    "qstart": qstart, "qend": qend, "sstart": sstart, "send": send,
                    "orientation": orientation,
                })
    return findings


def scan_getorganelle_log(log_path, patterns=None):
    """Grep a GetOrganelle log for known phrases indicating it could not
    cleanly disentangle a circular path (a common symptom of a genuine
    repeat-mediated locus that a forced-linear fallback then masks from the
    downstream per-node depth screen). Returns the matching lines, deduped,
    in the order they first appear."""
    if patterns is None:
        patterns = [
            r"Incomplete/Complicated graph",
            r"Disentangling unsuccessful",
            r"may not be as accurate",
        ]
    if not log_path or not os.path.exists(log_path):
        return []
    compiled = [re.compile(p, re.IGNORECASE) for p in patterns]
    seen = []
    with open(log_path, errors="replace") as f:
        for line in f:
            for c in compiled:
                if c.search(line) and line.strip() not in seen:
                    seen.append(line.strip())
    return seen


def parse_gfa(gfa_path):
    nodes = {}
    edges = []
    with open(gfa_path) as f:
        for line in f:
            if line.startswith("S\t"):
                parts = line.rstrip("\n").split("\t")
                name, seq = parts[1], parts[2]
                tags = {}
                for t in parts[3:]:
                    k, typ, v = t.split(":", 2)
                    tags[k] = v
                ln = int(tags.get("LN", len(seq)))
                rc = int(tags.get("RC", 0))
                cov = rc / ln if ln > 0 else 0.0
                nodes[name] = {"seq": seq, "len": ln, "cov": cov}
            elif line.startswith("L\t"):
                parts = line.rstrip("\n").split("\t")
                edges.append((parts[1], parts[3]))
    return nodes, edges


def connected_components(nodes, edges):
    adj = defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    seen = set()
    comps = []
    for n in nodes:
        if n in seen:
            continue
        stack = [n]
        comp = set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack.extend(adj[x] - comp)
        seen |= comp
        comps.append(comp)
    return comps


def merge_intervals(ivs):
    ivs = sorted(ivs)
    merged = []
    for lo, hi in ivs:
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--accumulated", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--depth-fold-threshold", type=float, default=2.0,
                     help="Flag nodes with coverage >= this many times the graph median")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--cross-lineage-reference",
                     help="DUI dual-lineage mode: FASTA of the SIBLING lineage's currently accumulated "
                          "segment set (e.g. the M-type track's accumulated_nonredundant.fasta, when this "
                          "call is processing the F-type track, or vice versa). When given, every graph "
                          "node is also screened against it: strong DIVERGED homology is expected and "
                          "never flagged, but NEAR-IDENTITY (>= --cross-lineage-max-identity) covering a "
                          "substantial fraction of the node (>= --cross-lineage-min-coverage-frac) is "
                          "flagged as suspect cross-lineage contamination -- index-hopping, a nuclear "
                          "mitochondrial pseudogene (NUMT), or read cross-recruitment between the two "
                          "haplotypes -- rather than genuine independent divergence. A full per-node "
                          "identity/coverage report is always written, flagged or not, so observed F/M "
                          "divergence can be checked against what's expected for the taxon.")
    ap.add_argument("--cross-lineage-max-identity", type=float, default=98.0,
                     help="Coverage-weighted percent identity to the sibling lineage's accumulated set, "
                          "at or above which a node is flagged as suspect cross-lineage contamination "
                          "rather than genuine DUI divergence (default 98.0).")
    ap.add_argument("--cross-lineage-min-coverage-frac", type=float, default=0.5,
                     help="Fraction (0-1) of a node's length that must be covered by near-identical "
                          "(>= --cross-lineage-max-identity) cross-lineage hits to flag it (default 0.5).")
    ap.add_argument("--cross-lineage-word-size", type=int, default=11,
                     help="BLASTn word size used specifically for the cross-lineage screen (default 11, "
                          "vs. word size 20 used elsewhere in this script for near-identical matching). "
                          "A smaller word size is needed here because the point of this screen is to "
                          "detect substantially DIVERGED homology between the two mitochondrial lineages, "
                          "not just near-identical matches.")
    ap.add_argument("--self-alignment-tsv", nargs="*", default=[],
                     help="Zero or more Stage A redundancy_and_consensus.py all_vs_all.tsv files (one per "
                          "round; typically found under consensus/round_*/all_vs_all.tsv). Screened for "
                          "non-trivial internal (self-vs-self) repeats that Stage A's own redundancy "
                          "clustering already computed, as a cross-check independent of this script's "
                          "own per-node graph-depth signal.")
    ap.add_argument("--getorganelle-log",
                     help="Path to this round's GetOrganelle log (round_logs/getorganelle.log). Scanned "
                          "for known 'could not disentangle a clean circular path' phrases, since that "
                          "situation can mask a genuine repeat's depth signature once GetOrganelle forces "
                          "a single linear fallback path through it.")
    ap.add_argument("--internal-repeat-min-len", type=int, default=50)
    ap.add_argument("--internal-repeat-min-pident", type=float, default=90.0)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    nodes, edges = parse_gfa(args.gfa)
    if not nodes:
        print("No nodes found in GFA.")
        return

    covs = [v["cov"] for v in nodes.values() if v["cov"] > 0]
    median_cov = statistics.median(covs) if covs else 0.0

    comps = connected_components(nodes, edges)
    comps_sorted = sorted(comps, key=lambda c: -sum(nodes[n]["len"] for n in c))
    main_comp = comps_sorted[0] if comps_sorted else set()

    node_fasta = os.path.join(args.outdir, "all_graph_nodes.fasta")
    with open(node_fasta, "w") as f:
        for name, v in nodes.items():
            f.write(f">{name}\n{v['seq']}\n")

    known_cov = {}
    if os.path.exists(args.accumulated) and os.path.getsize(args.accumulated) > 0:
        tsv = os.path.join(args.outdir, "nodes_vs_accumulated.tsv")
        sh(f"makeblastdb -in {args.accumulated} -dbtype nucl -out {args.accumulated}.db 2>/dev/null")
        sh(
            f"blastn -query {node_fasta} -db {args.accumulated}.db "
            f'-outfmt "6 qseqid sseqid pident length" '
            f"-evalue 1e-20 -word_size 20 -num_threads {args.threads} -out {tsv}"
        )
        best_cov = defaultdict(int)
        with open(tsv) as f:
            for line in f:
                q, s, pident, length = line.strip().split("\t")
                if float(pident) < 90:
                    continue
                best_cov[q] = max(best_cov[q], int(length))
        for name in nodes:
            known_cov[name] = 100.0 * best_cov.get(name, 0) / nodes[name]["len"] if nodes[name]["len"] else 0.0

    # DUI dual-lineage safety net: does this node's length mostly correspond to
    # the SIBLING lineage's currently accumulated set at NEAR-IDENTITY (not just
    # ordinary diverged homology, which is the expected, unflagged signal of two
    # genuinely independent but homologous F/M mitochondrial lineages)? An
    # anomalous depth signal together with near-identity to the sibling lineage
    # is more likely index-hopping / a NUMT / cross-recruited reads than a
    # genuine repeat unit of this lineage. Mirrors the cross-lineage screen in
    # redundancy_and_consensus.py, applied here at the graph-node level as the
    # third and final of the three cross-lineage screening layers (raw contig,
    # cluster consensus, graph node) described in the README.
    sibling_report = {}
    if args.cross_lineage_reference:
        if os.path.exists(args.cross_lineage_reference) and os.path.getsize(args.cross_lineage_reference) > 0:
            sib_tsv = os.path.join(args.outdir, "nodes_vs_sibling_lineage.tsv")
            sh(f"makeblastdb -in {args.cross_lineage_reference} -dbtype nucl -out {args.cross_lineage_reference}.db 2>/dev/null")
            sh(
                f"blastn -query {node_fasta} -db {args.cross_lineage_reference}.db "
                f'-outfmt "6 qseqid sseqid pident length qstart qend" '
                f"-evalue 1e-5 -word_size {args.cross_lineage_word_size} -num_threads {args.threads} -out {sib_tsv}"
            )
            hits_by_q = defaultdict(list)
            with open(sib_tsv) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 6:
                        continue
                    q, s, pident, length, qstart, qend = parts
                    length = int(length)
                    if length < 30:
                        continue
                    lo, hi = sorted([int(qstart), int(qend)])
                    hits_by_q[q].append((float(pident), length, lo, hi))

            for name, v in nodes.items():
                hlist = hits_by_q.get(name, [])
                merged = merge_intervals([(lo, hi) for _, _, lo, hi in hlist])
                covered_bp = sum(hi - lo + 1 for lo, hi in merged)
                pct_covered = 100.0 * covered_bp / v["len"] if v["len"] else 0.0
                total_len = sum(length for _, length, _, _ in hlist)
                weighted_identity = (
                    sum(pident * length for pident, length, _, _ in hlist) / total_len if total_len else 0.0
                )
                sibling_report[name] = (weighted_identity, pct_covered)
        else:
            # Expected on round 1 of whichever lineage runs first: the sibling
            # hasn't accumulated anything yet. Nothing to screen against this round.
            print(f"NOTE: --cross-lineage-reference given ({args.cross_lineage_reference}) but it does not "
                  f"exist yet or is empty (expected on the first round of whichever lineage runs first, "
                  f"before the sibling has accumulated anything) -- skipping cross-lineage screen.",
                  file=sys.stderr)

    out_tsv = os.path.join(args.outdir, "repeat_candidates.tsv")
    with open(out_tsv, "w") as out:
        out.write("node\tlength_bp\tcoverage\tfold_over_median\tconnected_to_main_graph\t"
                   "pct_matches_known_content\tflag_repeat_candidate\t"
                   "pident_vs_sibling_lineage\t"
                   "pct_matches_sibling_lineage\tflag_suspect_cross_lineage_contamination\n")
        for name, v in nodes.items():
            fold = v["cov"] / median_cov if median_cov > 0 else 0.0
            connected = name in main_comp
            pct_known = known_cov.get(name, None)
            flag = (fold >= args.depth_fold_threshold) or (not connected and v["len"] > 1000)

            sib_pident, sib_pct = sibling_report.get(name, (None, None))
            flag_cross_lineage = (
                sib_pident is not None and sib_pct is not None
                and sib_pident >= args.cross_lineage_max_identity
                and sib_pct >= args.cross_lineage_min_coverage_frac * 100
            )

            out.write(
                f"{name}\t{v['len']}\t{v['cov']:.2f}\t{fold:.2f}\t{connected}\t"
                f"{'' if pct_known is None else f'{pct_known:.1f}'}\t"
                f"{flag}\t"
                f"{'' if sib_pident is None else f'{sib_pident:.1f}'}\t"
                f"{'' if sib_pct is None else f'{sib_pct:.1f}'}\t{flag_cross_lineage}\n"
            )
            if flag_cross_lineage:
                print(f"WARNING: node {name} ({v['len']} bp, coverage {v['cov']:.2f}, {fold:.2f}x median) "
                      f"is {sib_pct:.1f}% covered by near-identical ({sib_pident:.1f}% identity) matches to "
                      f"the sibling lineage's accumulated set -- likely suspect cross-lineage contamination "
                      f"(index-hopping, a NUMT, or reads cross-recruited between the two haplotypes), not "
                      f"genuine independent F/M divergence. Recommend reviewing before treating this node "
                      f"as confirmed for this lineage.")

    # Cross-stage corroboration: independent evidence of a repeat-mediated locus
    # from Stage A's own self-alignment output and/or GetOrganelle's log, checked
    # even when (as validated against a real external dataset) the per-node table
    # above shows nothing anomalous because GetOrganelle's forced-linear fallback
    # erased the depth signature.
    internal_repeats = find_internal_repeats(
        args.self_alignment_tsv,
        min_len=args.internal_repeat_min_len,
        min_pident=args.internal_repeat_min_pident,
    )
    disentangle_warnings = scan_getorganelle_log(args.getorganelle_log)

    with open(out_tsv, "a") as out:
        out.write("#\n")
        out.write("# --- Cross-stage repeat corroboration (pipeline-level; not tied to a single node "
                   "above -- see script docstring) ---\n")
        if internal_repeats:
            out.write(f"# STAGE_A_SELF_ALIGNMENT_REPEAT_DETECTED\ttrue\t{len(internal_repeats)} hit(s)\n")
            for h in internal_repeats:
                out.write(
                    f"#   {os.path.basename(h['source'])}: {h['seqid']} {h['length']} bp, "
                    f"{h['pident']:.1f}% identity, {h['orientation']}, "
                    f"query {h['qstart']}-{h['qend']} vs subject {h['sstart']}-{h['send']}\n"
                )
        else:
            out.write("# STAGE_A_SELF_ALIGNMENT_REPEAT_DETECTED\tfalse\n")
        if disentangle_warnings:
            out.write(f"# GETORGANELLE_DISENTANGLE_WARNING\ttrue\t{len(disentangle_warnings)} line(s)\n")
            for w in disentangle_warnings:
                out.write(f"#   {w}\n")
        else:
            out.write("# GETORGANELLE_DISENTANGLE_WARNING\tfalse\n")

    print(f"Graph median coverage: {median_cov:.2f}")
    print(f"Repeat candidate table written to: {out_tsv}")

    if internal_repeats or disentangle_warnings:
        print("WARNING: cross-stage repeat corroboration found evidence of a repeat-mediated locus that "
              "the per-node graph-depth table above may NOT have flagged (this can happen when "
              "GetOrganelle forces a single linear path through a locus it could not fully disentangle, "
              "which erases the depth signature this script otherwise relies on). See the "
              "'Cross-stage repeat corroboration' block appended to repeat_candidates.tsv, and consider "
              "visually inspecting the assembly graph (Bandage) and/or Stage A's consensus contig around "
              "the reported coordinates before concluding this genome is repeat-free.")
        if internal_repeats:
            for h in internal_repeats:
                print(f"  - Stage A self-alignment: {h['seqid']} has a {h['length']} bp internal repeat "
                      f"({h['pident']:.1f}% identity, {h['orientation']}) at "
                      f"{h['qstart']}-{h['qend']} vs {h['sstart']}-{h['send']}")
        if disentangle_warnings:
            for w in disentangle_warnings:
                print(f"  - GetOrganelle log: {w}")


if __name__ == "__main__":
    main()
