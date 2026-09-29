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
"""
import argparse
import os
import subprocess
import statistics
import sys
from collections import defaultdict


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


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
    ap.add_argument("--chloroplast",
                     help="FASTA of the assembled/reference chloroplast genome. When given, every graph "
                          "node is also screened against it (final safety net, in addition to the Stage A "
                          "NOVOPlasty-level and redundancy_and_consensus.py-level checks): a node whose "
                          "length is mostly (>= --plastid-frac-threshold) explained by near-identical "
                          "chloroplast matches is flagged, since an anomalous coverage-depth signal on such "
                          "a node may reflect plastid-read cross-mapping rather than a genuine mitochondrial "
                          "repeat unit.")
    ap.add_argument("--plastid-frac-threshold", type=float, default=0.5)
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

    # Final safety-net organelle-identity screen: does this node's length mostly
    # correspond to the chloroplast genome? An anomalous depth signal on a
    # mostly-plastid node is more likely plastid-read cross-mapping than a
    # genuine mitochondrial repeat.
    plastid_cov = {}
    if args.chloroplast and os.path.exists(args.chloroplast) and os.path.getsize(args.chloroplast) > 0:
        cp_tsv = os.path.join(args.outdir, "nodes_vs_chloroplast.tsv")
        sh(f"makeblastdb -in {args.chloroplast} -dbtype nucl -out {args.chloroplast}.db 2>/dev/null")
        sh(
            f"blastn -query {node_fasta} -db {args.chloroplast}.db "
            f'-outfmt "6 qseqid sseqid pident length qstart qend" '
            f"-evalue 1e-20 -word_size 20 -num_threads {args.threads} -out {cp_tsv}"
        )
        covered_q = defaultdict(list)
        with open(cp_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue
                q, s, pident, length, qstart, qend = parts
                if float(pident) < 90:
                    continue
                lo, hi = sorted([int(qstart), int(qend)])
                covered_q[q].append((lo, hi))

        for name, v in nodes.items():
            merged = merge_intervals(covered_q.get(name, []))
            covered_bp = sum(hi - lo + 1 for lo, hi in merged)
            plastid_cov[name] = 100.0 * covered_bp / v["len"] if v["len"] else 0.0

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
                   "pct_matches_known_content\tpct_matches_chloroplast\tflag_repeat_candidate\t"
                   "flag_likely_plastid_contamination\tpident_vs_sibling_lineage\t"
                   "pct_matches_sibling_lineage\tflag_suspect_cross_lineage_contamination\n")
        for name, v in nodes.items():
            fold = v["cov"] / median_cov if median_cov > 0 else 0.0
            connected = name in main_comp
            pct_known = known_cov.get(name, None)
            pct_plastid = plastid_cov.get(name, None)
            flag = (fold >= args.depth_fold_threshold) or (not connected and v["len"] > 1000)
            flag_plastid = (pct_plastid is not None) and (pct_plastid >= args.plastid_frac_threshold * 100)

            sib_pident, sib_pct = sibling_report.get(name, (None, None))
            flag_cross_lineage = (
                sib_pident is not None and sib_pct is not None
                and sib_pident >= args.cross_lineage_max_identity
                and sib_pct >= args.cross_lineage_min_coverage_frac * 100
            )

            out.write(
                f"{name}\t{v['len']}\t{v['cov']:.2f}\t{fold:.2f}\t{connected}\t"
                f"{'' if pct_known is None else f'{pct_known:.1f}'}\t"
                f"{'' if pct_plastid is None else f'{pct_plastid:.1f}'}\t{flag}\t{flag_plastid}\t"
                f"{'' if sib_pident is None else f'{sib_pident:.1f}'}\t"
                f"{'' if sib_pct is None else f'{sib_pct:.1f}'}\t{flag_cross_lineage}\n"
            )
            if flag_plastid:
                print(f"WARNING: node {name} ({v['len']} bp, coverage {v['cov']:.2f}, "
                      f"{fold:.2f}x median) is {pct_plastid:.1f}% explained by chloroplast matches -- "
                      f"likely plastid contamination, not a mitochondrial repeat unit. Recommend excluding "
                      f"it from the final mitochondrial genome/table before annotation.")
            if flag_cross_lineage:
                print(f"WARNING: node {name} ({v['len']} bp, coverage {v['cov']:.2f}, {fold:.2f}x median) "
                      f"is {sib_pct:.1f}% covered by near-identical ({sib_pident:.1f}% identity) matches to "
                      f"the sibling lineage's accumulated set -- likely suspect cross-lineage contamination "
                      f"(index-hopping, a NUMT, or reads cross-recruited between the two haplotypes), not "
                      f"genuine independent F/M divergence. Recommend reviewing before treating this node "
                      f"as confirmed for this lineage.")

    print(f"Graph median coverage: {median_cov:.2f}")
    print(f"Repeat candidate table written to: {out_tsv}")


if __name__ == "__main__":
    main()
