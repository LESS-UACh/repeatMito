#!/usr/bin/env python3
"""
redundancy_and_consensus.py

Given a fresh set of NOVOPlasty contigs and the segment set accumulated from
previous rounds, this script:
  1. Runs all-vs-all BLASTn to detect redundant/overlapping contigs.
  2. Reorients redundant contigs to a common strand and builds a majority-rule
     consensus (minimap2 + pileup) for the largest connected cluster.
  3. Masks all content already present in the accumulated segment set and
     extracts genuinely novel, non-redundant sequence.

This formalizes the manual workflow used to characterize repeat-mediated
assembly breaks in animal mitogenomes (e.g. the inverted repeats in the
control region of Potamopyrgus antipodarum).

DUI / dual-lineage mode
------------------------
When --cross-lineage-reference is given (the sibling lineage's currently
accumulated segment set, e.g. the M-type track's accumulated fasta when this
invocation is processing the F-type track, or vice versa, in a Doubly
Uniparental Inheritance dual-seed run), an additional screen runs: strong,
DIVERGED homology to the sibling lineage is EXPECTED (F and M mitochondrial
genomes are homologous but genuinely divergent), so only NEAR-IDENTITY to the
sibling lineage is flagged, as that is the signature of index-hopping, a
NUMT, or read cross-recruitment between the two haplotypes during assembly --
not real independent divergence.
"""
import argparse
import os
import subprocess
import sys
from collections import defaultdict, Counter


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


def load_fasta(path):
    seqs = {}
    name = None
    seq = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return seqs
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            seqs[name] = "".join(seq)
    return seqs


def write_fasta(seqs, path):
    with open(path, "w") as f:
        for k, v in seqs.items():
            f.write(f">{k}\n")
            for i in range(0, len(v), 70):
                f.write(v[i:i + 70] + "\n")


def revcomp(s):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(comp)[::-1]


def run_blast(query, db_fasta, out_tsv, identity, minlen, threads):
    sh(f"makeblastdb -in {db_fasta} -dbtype nucl -out {db_fasta}.db 2>/dev/null")
    sh(
        f"blastn -query {query} -db {db_fasta}.db "
        f'-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" '
        f"-evalue 1e-20 -word_size 20 -num_threads {threads} -out {out_tsv}"
    )


def run_blast_sensitive(query, db_fasta, out_tsv, word_size, threads):
    """BLASTn tuned to detect DIVERGED homology (smaller word size, looser evalue) rather than
    near-identical matches. Used only for the cross-lineage (DUI F-vs-M) divergence screen, where
    genuine, substantially diverged homology is exactly the signal of interest -- not noise to
    exclude the way it is for the near-identical redundancy screen above."""
    sh(f"makeblastdb -in {db_fasta} -dbtype nucl -out {db_fasta}.db 2>/dev/null")
    sh(
        f"blastn -query {query} -db {db_fasta}.db "
        f'-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" '
        f"-evalue 1e-5 -word_size {word_size} -num_threads {threads} -out {out_tsv}"
    )


def parse_hits(tsv, min_identity, min_len):
    hits = []
    if not os.path.exists(tsv):
        return hits
    with open(tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 13:
                continue
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen = parts[:13]
            if qseqid == sseqid:
                continue
            if float(pident) < min_identity or int(length) < min_len:
                continue
            hits.append((qseqid, sseqid, int(length), int(sstart), int(send)))
    return hits


def find_redundant_clusters(seqs, hits, coverage_frac=0.15):
    lens = {k: len(v) for k, v in seqs.items()}
    best = defaultdict(int)
    for q, s, length, sstart, send in hits:
        key = tuple(sorted([q, s]))
        best[key] = max(best[key], length)

    adj = defaultdict(set)
    for (a, b), length in best.items():
        frac = length / min(lens[a], lens[b])
        if frac > coverage_frac:
            adj[a].add(b)
            adj[b].add(a)

    seen = set()
    clusters = []
    for node in seqs:
        if node in seen:
            continue
        stack = [node]
        comp = set()
        while stack:
            n = stack.pop()
            if n in comp:
                continue
            comp.add(n)
            stack.extend(adj[n] - comp)
        seen |= comp
        clusters.append(comp)
    return clusters


def build_consensus(seqs, cluster, hits, outdir, threads):
    """Reorient cluster members to the longest member's strand and build a
    majority-rule consensus via minimap2 + pileup."""
    members = sorted(cluster, key=lambda k: -len(seqs[k]))
    ref_name = members[0]
    ref_seq = seqs[ref_name]

    orient = {ref_name: "+"}
    for m in members[1:]:
        best_hit = None
        for q, s, length, sstart, send in hits:
            if q == m and s == ref_name:
                if best_hit is None or length > best_hit[0]:
                    best_hit = (length, sstart, send)
        if best_hit:
            _, sstart, send = best_hit
            orient[m] = "+" if sstart < send else "-"
        else:
            orient[m] = "+"  # fallback, no direct hit to reference

    reoriented = {}
    for m in members:
        s = seqs[m]
        if orient.get(m) == "-":
            s = revcomp(s)
        reoriented[m] = s

    ref_path = os.path.join(outdir, "cluster_reference.fasta")
    reads_path = os.path.join(outdir, "cluster_reads.fasta")
    write_fasta({ref_name: reoriented[ref_name]}, ref_path)
    write_fasta({m: reoriented[m] for m in members[1:]}, reads_path)

    if len(members) == 1:
        return reoriented[ref_name]

    sam = os.path.join(outdir, "aln.sam")
    bam = os.path.join(outdir, "aln.sorted.bam")
    sh(f"minimap2 -ax asm5 --secondary=no -t {threads} {ref_path} {reads_path} > {sam} 2>/dev/null")
    sh(f"samtools sort -o {bam} {sam} 2>/dev/null")
    sh(f"samtools index {bam} 2>/dev/null")
    pileup = os.path.join(outdir, "pileup.txt")
    sh(f"samtools mpileup -f {ref_path} -a {bam} > {pileup} 2>/dev/null")

    def parse_pileup_bases(bases, ref_base):
        calls = []
        i, n = 0, len(bases)
        while i < n:
            c = bases[i]
            if c == "^":
                i += 2
                continue
            if c == "$":
                i += 1
                continue
            if c in "+-":
                i += 1
                num = ""
                while i < n and bases[i].isdigit():
                    num += bases[i]
                    i += 1
                i += int(num)
                continue
            if c in ".,":
                calls.append(ref_base.upper())
                i += 1
                continue
            if c in "ACGTNacgtn":
                calls.append(c.upper())
                i += 1
                continue
            if c == "*":
                i += 1
                continue
            i += 1
        return calls

    consensus = list(reoriented[ref_name])
    with open(pileup) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            _, pos, ref_base, depth, bases, _ = parts[:6]
            pos = int(pos)
            calls = parse_pileup_bases(bases, ref_base)
            votes = calls + [ref_base.upper()]
            best_base, _ = Counter(votes).most_common(1)[0]
            if best_base in "ACGT":
                consensus[pos - 1] = best_base
    return "".join(consensus)


def merge_intervals(ivs):
    ivs = sorted(ivs)
    merged = []
    for lo, hi in ivs:
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def screen_cross_lineage(seqs, reference_fasta, outdir, max_identity, min_coverage_frac,
                          word_size, threads, tag):
    """DUI dual-lineage screen: compare `seqs` against the SIBLING lineage's currently
    accumulated segment set (e.g. F contigs vs. the M-type accumulated set, or vice versa).

    Strong, DIVERGED
    homology to the sibling lineage is the EXPECTED signal (F and M mitochondrial genomes
    are homologous but often substantially diverged, even in conserved protein-coding
    genes), so it is never itself a reason to exclude anything. Only NEAR-IDENTITY,
    covering a substantial fraction of the sequence, is flagged -- that combination is the
    signature of index-hopping, a nuclear mitochondrial pseudogene (NUMT), or a chimeric
    contig produced by cross-recruitment of reads between the two haplotypes, not genuine
    independent divergence.

    Returns (clean_seqs, flagged_seqs, report: dict[name] -> (weighted_pident, pct_covered)).
    Every input sequence appears in `report`, whether or not it was flagged, so the observed
    F/M divergence can be sanity-checked against what is expected for the taxon.
    """
    if not reference_fasta or not seqs:
        return seqs, {}, {}
    if not os.path.exists(reference_fasta) or os.path.getsize(reference_fasta) == 0:
        # Expected on round 1 of whichever lineage runs first: the sibling hasn't
        # accumulated anything yet. Nothing to screen against this round.
        return seqs, {}, {}

    tmp_fasta = os.path.join(outdir, f"_cross_lineage_query_{tag}.fasta")
    write_fasta(seqs, tmp_fasta)
    tsv = os.path.join(outdir, f"vs_sibling_lineage_{tag}.tsv")
    run_blast_sensitive(tmp_fasta, reference_fasta, tsv, word_size, threads)

    hits_by_q = defaultdict(list)
    with open(tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 13:
                continue
            q, s, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen = parts[:13]
            length = int(length)
            if length < 30:
                continue
            lo, hi = sorted([int(qstart), int(qend)])
            hits_by_q[q].append((float(pident), length, lo, hi))

    clean, flagged, report = {}, {}, {}
    for name, seq in seqs.items():
        hlist = hits_by_q.get(name, [])
        merged = merge_intervals([(lo, hi) for _, _, lo, hi in hlist])
        covered_bp = sum(hi - lo + 1 for lo, hi in merged)
        pct_covered = 100.0 * covered_bp / len(seq) if len(seq) else 0.0
        total_len = sum(length for _, length, _, _ in hlist)
        weighted_identity = (
            sum(pident * length for pident, length, _, _ in hlist) / total_len if total_len else 0.0
        )
        report[name] = (weighted_identity, pct_covered)
        if weighted_identity >= max_identity and pct_covered >= min_coverage_frac * 100:
            flagged[name] = seq
        else:
            clean[name] = seq
    return clean, flagged, report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--accumulated", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--identity", type=float, default=90)
    ap.add_argument("--minlen", type=int, default=200)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--cross-lineage-reference",
                     help="DUI dual-lineage mode: FASTA of the SIBLING lineage's currently accumulated "
                          "segment set (e.g. the M-type track's accumulated_nonredundant.fasta, when this "
                          "call is processing the F-type track, or vice versa). When given, every raw "
                          "contig and every cluster consensus is screened against it; strong DIVERGED "
                          "homology is expected and never flagged, but NEAR-IDENTITY (>= "
                          "--cross-lineage-max-identity) covering a substantial fraction of the sequence "
                          "(>= --cross-lineage-min-coverage-frac) is flagged as suspect cross-lineage "
                          "contamination -- index-hopping, a NUMT, or read cross-recruitment between the "
                          "two haplotypes -- and excluded from this lineage's accumulated set. A full "
                          "per-sequence identity/coverage report is always written, flagged or not, so "
                          "observed F/M divergence can be checked against what's expected for the taxon.")
    ap.add_argument("--cross-lineage-max-identity", type=float, default=98.0,
                     help="Coverage-weighted percent identity to the sibling lineage's accumulated set, "
                          "at or above which a sequence is flagged as suspect cross-lineage contamination "
                          "rather than genuine DUI divergence (default 98.0). Genuine F/M divergence, even "
                          "in conserved protein-coding genes, is normally well below this in DUI bivalves; "
                          "near-identity is the signature of contamination or NUMT capture, not biology.")
    ap.add_argument("--cross-lineage-min-coverage-frac", type=float, default=0.5,
                     help="Fraction (0-1) of a sequence's length that must be covered by near-identical "
                          "(>= --cross-lineage-max-identity) cross-lineage hits to flag it (default 0.5). "
                          "A short near-identical patch (e.g. a real recombination breakpoint or a small "
                          "NUMT) covering only a small fraction of a much longer contig will not by itself "
                          "trigger the flag; it still appears in the per-sequence report for manual review.")
    ap.add_argument("--raw-cross-lineage-min-coverage-frac", type=float, default=0.3,
                     help="Same idea as --cross-lineage-min-coverage-frac, applied to RAW contigs before "
                          "clustering (default 0.3, stricter): single-linkage clustering can transitively "
                          "pull a contaminated/chimeric raw contig into a cluster with genuinely clean "
                          "contigs before the post-consensus screen ever gets a chance to see it.")
    ap.add_argument("--cross-lineage-word-size", type=int, default=11,
                     help="BLASTn word size used specifically for the cross-lineage screen (default 11, "
                          "vs. word size 20 used elsewhere in this script for near-identical matching). "
                          "A smaller word size is needed here because the point of this screen is to "
                          "detect substantially DIVERGED homology between the two mitochondrial lineages, "
                          "not just near-identical matches.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    seqs = load_fasta(args.contigs)
    if not seqs:
        write_fasta({}, os.path.join(args.outdir, "novel_segment.fasta"))
        return

    # 0) DUI dual-lineage pre-screen on RAW contigs: exclude anything suspiciously
    # near-identical to the sibling lineage's accumulated set BEFORE clustering, since
    # single-linkage clustering could otherwise transitively merge a contaminated raw
    # contig into a cluster with genuinely clean ones via shared homology.
    if args.cross_lineage_reference:
        seqs, flagged_cross_raw, cross_raw_report = screen_cross_lineage(
            seqs, args.cross_lineage_reference, args.outdir,
            args.cross_lineage_max_identity, args.raw_cross_lineage_min_coverage_frac,
            args.cross_lineage_word_size, args.threads, "raw",
        )
        if cross_raw_report:
            cross_raw_report_path = os.path.join(args.outdir, "cross_lineage_screen_report_RAW_contigs.tsv")
            with open(cross_raw_report_path, "w") as rf:
                rf.write("raw_contig\tcoverage_weighted_pident_vs_sibling_lineage\tpct_length_covered\tverdict\n")
                for name, (pid, pct) in sorted(cross_raw_report.items(), key=lambda x: -x[1][1]):
                    verdict = ("FLAGGED_SUSPECT_CROSS_LINEAGE_excluded_before_clustering"
                               if name in flagged_cross_raw else "ok")
                    rf.write(f"{name}\t{pid:.1f}\t{pct:.1f}\t{verdict}\n")
            if flagged_cross_raw:
                write_fasta(flagged_cross_raw, os.path.join(args.outdir, "flagged_cross_lineage_raw_contigs.fasta"))
                print(f"WARNING: {len(flagged_cross_raw)} RAW contig(s) excluded BEFORE clustering as "
                      f"suspect cross-lineage contamination (>= {args.cross_lineage_max_identity:.0f}% "
                      f"identity and >= {args.raw_cross_lineage_min_coverage_frac*100:.0f}% coverage vs. "
                      f"the sibling lineage's accumulated set): {', '.join(flagged_cross_raw)}. "
                      f"See {cross_raw_report_path}.")
        if not seqs:
            write_fasta({}, os.path.join(args.outdir, "novel_segment.fasta"))
            print("All raw contigs this round were flagged as suspect cross-lineage contamination; "
                  "nothing left to cluster.")
            return

    # 1) all-vs-all redundancy within this round's contigs (post cross-lineage pre-screen)
    clean_contigs_fasta = os.path.join(args.outdir, "_raw_contigs_post_screens.fasta")
    write_fasta(seqs, clean_contigs_fasta)
    aa_tsv = os.path.join(args.outdir, "all_vs_all.tsv")
    run_blast(clean_contigs_fasta, clean_contigs_fasta, aa_tsv, args.identity, args.minlen, args.threads)
    hits = parse_hits(aa_tsv, args.identity, args.minlen)
    clusters = find_redundant_clusters(seqs, hits)

    consensus_seqs = {}
    for i, cluster in enumerate(clusters):
        cdir = os.path.join(args.outdir, f"cluster_{i}")
        os.makedirs(cdir, exist_ok=True)
        if len(cluster) == 1:
            name = next(iter(cluster))
            consensus_seqs[f"cluster{i}_{name}"] = seqs[name]
        else:
            cons = build_consensus(seqs, cluster, hits, cdir, args.threads)
            consensus_seqs[f"cluster{i}_consensus_n{len(cluster)}"] = cons

    consolidated_round = os.path.join(args.outdir, "round_consolidated.fasta")
    write_fasta(consensus_seqs, consolidated_round)

    # 1b) DUI dual-lineage safeguard: screen every cluster consensus against the sibling
    # lineage's accumulated set. Ordinary strong homology here is expected (F and M are
    # the same organelle type, divergent) -- only near-identity is excluded, as it
    # signals contamination/NUMT rather than biology.
    if args.cross_lineage_reference:
        consensus_seqs, flagged_cross, cross_report = screen_cross_lineage(
            consensus_seqs, args.cross_lineage_reference, args.outdir,
            args.cross_lineage_max_identity, args.cross_lineage_min_coverage_frac,
            args.cross_lineage_word_size, args.threads, "consensus",
        )
        if cross_report:
            cross_report_path = os.path.join(args.outdir, "cross_lineage_divergence_report.tsv")
            with open(cross_report_path, "w") as rf:
                rf.write("sequence\tcoverage_weighted_pident_vs_sibling_lineage\tpct_length_covered\tverdict\n")
                for name, (pid, pct) in sorted(cross_report.items(), key=lambda x: -x[1][1]):
                    verdict = ("FLAGGED_SUSPECT_CROSS_LINEAGE" if name in flagged_cross
                               else "ok_consistent_with_independent_divergence")
                    rf.write(f"{name}\t{pid:.1f}\t{pct:.1f}\t{verdict}\n")
            if flagged_cross:
                write_fasta(flagged_cross, os.path.join(args.outdir, "flagged_cross_lineage_contamination.fasta"))
                print(f"WARNING: {len(flagged_cross)} cluster consensus sequence(s) flagged as suspect "
                      f"cross-lineage contamination (>= {args.cross_lineage_max_identity:.0f}% identity, "
                      f">= {args.cross_lineage_min_coverage_frac*100:.0f}% coverage, vs. the sibling "
                      f"lineage's accumulated set) and EXCLUDED from this lineage's accumulated set. "
                      f"See {cross_report_path} and flagged_cross_lineage_contamination.fasta.")
        elif args.cross_lineage_reference and (not os.path.exists(args.cross_lineage_reference)
                                                or os.path.getsize(args.cross_lineage_reference) == 0):
            print(f"NOTE: --cross-lineage-reference given ({args.cross_lineage_reference}) but it does not "
                  f"exist yet or is empty (expected on the first round of whichever lineage runs first, "
                  f"before the sibling has accumulated anything) -- skipping cross-lineage screen this round.",
                  file=sys.stderr)

    # 2) mask content already present in the accumulated segment set -> novel segment(s)
    if os.path.exists(args.accumulated) and os.path.getsize(args.accumulated) > 0:
        vs_acc_tsv = os.path.join(args.outdir, "vs_accumulated.tsv")
        run_blast(consolidated_round, args.accumulated, vs_acc_tsv, 90, 100, args.threads)
        known_hits = parse_hits(vs_acc_tsv, 90, 100)
        covered = defaultdict(list)
        for q, s, length, sstart, send in known_hits:
            pass
        # need query-side coordinates, not subject; rerun parse with qstart/qend
        covered_q = defaultdict(list)
        with open(vs_acc_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 13:
                    continue
                q, s, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen = parts[:13]
                if float(pident) < 90 or int(length) < 100:
                    continue
                lo, hi = sorted([int(qstart), int(qend)])
                covered_q[q].append((lo, hi))

        def merge(ivs):
            ivs = sorted(ivs)
            merged = []
            for lo, hi in ivs:
                if merged and lo <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
                else:
                    merged.append((lo, hi))
            return merged

        novel_parts = {}
        for name, seq in consensus_seqs.items():
            L = len(seq)
            merged = merge(covered_q.get(name, []))
            pos = 1
            gaps = []
            for lo, hi in merged:
                if lo > pos:
                    gaps.append((pos, lo - 1))
                pos = max(pos, hi + 1)
            if pos <= L:
                gaps.append((pos, L))
            for gi, (lo, hi) in enumerate(gaps):
                if hi - lo + 1 < 200:
                    continue
                novel_parts[f"{name}_novel_{lo}_{hi}"] = seq[lo - 1:hi]
    else:
        novel_parts = consensus_seqs

    write_fasta(novel_parts, os.path.join(args.outdir, "novel_segment.fasta"))
    print(f"Round complete: {len(clusters)} redundancy cluster(s), "
          f"{sum(len(v) for v in novel_parts.values())} bp novel sequence extracted.")


if __name__ == "__main__":
    main()
