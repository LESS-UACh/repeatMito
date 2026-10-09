#!/usr/bin/env python3
"""Pre-flight check: does a seed sequence share exact k-mers with the reads?

NOVOPlasty can only start extending a seed if (part of) it matches reads
EXACTLY over its k-mer length (repeatMito's config uses K-mer = 39). A seed
that is too divergent from the sequenced species (e.g. a congeneric
reference for a different genus, or a DUI seed of the wrong lineage) yields
zero extension and a silent "0 bp recovered" run. This script samples the
first N reads and reports how many contain at least one k-mer shared with
the seed (either strand), so that situation is flagged BEFORE the (slow)
assembly starts.

Usage:
    check_seed.py --seed seed.fasta --reads R1.fq.gz [R2.fq.gz ...] \
        [--k 39] [--max-reads 500000] [--label F]

Exit status is always 0 (advisory only). The result line is machine-readable:
    SEED_CHECK label=F k=39 reads_sampled=500000 reads_with_hit=123 status=OK|LOW|NONE
"""
import argparse
import gzip
import sys

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def read_fasta_seq(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--reads", nargs="+", required=True)
    ap.add_argument("--k", type=int, default=39)
    ap.add_argument("--max-reads", type=int, default=500000,
                    help="reads to sample from the START of the first file (default 500000)")
    ap.add_argument("--min-hits", type=int, default=10,
                    help="fewer reads than this with a shared k-mer is reported as LOW (default 10)")
    ap.add_argument("--label", default="seed")
    a = ap.parse_args()

    seed = read_fasta_seq(a.seed)
    k = a.k
    if len(seed) < k:
        print(f"SEED_CHECK label={a.label} k={k} reads_sampled=0 reads_with_hit=0 status=NONE")
        print(f"WARNING: seed is shorter ({len(seed)} bp) than the k-mer size ({k}); NOVOPlasty cannot use it.", file=sys.stderr)
        return 0
    kmers = {seed[i:i + k] for i in range(len(seed) - k + 1)}
    kmers |= {s[::-1].translate(COMP) for s in kmers}

    n = hit = 0
    with opener(a.reads[0]) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            n += 1
            r = line.strip().upper()
            for j in range(len(r) - k + 1):
                if r[j:j + k] in kmers:
                    hit += 1
                    break
            if n >= a.max_reads:
                break

    status = "NONE" if hit == 0 else ("LOW" if hit < a.min_hits else "OK")
    print(f"SEED_CHECK label={a.label} k={k} reads_sampled={n} reads_with_hit={hit} status={status}")
    if status != "OK":
        print(f"WARNING [{a.label}]: only {hit} of {n} sampled reads share an exact {k}-mer with the seed ({a.seed}).", file=sys.stderr)
        print("  NOVOPlasty will probably fail to extend this seed and the run may recover 0 bp.", file=sys.stderr)
        print("  Use a more conserved region (e.g. cox1) or a seed from a closer relative; for DUI make sure", file=sys.stderr)
        print("  the F- and M-type seeds each belong to the right lineage. If the dataset is huge, mitochondrial", file=sys.stderr)
        print("  reads may be rare in the sampled head of the file: re-run with a larger --max-reads.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
