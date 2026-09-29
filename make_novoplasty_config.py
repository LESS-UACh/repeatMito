#!/usr/bin/env python3
"""Generate a NOVOPlasty config.txt for a given seed, supporting both
paired-end and single-end read inputs.

Used internally by repeatmito.sh; can also be called standalone.

Paired-end usage:
    make_novoplasty_config.py --seed seed.fasta --mode paired \
        --r1 R1.fq.gz --r2 R2.fq.gz --out config.txt --project round1

Single-end usage:
    make_novoplasty_config.py --seed seed.fasta --mode single \
        --reads reads.fq.gz --out config.txt --project round1
"""
import argparse
import os

TEMPLATE_PE = """Project:
-----------------------
Project name          = {project}
Type                   = {novo_type}
Genome Range           = {genome_range}
K-mer                  = 39
Max memory             = 20
Extended log           = 0
Save assembled reads   = no
Seed Input             = {seed}
Extend seed directly   = {extend_directly}
Reference sequence     =
Variance detection     =
Chloroplast sequence   =

Dataset 1:
-----------------------
Read Length            = 150
Insert size             = 300
Platform                = illumina
Single/Paired           = PE
Combined reads          =
Forward reads           = {r1}
Reverse reads           = {r2}
Store Hash              =

Optional:
-----------------------
Insert size auto        = yes
Use Quality Scores      = no
"""

TEMPLATE_SE = """Project:
-----------------------
Project name          = {project}
Type                   = {novo_type}
Genome Range           = {genome_range}
K-mer                  = 39
Max memory             = 20
Extended log           = 0
Save assembled reads   = no
Seed Input             = {seed}
Extend seed directly   = {extend_directly}
Reference sequence     =
Variance detection     =
Chloroplast sequence   =

Dataset 1:
-----------------------
Read Length            = 150
Insert size             =
Platform                = illumina
Single/Paired           = SE
Combined reads          = {reads}
Forward reads           =
Reverse reads           =
Store Hash              =

Optional:
-----------------------
Insert size auto        = yes
Use Quality Scores      = no
"""


def seed_length(seed_path):
    seed_len = 0
    seq = []
    with open(seed_path) as f:
        for line in f:
            if line.startswith(">"):
                if seq:
                    seed_len = max(seed_len, len("".join(seq)))
                seq = []
            else:
                seq.append(line.strip())
    if seq:
        seed_len = max(seed_len, len("".join(seq)))
    return seed_len


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", required=True)
    ap.add_argument("--mode", choices=["paired", "single"], default="paired",
                     help="Read layout: 'paired' (default) requires --r1/--r2; "
                          "'single' requires --reads")
    ap.add_argument("--r1", help="Forward reads (paired mode)")
    ap.add_argument("--r2", help="Reverse reads (paired mode)")
    ap.add_argument("--reads", help="Combined/unpaired reads file (single mode)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--project", default="assembly")
    ap.add_argument("--organism-type", choices=["plant", "animal"], default="plant",
                     help="Maps to NOVOPlasty's 'Type' field: 'plant' -> mito_plant (default), "
                          "'animal' -> mito. repeatMito is scoped to repeat-mediated assembly breaks and "
                          "DUI dual-lineage recovery; it does not perform chloroplast-aware assembly "
                          "(NOVOPlasty's 'Chloroplast sequence' field is always left blank -- for joint "
                          "mitochondrial + chloroplast/plastid genome assembly, use duOrganelle instead).")
    ap.add_argument("--genome-range", default="12000-2000000",
                     help="NOVOPlasty 'Genome Range' field: expected min-max size (bp) of the target "
                          "organelle genome, used to bound the extension. Default 12000-2000000 (NOVOPlasty's "
                          "own generic default, deliberately loose). Narrowing it (e.g. 300000-800000 once "
                          "you have a size expectation from related species / prior assembly attempts) can "
                          "help NOVOPlasty stop extending into non-target sequence.")
    args = ap.parse_args()

    if args.mode == "paired":
        if not args.r1 or not args.r2:
            ap.error("--mode paired requires both --r1 and --r2")
    else:
        if not args.reads:
            ap.error("--mode single requires --reads")

    # Heuristic: if the seed fasta's single sequence is long (>5 kb),
    # NOVOPlasty should extend it directly rather than treat it as a short gene seed.
    seed_len = seed_length(args.seed)
    extend_directly = "yes" if seed_len > 5000 else "no"
    novo_type = "mito_plant" if args.organism_type == "plant" else "mito"

    if args.mode == "paired":
        cfg = TEMPLATE_PE.format(
            project=args.project,
            seed=os.path.abspath(args.seed),
            extend_directly=extend_directly,
            novo_type=novo_type,
            genome_range=args.genome_range,
            r1=os.path.abspath(args.r1),
            r2=os.path.abspath(args.r2),
        )
    else:
        cfg = TEMPLATE_SE.format(
            project=args.project,
            seed=os.path.abspath(args.seed),
            extend_directly=extend_directly,
            novo_type=novo_type,
            genome_range=args.genome_range,
            reads=os.path.abspath(args.reads),
        )

    with open(args.out, "w") as f:
        f.write(cfg)
    print(f"Wrote {args.out} (mode={args.mode}, seed_len={seed_len}bp, "
          f"Extend seed directly={extend_directly}, Type={novo_type})")

    if args.mode == "single":
        print("NOTE: NOVOPlasty's paired-read cross-validation of each seed "
              "extension is unavailable in single-end mode. The Stage A "
              "redundancy-detection and majority-vote consensus step becomes "
              "a more important quality-control safeguard as a result.")


if __name__ == "__main__":
    main()
