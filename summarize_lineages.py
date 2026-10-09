#!/usr/bin/env python3
"""Per-lineage run summary for repeatMito (single- or dual-lineage mode).

Reads the output tree written by repeatmito.sh and writes
<outdir>/report/lineage_summary.tsv, with one row per lineage:

    lineage, stageA_bp, stageA_n_seqs, stageA_circularized,
    stageB_bp, stageB_n_seqs, stageC_nodes, stageC_nodes_flagged_cross_lineage,
    stageC_max_pident_vs_sibling, stageC_repeat_candidates, status, recommendation

`status` / `recommendation` encode what the Perumytilus purpuratus validation
showed (see the manuscript, Section 3.3): in a DUI dataset where one lineage is
much more abundant, GetOrganelle's reassembly (Stage B) of the LESS abundant
lineage can converge on the MORE abundant lineage's genome. Stage C's
cross-lineage screen flags this (flag_suspect_cross_lineage_contamination), and in
that case the Stage A assembly, not the Stage B one, is the lineage-specific result.

Usage: summarize_lineages.py --outdir OUTDIR [--dual]
"""
import argparse
import csv
import glob
import os
import sys


def fasta_stats(paths):
    bp = n = 0
    for p in paths:
        try:
            with open(p) as fh:
                for line in fh:
                    if line.startswith(">"):
                        n += 1
                    else:
                        bp += len(line.strip())
        except OSError:
            pass
    return bp, n


def parse_stage_c(path):
    nodes = flagged = repeats = 0
    flagged_bp = total_bp = 0
    max_pid = None
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        rows = [l for l in fh if l.strip() and not l.startswith("#")]
    if not rows:
        return None
    for r in csv.DictReader(rows, delimiter="\t"):
        nodes += 1
        try:
            ln = int(float(r.get("length_bp") or 0))
        except ValueError:
            ln = 0
        total_bp += ln
        if str(r.get("flag_suspect_cross_lineage_contamination", "")).strip().lower() == "true":
            flagged += 1
            flagged_bp += ln
        if str(r.get("flag_repeat_candidate", "")).strip().lower() == "true":
            repeats += 1
        pid = r.get("pident_vs_sibling_lineage", "")
        try:
            v = float(pid)
            max_pid = v if max_pid is None else max(max_pid, v)
        except ValueError:
            pass
    return dict(nodes=nodes, flagged=flagged, repeats=repeats, flagged_bp=flagged_bp,
                total_bp=total_bp, max_pid=max_pid)


def summarize(root, label):
    acc = os.path.join(root, "consensus", "accumulated_nonredundant.fasta")
    a_bp, a_n = fasta_stats([acc])
    circ = bool(glob.glob(os.path.join(root, "novoplasty_runs", "round_*", "Circularized_assembly*")))
    b_files = sorted(glob.glob(os.path.join(root, "getorganelle", "*.path_sequence.fasta")))
    b_bp, b_n = fasta_stats(b_files)
    c = parse_stage_c(os.path.join(root, "bandage", "repeat_candidates.tsv"))

    if a_bp == 0:
        status, rec = "STAGE_A_EMPTY", ("Stage A recovered no sequence. Check the seed (see the SEED_CHECK line in run.log) "
                                        "and the NOVOPlasty log under round_logs/.")
    elif c is None or b_bp == 0:
        status, rec = "STAGE_B_OR_C_MISSING", "Stage B/C did not produce a graph/report; use the Stage A set and inspect round_logs/."
    elif c["total_bp"] and c["flagged_bp"] >= 0.5 * c["total_bp"]:
        status = "STAGE_B_COLLAPSED_ONTO_SIBLING"
        rec = ("Most of the Stage B graph is near-identical to the sibling lineage: GetOrganelle likely converged on the "
               "more abundant lineage. Treat the Stage A assembly (consensus/accumulated_nonredundant.fasta"
               + ("; circularized" if circ else "") + ") as this lineage's result and verify it manually.")
    elif c["flagged"] > 0:
        status = "PARTIAL_CROSS_LINEAGE_FLAG"
        rec = "Some Stage B nodes are near-identical to the sibling lineage; review them before accepting this lineage."
    else:
        status = "OK"
        rec = "No cross-lineage flag; compare Stage A and Stage B lengths and verify against references."
        if abs(a_bp - b_bp) > 0.01 * max(a_bp, b_bp) and not (a_n > 1):
            rec += f" NOTE: Stage A ({a_bp} bp) and Stage B ({b_bp} bp) differ by >1%."
    return dict(lineage=label, stageA_bp=a_bp, stageA_n_seqs=a_n, stageA_circularized=str(circ).lower(),
                stageB_bp=b_bp, stageB_n_seqs=b_n,
                stageC_nodes=(c or {}).get("nodes", ""),
                stageC_nodes_flagged_cross_lineage=(c or {}).get("flagged", ""),
                stageC_max_pident_vs_sibling=(c or {}).get("max_pid", ""),
                stageC_repeat_candidates=(c or {}).get("repeats", ""),
                status=status, recommendation=rec)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--dual", action="store_true", help="dual-lineage layout (outdir/F and outdir/M)")
    a = ap.parse_args()
    lineages = [("F", os.path.join(a.outdir, "F")), ("M", os.path.join(a.outdir, "M"))] if a.dual \
        else [("single", a.outdir)]
    rows = [summarize(root, label) for label, root in lineages]
    os.makedirs(os.path.join(a.outdir, "report"), exist_ok=True)
    out = os.path.join(a.outdir, "report", "lineage_summary.tsv")
    cols = list(rows[0].keys())
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Lineage summary: {out}")
    for r in rows:
        print(f"  [{r['lineage']}] Stage A {r['stageA_bp']} bp"
              f"{' (circularized)' if r['stageA_circularized']=='true' else ''}; Stage B {r['stageB_bp']} bp; "
              f"status={r['status']}")
        if r["status"] != "OK":
            print(f"      -> {r['recommendation']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
