#!/usr/bin/env python3

import csv
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, Set


@dataclass(frozen=True, order=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    variant_type: str
    indel_subtype: str


def load_variant_tsv(tsv_path: str) -> Set[Variant]:
    variants: Set[Variant] = set()

    with open(tsv_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"CHROM", "POS", "REF", "ALT", "TYPE", "INDEL_SUBTYPE"}

        missing = required_columns.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Missing required columns in {tsv_path}: {', '.join(sorted(missing))}"
            )

        for row in reader:
            variants.add(
                Variant(
                    chrom=row["CHROM"],
                    pos=int(row["POS"]),
                    ref=row["REF"],
                    alt=row["ALT"],
                    variant_type=row["TYPE"],
                    indel_subtype=row["INDEL_SUBTYPE"],
                )
            )

    return variants


def write_variant_table(path: str, variants: Iterable[Variant]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["CHROM", "POS", "REF", "ALT", "TYPE", "INDEL_SUBTYPE"])
        for v in sorted(variants):
            writer.writerow([v.chrom, v.pos, v.ref, v.alt, v.variant_type, v.indel_subtype])


def subset_by_indel_type(variants: Set[Variant], subtype: str) -> Set[Variant]:
    return {v for v in variants if v.indel_subtype == subtype}


def compute_summary(
    clara_indels: Set[Variant],
    dragen_indels: Set[Variant],
    shared_indels: Set[Variant],
    clara_only_indels: Set[Variant],
    dragen_only_indels: Set[Variant],
) -> Dict[str, float]:
    union_total = len(clara_indels | dragen_indels)

    shared_pct = (len(shared_indels) / union_total * 100.0) if union_total else 0.0
    clara_only_pct = (len(clara_only_indels) / union_total * 100.0) if union_total else 0.0
    dragen_only_pct = (len(dragen_only_indels) / union_total * 100.0) if union_total else 0.0

    shared_del = subset_by_indel_type(shared_indels, "DEL")
    shared_ins = subset_by_indel_type(shared_indels, "INS")

    clara_only_del = subset_by_indel_type(clara_only_indels, "DEL")
    clara_only_ins = subset_by_indel_type(clara_only_indels, "INS")

    dragen_only_del = subset_by_indel_type(dragen_only_indels, "DEL")
    dragen_only_ins = subset_by_indel_type(dragen_only_indels, "INS")

    return {
        "clara_total": len(clara_indels),
        "dragen_total": len(dragen_indels),
        "shared": len(shared_indels),
        "clara_only": len(clara_only_indels),
        "dragen_only": len(dragen_only_indels),
        "union_total": union_total,
        "shared_pct": shared_pct,
        "clara_only_pct": clara_only_pct,
        "dragen_only_pct": dragen_only_pct,
        "shared_del": len(shared_del),
        "shared_ins": len(shared_ins),
        "clara_only_del": len(clara_only_del),
        "clara_only_ins": len(clara_only_ins),
        "dragen_only_del": len(dragen_only_del),
        "dragen_only_ins": len(dragen_only_ins),
    }


def write_summary(path: str, summary: Dict[str, float]) -> None:
    with open(path, "w", encoding="utf-8") as out:
        out.write("category=INDELs\n")
        out.write(f"clara_total={summary['clara_total']}\n")
        out.write(f"dragen_total={summary['dragen_total']}\n")
        out.write(f"shared={summary['shared']}\n")
        out.write(f"clara_only={summary['clara_only']}\n")
        out.write(f"dragen_only={summary['dragen_only']}\n")
        out.write(f"union_total={summary['union_total']}\n")
        out.write(f"shared_pct={summary['shared_pct']:.2f}\n")
        out.write(f"clara_only_pct={summary['clara_only_pct']:.2f}\n")
        out.write(f"dragen_only_pct={summary['dragen_only_pct']:.2f}\n")
        out.write(f"shared_del={summary['shared_del']}\n")
        out.write(f"shared_ins={summary['shared_ins']}\n")
        out.write(f"clara_only_del={summary['clara_only_del']}\n")
        out.write(f"clara_only_ins={summary['clara_only_ins']}\n")
        out.write(f"dragen_only_del={summary['dragen_only_del']}\n")
        out.write(f"dragen_only_ins={summary['dragen_only_ins']}\n")


def print_summary(summary: Dict[str, float]) -> None:
    print("\n=========== INDEL COMPARISON SUMMARY ===========")
    print(">> Totals by pipeline (DEL + INS):")
    print(f"  CLARA:  {summary['clara_total']}")
    print(f"  DRAGEN: {summary['dragen_total']}")
    print()
    print(">> Shared InDels:")
    print(f"  Shared deletions:   {summary['shared_del']}")
    print(f"  Shared insertions:  {summary['shared_ins']}")
    print(f"  Total shared:       {summary['shared']}")
    print()
    print(">> CLARA-only InDels:")
    print(f"  Exclusive deletions CLARA:  {summary['clara_only_del']}")
    print(f"  Exclusive insertions CLARA: {summary['clara_only_ins']}")
    print(f"  Total exclusive CLARA:      {summary['clara_only']}")
    print()
    print(">> DRAGEN-only InDels:")
    print(f"  Exclusive deletions DRAGEN:  {summary['dragen_only_del']}")
    print(f"  Exclusive insertions DRAGEN: {summary['dragen_only_ins']}")
    print(f"  Total exclusive DRAGEN:      {summary['dragen_only']}")
    print()
    print(">> Percentages over total InDels considered:")
    print(f"  Total InDels considered: {summary['union_total']}")
    print(f"  Shared:             {summary['shared_pct']:.2f}%  ({summary['shared']} variants)")
    print(f"  Exclusive CLARA:    {summary['clara_only_pct']:.2f}%  ({summary['clara_only']} variants)")
    print(f"  Exclusive DRAGEN:   {summary['dragen_only_pct']:.2f}%  ({summary['dragen_only']} variants)")
    print("================================================\n")


def main() -> None:
    if len(sys.argv) != 4:
        print(
            "Usage: python3 src/04_compare_indels.py "
            "<clara/indels.tsv> <dragen/indels.tsv> <output_dir>"
        )
        sys.exit(1)

    clara_tsv = sys.argv[1]
    dragen_tsv = sys.argv[2]
    outdir = sys.argv[3]

    for path in (clara_tsv, dragen_tsv):
        if not os.path.isfile(path):
            print(f"ERROR: File not found: {path}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(outdir, exist_ok=True)

    print("[INFO] Loading CLARA InDels...")
    clara_indels = load_variant_tsv(clara_tsv)

    print("[INFO] Loading DRAGEN InDels...")
    dragen_indels = load_variant_tsv(dragen_tsv)

    shared_indels = clara_indels & dragen_indels
    clara_only_indels = clara_indels - dragen_indels
    dragen_only_indels = dragen_indels - clara_indels

    write_variant_table(os.path.join(outdir, "shared_indels.tsv"), shared_indels)
    write_variant_table(os.path.join(outdir, "clara_only_indels.tsv"), clara_only_indels)
    write_variant_table(os.path.join(outdir, "dragen_only_indels.tsv"), dragen_only_indels)

    summary = compute_summary(
        clara_indels=clara_indels,
        dragen_indels=dragen_indels,
        shared_indels=shared_indels,
        clara_only_indels=clara_only_indels,
        dragen_only_indels=dragen_only_indels,
    )

    write_summary(os.path.join(outdir, "summary.txt"), summary)
    print_summary(summary)

    print(f"[INFO] Results written to: {outdir}")
    print("[INFO] Step 4 completed successfully.")


if __name__ == "__main__":
    main()
