#!/usr/bin/env python3
"""
Inspect FASTQ, BAM, and VCF files with Python.

A minimal, dependency-light tour of genomics file formats using:
- Biopython for FASTQ
- pysam for BAM/SAM
- cyvcf2 for VCF/BCF

Outputs are JSON summaries and CSV/TSV tables written to `outputs/`.
"""

import json
import csv
from pathlib import Path
from collections import Counter

from Bio import SeqIO
import pysam
import cyvcf2

# CIGAR operation codes -> human-readable names
CIGAR_OP_NAMES = {
    pysam.CMATCH: "MATCH",
    pysam.CINS: "INS",
    pysam.CDEL: "DEL",
    pysam.CREF_SKIP: "REF_SKIP",
    pysam.CSOFT_CLIP: "SOFT_CLIP",
    pysam.CHARD_CLIP: "HARD_CLIP",
    pysam.CPAD: "PAD",
    pysam.CEQUAL: "EQUAL",
    pysam.CDIFF: "DIFF",
    pysam.CBACK: "BACK",
}

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "outputs"
OUT.mkdir(parents=True, exist_ok=True)


def write_json(name, data):
    path = OUT / f"{name}.json"
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"  -> {path}")


def write_csv(name, headers, rows):
    path = OUT / f"{name}.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(headers)
        w.writerows(rows)
    print(f"  -> {path}")


# ---------------------------------------------------------------------------
# FASTQ
# ---------------------------------------------------------------------------
def inspect_fastq():
    print("\n=== FASTQ ===")
    r1 = DATA / "NA12878_R1_subset.fastq"
    r2 = DATA / "NA12878_R2_subset.fastq"

    summary = {
        "files": [str(r1.name), str(r2.name)],
        "total_reads": 0,
        "length_stats": {"min": float("inf"), "max": 0, "sum": 0, "counts": []},
        "quality_stats": {"min": float("inf"), "max": 0, "sum": 0, "n": 0},
        "per_file": {},
    }

    per_file_rows = []

    for fastq in [r1, r2]:
        counts = 0
        lens = []
        quals = []
        for record in SeqIO.parse(str(fastq), "fastq"):
            counts += 1
            seq_len = len(record.seq)
            lens.append(seq_len)
            mean_q = sum(record.letter_annotations["phred_quality"]) / seq_len
            quals.append(mean_q)

        summary["total_reads"] += counts
        summary["length_stats"]["min"] = min(summary["length_stats"]["min"], min(lens))
        summary["length_stats"]["max"] = max(summary["length_stats"]["max"], max(lens))
        summary["length_stats"]["sum"] += sum(lens)
        summary["length_stats"]["counts"].extend(lens)

        summary["quality_stats"]["min"] = min(summary["quality_stats"]["min"], min(quals))
        summary["quality_stats"]["max"] = max(summary["quality_stats"]["max"], max(quals))
        summary["quality_stats"]["sum"] += sum(quals)
        summary["quality_stats"]["n"] += len(quals)

        summary["per_file"][fastq.name] = {
            "reads": counts,
            "length_min": min(lens),
            "length_max": max(lens),
            "length_mean": round(sum(lens) / len(lens), 2),
            "quality_mean": round(sum(quals) / len(quals), 2),
        }

        per_file_rows.append([fastq.name, counts, min(lens), max(lens), round(sum(lens) / len(lens), 2), round(sum(quals) / len(quals), 2)])

    total_n = summary["length_stats"]["counts"]
    summary["length_stats"]["mean"] = round(summary["length_stats"]["sum"] / len(total_n), 2)
    summary["quality_stats"]["mean"] = round(summary["quality_stats"]["sum"] / summary["quality_stats"]["n"], 2)

    # Build length distribution (rounded to nearest 5 bp)
    length_dist = Counter((l // 5) * 5 for l in total_n)
    summary["length_distribution"] = {k: v for k, v in sorted(length_dist.items())}

    write_json("fastq_summary", summary)
    write_csv("fastq_per_file", ["file", "reads", "length_min", "length_max", "length_mean", "quality_mean"], per_file_rows)

    dist_rows = [["length_bin", "count"]]
    for k, v in sorted(length_dist.items()):
        dist_rows.append([f"{k}-{k+4}", v])
    write_csv("fastq_length_distribution", ["length_bin", "count"], dist_rows[1:])

    print(f"Total reads: {summary['total_reads']}")
    print(f"Read length: {summary['length_stats']['min']}-{summary['length_stats']['max']} bp (mean {summary['length_stats']['mean']})")
    print(f"Mean Phred quality: {summary['quality_stats']['mean']:.2f}")


# ---------------------------------------------------------------------------
# BAM
# ---------------------------------------------------------------------------
def inspect_bam():
    print("\n=== BAM ===")
    bam_path = DATA / "sample.bam"

    total = 0
    mapped = 0
    unmapped = 0
    paired = 0
    proper_pair = 0
    duplicate = 0
    qcfail = 0
    lengths = []
    mapqs = []
    refs = Counter()
    cigar_ops = Counter()
    insert_sizes = []

    # Iterate reads
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            total += 1
            if not read.is_unmapped:
                mapped += 1
                refs[read.reference_name] += 1
                mapqs.append(read.mapping_quality)
                lengths.append(read.query_length)
                if read.cigar:
                    for op, length in read.cigar:
                        cigar_ops[op] += length
            else:
                unmapped += 1

            if read.is_paired:
                paired += 1
            if read.is_proper_pair:
                proper_pair += 1
            if read.is_duplicate:
                duplicate += 1
            if read.is_qcfail:
                qcfail += 1

            if read.is_proper_pair and read.template_length > 0:
                insert_sizes.append(read.template_length)

    summary = {
        "file": bam_path.name,
        "total_reads": total,
        "mapped": mapped,
        "unmapped": unmapped,
        "paired": paired,
        "proper_pair": proper_pair,
        "duplicate": duplicate,
        "qcfail": qcfail,
        "mapping_quality_mean": round(sum(mapqs) / len(mapqs), 2) if mapqs else 0,
        "mapping_quality_distribution": dict(Counter(mapqs).most_common(20)),
        "read_length_mean": round(sum(lengths) / len(lengths), 2) if lengths else 0,
        "reference_distribution": dict(refs),
        "cigar_operations": {CIGAR_OP_NAMES.get(op, f"op_{op}"): count for op, count in cigar_ops.items()},
        "insert_size_mean": round(sum(insert_sizes) / len(insert_sizes), 2) if insert_sizes else 0,
        "insert_size_distribution": dict(Counter((x // 10) * 10 for x in insert_sizes)),
    }

    write_json("bam_summary", summary)

    ref_rows = [[ref, count] for ref, count in refs.most_common()]
    write_csv("bam_reference_distribution", ["reference", "count"], ref_rows)

    mapq_rows = [[q, c] for q, c in sorted(Counter(mapqs).items())]
    write_csv("bam_mapping_quality", ["mapq", "count"], mapq_rows)

    print(f"Total reads: {total}")
    print(f"Mapped: {mapped} | Unmapped: {unmapped}")
    print(f"References: {dict(refs)}")
    print(f"Mean MAPQ: {summary['mapping_quality_mean']}")


# ---------------------------------------------------------------------------
# VCF
# ---------------------------------------------------------------------------
def inspect_vcf():
    print("\n=== VCF ===")

    for vcf_name in ["basic.vcf", "basic_multisample.vcf", "indexed_tbi.vcf.gz"]:
        vcf_path = DATA / vcf_name
        print(f"\n  {vcf_name}")

        variants = 0
        snps = 0
        indels = 0
        filters = Counter()
        qualities = []
        allele_counts = Counter()
        sample_genotypes = Counter()
        chrom_counts = Counter()

        vcf = cyvcf2.VCF(str(vcf_path))
        samples = vcf.samples

        for variant in vcf:
            variants += 1
            chrom_counts[variant.CHROM] += 1
            filters[variant.FILTER] += 1
            qualities.append(variant.QUAL)

            if variant.is_snp:
                snps += 1
            elif variant.is_indel:
                indels += 1

            ac = variant.INFO.get("AC")
            if ac is not None:
                if isinstance(ac, (list, tuple)):
                    for a in ac:
                        allele_counts[str(a)] += 1
                else:
                    allele_counts[str(ac)] += 1

            if variant.genotypes:
                for gt in variant.genotypes:
                    # genotype is [allele1, allele2, phased]
                    g = f"{gt[0]}/{gt[1]}"
                    sample_genotypes[g] += 1

        summary = {
            "file": vcf_name,
            "samples": samples,
            "num_variants": variants,
            "snps": snps,
            "indels": indels,
            "filters": dict(filters),
            "chromosomes": dict(chrom_counts),
            "quality_mean": round(sum(qualities) / len(qualities), 2) if qualities else 0,
            "quality_distribution": dict(Counter((int(q) // 5) * 5 for q in qualities if q is not None)),
            "allele_counts": dict(allele_counts),
            "genotype_counts": dict(sample_genotypes),
        }

        safe_name = vcf_name.replace(".", "_")
        write_json(f"vcf_summary_{safe_name}", summary)

        filter_rows = [[f, c] for f, c in filters.items()]
        write_csv(f"vcf_filters_{safe_name}", ["filter", "count"], filter_rows)

        chrom_rows = [[c, n] for c, n in chrom_counts.items()]
        write_csv(f"vcf_chromosomes_{safe_name}", ["chrom", "count"], chrom_rows)

        print(f"Variants: {variants} | SNPs: {snps} | Indels: {indels}")
        print(f"Samples: {len(samples)}")
        print(f"Filters: {dict(filters)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("Inspecting genomics formats in", DATA)
    inspect_fastq()
    inspect_bam()
    inspect_vcf()
    print("\nDone. Outputs in", OUT)


if __name__ == "__main__":
    main()
