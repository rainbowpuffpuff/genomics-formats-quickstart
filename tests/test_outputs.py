"""Basic smoke tests for the inspection outputs."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"


def test_fastq_summary():
    path = OUT / "fastq_summary.json"
    assert path.exists()
    data = json.loads(path.read_text())
    assert data["total_reads"] == 200
    assert data["length_stats"]["mean"] == 36.0


def test_bam_summary():
    path = OUT / "bam_summary.json"
    assert path.exists()
    data = json.loads(path.read_text())
    assert data["total_reads"] == 79
    assert data["mapped"] == 79
    assert data["reference_distribution"]["11"] == 79


def test_vcf_summary():
    path = OUT / "vcf_summary_basic_vcf.json"
    assert path.exists()
    data = json.loads(path.read_text())
    assert data["num_variants"] == 48
    assert data["snps"] == 43
    assert data["indels"] == 5
