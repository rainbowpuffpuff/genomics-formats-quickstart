#!/usr/bin/env python3
"""Generate an interactive static dashboard for project A in dashboard/index.html."""

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"
DASH = ROOT / "dashboard"
DASH.mkdir(parents=True, exist_ok=True)


def load_json(name):
    with open(OUT / name) as f:
        return json.load(f)


def load_csv(name):
    path = OUT / name
    if path.exists():
        return pd.read_csv(path).to_dict(orient="records")
    return []


def main():
    fastq = load_json("fastq_summary.json")
    bam = load_json("bam_summary.json")
    vcf_basic = load_json("vcf_summary_basic_vcf.json")
    vcf_multi = load_json("vcf_summary_basic_multisample_vcf.json")
    vcf_indexed = load_json("vcf_summary_indexed_tbi_vcf_gz.json")

    fastq_per_file = load_csv("fastq_per_file.csv")
    fastq_len = load_csv("fastq_length_distribution.csv")
    bam_mq = load_csv("bam_mapping_quality.csv")
    bam_ref = load_csv("bam_reference_distribution.csv")
    vcf_chrom = load_csv("vcf_chromosomes_basic_vcf.csv")
    vcf_filters = load_csv("vcf_filters_basic_vcf.csv")

    data = {
        "title": "Genomics Formats Quick Tour",
        "fastq": {
            "total_reads": fastq["total_reads"],
            "length_mean": fastq["length_stats"]["mean"],
            "quality_mean": fastq["quality_stats"]["mean"],
            "per_file": fastq_per_file,
            "length_distribution": fastq_len,
        },
        "bam": {
            "total_reads": bam["total_reads"],
            "mapped": bam["mapped"],
            "unmapped": bam["unmapped"],
            "mapping_quality_mean": bam["mapping_quality_mean"],
            "mapping_quality": bam_mq,
            "reference": bam_ref,
            "insert_size": bam["insert_size_distribution"],
        },
        "vcf": {
            "basic": vcf_basic,
            "multisample": vcf_multi,
            "indexed": vcf_indexed,
            "chromosomes": vcf_chrom,
            "filters": vcf_filters,
        },
    }

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{data['title']}</title>
<script src="https://cdn.tailwindcss.com"></script>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  body {{background:#0f172a; color:#e2e8f0;}}
  .card {{background:#1e293b; border-radius:1rem; padding:1.5rem;}}
</style>
</head>
<body class="font-sans">
<div class="max-w-7xl mx-auto p-6">
  <h1 class="text-3xl font-bold mb-2">{data['title']}</h1>
  <p class="text-slate-400 mb-6">FASTQ, BAM and VCF summaries from tiny sample files</p>

  <h2 class="text-xl font-semibold mb-4">FASTQ</h2>
  <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
    <div class="card"><div class="text-sm text-slate-400">Total reads</div><div class="text-xl font-semibold">{data['fastq']['total_reads']}</div></div>
    <div class="card"><div class="text-sm text-slate-400">Mean length</div><div class="text-xl font-semibold">{data['fastq']['length_mean']}</div></div>
    <div class="card"><div class="text-sm text-slate-400">Mean quality</div><div class="text-xl font-semibold">{data['fastq']['quality_mean']:.2f}</div></div>
    <div class="card"><div class="text-sm text-slate-400">Files</div><div class="text-xl font-semibold">{len(data['fastq']['per_file'])}</div></div>
  </div>
  <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
    <div class="card"><div id="fastq-files" class="h-72"></div></div>
    <div class="card"><div id="fastq-len" class="h-72"></div></div>
  </div>

  <h2 class="text-xl font-semibold mb-4">BAM</h2>
  <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
    <div class="card"><div class="text-sm text-slate-400">Total reads</div><div class="text-xl font-semibold">{data['bam']['total_reads']}</div></div>
    <div class="card"><div class="text-sm text-slate-400">Mapped</div><div class="text-xl font-semibold">{data['bam']['mapped']}</div></div>
    <div class="card"><div class="text-sm text-slate-400">Unmapped</div><div class="text-xl font-semibold">{data['bam']['unmapped']}</div></div>
    <div class="card"><div class="text-sm text-slate-400">Mean MAPQ</div><div class="text-xl font-semibold">{data['bam']['mapping_quality_mean']:.2f}</div></div>
  </div>
  <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
    <div class="card"><div id="bam-mq" class="h-72"></div></div>
    <div class="card"><div id="bam-insert" class="h-72"></div></div>
  </div>

  <h2 class="text-xl font-semibold mb-4">VCF</h2>
  <div class="grid grid-cols-1 md:grid-cols-3 gap-6 mb-6">
    <div class="card"><div id="vcf-basic" class="h-72"></div></div>
    <div class="card"><div id="vcf-multi" class="h-72"></div></div>
    <div class="card"><div id="vcf-indexed" class="h-72"></div></div>
  </div>
</div>

<script>
const data = {json.dumps(data, indent=2)};

const layout = (title) => ({{
  title: {{text: title, font: {{color: '#e2e8f0'}}}},
  paper_bgcolor: 'rgba(0,0,0,0)',
  plot_bgcolor: 'rgba(0,0,0,0)',
  font: {{color: '#94a3b8'}},
  margin: {{t: 40, r: 40, b: 50, l: 50}}
}});

Plotly.newPlot('fastq-files', [{{
  type: 'bar',
  x: data.fastq.per_file.map(p => p.file),
  y: data.fastq.per_file.map(p => p.reads),
  marker: {{color: '#38bdf8'}}
}}], layout('FASTQ reads per file'));

Plotly.newPlot('fastq-len', [{{
  type: 'bar',
  x: data.fastq.length_distribution.map(p => p.length_bin),
  y: data.fastq.length_distribution.map(p => p.count),
  marker: {{color: '#34d399'}}
}}], layout('FASTQ length distribution'));

Plotly.newPlot('bam-mq', [{{
  type: 'bar',
  x: data.bam.mapping_quality.map(p => p.mapq),
  y: data.bam.mapping_quality.map(p => p.count),
  marker: {{color: '#f472b6'}}
}}], layout('BAM mapping quality'));

const insertSizes = Object.entries(data.bam.insert_size).sort((a,b) => parseInt(a[0]) - parseInt(b[0]));
Plotly.newPlot('bam-insert', [{{
  type: 'bar',
  x: insertSizes.map(([k,v]) => k),
  y: insertSizes.map(([k,v]) => v),
  marker: {{color: '#a78bfa'}}
}}], layout('BAM insert size distribution'));

function vcfPie(summary, elemId, title) {{
  const labels = ['SNPs', 'Indels'];
  const values = [summary.snps, summary.indels];
  Plotly.newPlot(elemId, [{{
    type: 'pie',
    labels: labels,
    values: values,
    hole: 0.4,
    marker: {{colors: ['#38bdf8', '#f87171']}},
    textinfo: 'label+percent'
  }}], layout(title + ' — ' + summary.num_variants + ' variants'));
}}

vcfPie(data.vcf.basic, 'vcf-basic', 'VCF basic');
vcfPie(data.vcf.multisample, 'vcf-multi', 'VCF multisample');
vcfPie(data.vcf.indexed, 'vcf-indexed', 'VCF indexed');
</script>
</body>
</html>
"""
    (DASH / "index.html").write_text(html)
    print("Dashboard written to", DASH / "index.html")


if __name__ == "__main__":
    main()
