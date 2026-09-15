# CUT&Tag Snakemake Upstream Analysis Pipeline Template

This template allows laboratory members to reproduce the same upstream analysis from paired-end FASTQ files using the **Trim Galore → Bowtie2 → SEACR / deepTools** pipeline. It requires **Linux and Bash**. Run all commands from the project root directory.

The workflow covers raw FASTQ quality control, adapter and quality trimming, post-trimming quality control, target genome alignment, optional spike-in alignment, BAM sorting and indexing, fragment coverage, peak calling, blacklist filtering, bigWig generation, fragment length plots, and signal plots around peaks. Data downloading, differential peak analysis, gene annotation, enrichment analysis, merged group tracks, and other project-specific analyses are outside its scope.

## 1. Assumptions and requirements

- Replace every `[path_to_your_...]` placeholder with an actual absolute path, and replace `sample_01` and `sample_02` with your sample IDs. The placeholders must be replaced before running the workflow.
- Input files must be **paired-end, Phred+33, gzip-compressed FASTQ files**. Provide one pair per sample. For multiple sequencing lanes, concatenate R1 and R2 files in the same lane order and retain a record of the merge. Use only letters, numbers, underscores, and hyphens in sample IDs; avoid spaces in project paths.
- The target reference index, chromosome sizes, and blacklist must use the same genome assembly and chromosome naming convention. `chrom.sizes` is a two-column, tab-delimited file containing chromosome names and lengths; it is not an effective genome size value.
- The default **SEACR numeric threshold mode** follows the original workflow and does not use an IgG control for peak calling. A value of `0.01` selects approximately the top 1% of regions ranked by total signal; **it is not a P value or FDR**. If the experiment requires IgG control mode, define the sample-to-control assignments and modify the SEACR inputs first. Do not describe numeric-threshold results as control-adjusted peaks. [SEACR input documentation](https://github.com/FredHutch/SEACR#description-of-input-fields)
- bigWig files use **BPM** normalization by default. If the experimental design supports spike-in calibration, switch the entire batch to `spikein`. This enables spike-in alignment and applies `10000 / number of valid spike-in fragments` to both bigWig and bedGraph signals. A zero spike-in count stops the workflow without falling back to BPM. Signal values from different normalization modes should not be compared directly.
- As in the original workflow, **PCR duplicates are retained, and no additional MAPQ threshold is applied**. Properly paired primary alignments reported by Bowtie2 are retained, including some multimapping fragments. The resulting data are therefore not restricted to unique, deduplicated alignments. Any stricter filtering must be justified by the experimental design and applied consistently across relevant samples.

### Changes from the original workflow record

This template adds FastQC and MultiQC, explicitly records parameters and index/output dependencies, and removes `touch` commands that could create missing FASTQ or peak outputs. Spike-in fragments are counted from read 1 records in BAM files instead of fixed Bowtie2 log line numbers. The Bowtie2 random seed is fixed, bigWig coverage counts paired-end fragments, and basic plots use peaks remaining after blacklist filtering. The target reference is named `target` throughout so that the template can be used with different species.

## 2. Directory and software setup

### 2.1 Project directory

```bash
mkdir -p "[path_to_your_dir]/CUT_TAG_project"
cd "[path_to_your_dir]/CUT_TAG_project"
mkdir -p config workflow data/raw data/clean logs results
```

Prepare the input structure below and save the code from Sections 3–5 to the indicated files. Copy only the contents of each code block, excluding the Markdown fences.

```text
CUT_TAG_project/
├── config/
│   └── config.yaml
├── workflow/
│   ├── Snakefile
│   └── all_rules.smk
├── data/raw/
│   ├── sample_01/
│   │   ├── sample_01_R1.fq.gz
│   │   └── sample_01_R2.fq.gz
│   └── sample_02/
│       ├── sample_02_R1.fq.gz
│       └── sample_02_R2.fq.gz
├── logs/
└── results/
```

If your FASTQ files use different names, update the raw input patterns in both `fastqc_raw` and `trim_galore`, or create symbolic links matching this naming convention. There is no need to move or overwrite the original data.

### 2.2 Software environment

Prepare a Linux environment validated by your laboratory, containing at least the following software:

| Software | Purpose |
| --- | --- |
| Snakemake, Python, PyYAML | Workflow scheduling and configuration parsing |
| FastQC, Java, MultiQC | FASTQ quality control and report aggregation |
| Trim Galore, Cutadapt, Perl, pigz | Paired-end trimming and compression |
| Bowtie2, SAMtools, BEDTools | Alignment, BAM processing, and fragment coverage |
| deepTools | `bamCoverage`, `bamPEFragmentSize`, `computeMatrix`, `plotHeatmap`, `plotProfile` |
| SEACR 1.3, R / Rscript | Peak calling |
| Bash, awk, sort, coreutils | Shell pipelines, text processing, and file operations |

This template uses Snakemake 8.x syntax. After upgrading software, repeat the dry run and validate a representative sample. Set `software.bin` to the environment's `bin/` directory; Trim Galore and SEACR may be stored separately. Keep the SEACR `.sh` script and matching `.R` file in the same directory, and make both `Rscript` and `bedtools` available through PATH. [SEACR requirements](https://github.com/FredHutch/SEACR)

The template uses one existing environment and does not declare per-rule `conda:` environments. After the first successful run, record the exact environment versions as described in Section 6 so that others can reuse them.

## 3. Configuration file: `config/config.yaml`

```yaml
# These values are starting points. Review and freeze the configuration before analysis.
samples:
  - sample_01
  - sample_02

ref:
  # Bowtie2 index prefix, without .1.bt2; provide the prefix rather than the directory.
  target_index: "[path_to_your_reference]/target/bowtie2/genome"
  target_index_ext: "bt2"  # Use bt2l for large indexes; each index consists of six files.
  chrom_sizes: "[path_to_your_reference]/target/genome.chrom.sizes"
  blacklist: "[path_to_your_reference]/target/blacklist.bed"
  # The following index is required only when normalization.mode is spikein.
  spike_index: "[path_to_your_reference]/spikein/bowtie2/genome"
  spike_index_ext: "bt2"

software:
  bin: "[path_to_your_software]/bin"
  trim_galore: "[path_to_your_software]/TrimGalore/trim_galore"
  seacr: "[path_to_your_software]/SEACR/SEACR_1.3.sh"

threads:
  fastqc: 2      # Process two FASTQ files in parallel.
  trim: 8        # Budget for the entire Trim Galore job; the rule fixes --cores at 1.
  mapping: 8     # Bowtie2 uses n-1 threads; reserve one for SAMtools view.
  samtools: 4   # Total sorting/indexing budget; -@ uses n-1 additional threads.
  deeptools: 4

trim:
  quality: 20    # Trim terminal bases below this Phred quality threshold.
  length: 20     # Discard the read pair if either mate is too short after trimming.
  stringency: 1  # Minimum overlap with the adapter sequence, in bp.
  # Use Trim Galore adapter autodetection; specify nonstandard adapters separately.

mapping:
  min_insert: 10   # Bowtie2 -I: minimum permitted paired-end fragment length.
  max_insert: 700  # Bowtie2 -X: maximum fragment length; check against the library design.
  seed: 42        # Fixed seed; reproducibility also requires fixed software and inputs.
  # Spike-in alignment retains --no-overlap --no-dovetail from the original workflow.
  # These restrictions may exclude valid short spike-in fragments; check the library design.

normalization:
  mode: "BPM"      # Allowed values: BPM or spikein. Use one mode for the entire batch.
  spike_constant: 10000  # Scaling constant for spikein mode, not a sequencing depth cutoff.

coverage:
  bin_size: 25     # bigWig bin size in bp; smaller bins provide higher resolution.

peaks:
  top_fraction: 0.01  # SEACR numeric threshold: approximately the top 1% of regions, not FDR.
  # Use non stringent, as in the original workflow.

plot:
  upstream: 3000   # Distance upstream of the peak maximum-signal region center, in bp.
  downstream: 3000 # Downstream distance, in bp.
  bin_size: 25     # Bin size used by computeMatrix.
  dpi: 300        # Output image resolution.
```

`Trim Galore --cores` does not specify the total number of threads used by the entire job: Cutadapt and compression programs also consume threads. This template fixes `--cores 1` and reserves eight cores per trimming job to control concurrent trimming tasks. Before increasing parallelism, check the Trim Galore version and the scheduler's resource budget. [Trim Galore parameter documentation](https://github.com/FelixKrueger/TrimGalore)

If no suitable blacklist is available, explicitly create an empty BED file, point `ref.blacklist` to it, and record that no regions were excluded. Do not substitute a blacklist from another species or genome assembly.

## 4. Main workflow file: `workflow/Snakefile`

```python
import os
import re
from pathlib import Path

configfile: "config/config.yaml"

SAMPLES = config["samples"]
MODE = config["normalization"]["mode"]
SW = config["software"]
REF = config["ref"]
GENOMES = ["target"] + (["spike"] if MODE == "spikein" else [])

if not SAMPLES or len(SAMPLES) != len(set(SAMPLES)):
    raise ValueError("samples must be non-empty and unique")
if any(not re.fullmatch(r"[A-Za-z0-9_-]+", s) for s in SAMPLES):
    raise ValueError("Use only letters, digits, underscores and hyphens in sample IDs")
if MODE not in ("BPM", "spikein"):
    raise ValueError("normalization.mode must be BPM or spikein")
if not 0 < config["peaks"]["top_fraction"] < 1:
    raise ValueError("peaks.top_fraction must be between 0 and 1")
if config["normalization"]["spike_constant"] <= 0:
    raise ValueError("spike_constant must be positive")

# Use the same environment for external tools and their subprocesses to avoid mixing versions.
os.environ["PATH"] = SW["bin"] + os.pathsep + os.environ["PATH"]
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; export LC_ALL=C; ")

# Track all six Bowtie2 index files as dependencies so reference changes can trigger reruns.
INDEX_PARTS = ["1", "2", "3", "4", "rev.1", "rev.2"]
INDEX_FILES = {
    g: [f"{REF[g + '_index']}.{part}.{REF[g + '_index_ext']}"
        for part in INDEX_PARTS]
    for g in GENOMES
}

wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLES),
    genome="|".join(GENOMES)

rule all:
    input:
        "results/qc/multiqc_report.html",
        expand("results/{genome}/{sample}.bam.bai", genome=GENOMES, sample=SAMPLES),
        expand("results/qc/alignment/{genome}/{sample}.flagstat.txt",
               genome=GENOMES, sample=SAMPLES),
        expand("results/qc/fragment_length/{sample}.pdf", sample=SAMPLES),
        expand("results/qc/fragment_length/{sample}.tsv", sample=SAMPLES),
        expand("results/normalization/{sample}.tsv", sample=SAMPLES),
        expand("results/filtered_peaks/{sample}.bed", sample=SAMPLES),
        expand("results/bw/{sample}.bw", sample=SAMPLES),
        expand("results/plots/peaks/{sample}.heatmap.pdf", sample=SAMPLES),
        expand("results/plots/peaks/{sample}.profile.pdf", sample=SAMPLES)

include: "all_rules.smk"
```

`rule all` lists the final deliverables, allowing Snakemake to identify the required upstream jobs. The `include` path is relative to the Snakefile, while input and output paths are relative to the project working directory. For this reason, Section 6 explicitly specifies `--snakefile workflow/Snakefile`. [Snakemake rule documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)

## 5. Workflow rules: `workflow/all_rules.smk`

**Concatenate the Python code blocks in Sections 5.1–5.8, in order, into this single file.**

### 5.1 Raw quality control, trimming, and post-trimming quality control

FastQC generates HTML and ZIP reports for each mate. Trim Galore's `--paired` option preserves read pairing, while `--basename` fixes the output names. Success is determined by the program exit status and actual output files.

```python
rule fastqc_raw:
    input:
        r1="data/raw/{sample}/{sample}_R1.fq.gz",
        r2="data/raw/{sample}/{sample}_R2.fq.gz"
    output:
        html1="results/qc/fastqc_raw/{sample}_R1_fastqc.html",
        zip1="results/qc/fastqc_raw/{sample}_R1_fastqc.zip",
        html2="results/qc/fastqc_raw/{sample}_R2_fastqc.html",
        zip2="results/qc/fastqc_raw/{sample}_R2_fastqc.zip"
    threads: config["threads"]["fastqc"]
    log: "logs/fastqc_raw/{sample}.log"
    shell:
        r"""
        mkdir -p results/qc/fastqc_raw
        fastqc --threads {threads} --outdir results/qc/fastqc_raw \
            {input.r1:q} {input.r2:q} > {log:q} 2>&1
        """


rule trim_galore:
    input:
        r1="data/raw/{sample}/{sample}_R1.fq.gz",
        r2="data/raw/{sample}/{sample}_R2.fq.gz"
    output:
        r1="data/clean/{sample}_val_1.fq.gz",
        r2="data/clean/{sample}_val_2.fq.gz",
        report1="data/clean/{sample}_R1.fq.gz_trimming_report.txt",
        report2="data/clean/{sample}_R2.fq.gz_trimming_report.txt"
    threads: config["threads"]["trim"]
    params:
        tool=SW["trim_galore"],
        quality=config["trim"]["quality"],
        length=config["trim"]["length"],
        stringency=config["trim"]["stringency"]
    log: "logs/trim/{sample}.log"
    shell:
        r"""
        mkdir -p data/clean
        {params.tool:q} --paired --phred33 --cores 1 \
            --quality {params.quality} --length {params.length} \
            --stringency {params.stringency} --gzip \
            --basename {wildcards.sample:q} --output_dir data/clean \
            {input.r1:q} {input.r2:q} > {log:q} 2>&1
        """


rule fastqc_clean:
    input:
        r1="data/clean/{sample}_val_1.fq.gz",
        r2="data/clean/{sample}_val_2.fq.gz"
    output:
        html1="results/qc/fastqc_clean/{sample}_val_1_fastqc.html",
        zip1="results/qc/fastqc_clean/{sample}_val_1_fastqc.zip",
        html2="results/qc/fastqc_clean/{sample}_val_2_fastqc.html",
        zip2="results/qc/fastqc_clean/{sample}_val_2_fastqc.zip"
    threads: config["threads"]["fastqc"]
    log: "logs/fastqc_clean/{sample}.log"
    shell:
        r"""
        mkdir -p results/qc/fastqc_clean
        fastqc --threads {threads} --outdir results/qc/fastqc_clean \
            {input.r1:q} {input.r2:q} > {log:q} 2>&1
        """
```

### 5.2 Target genome / spike-in alignment, BAM sorting, and indexing

`--local --very-sensitive-local` enables highly sensitive local alignment. `--no-mixed --no-discordant` disables single-mate rescue and discordant pairing. `samtools view -f 2 -F 2308` retains properly paired alignments and excludes unmapped, secondary, and supplementary records. Spike-in alignment additionally disallows overlapping pairs, following the original settings. The two references are aligned separately rather than competitively; assess potential cross-mapping between the spike-in and target references before use. [Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), [SAMtools view parameters](https://www.htslib.org/doc/samtools-view.html)

```python
rule bowtie2_mapping:
    input:
        r1="data/clean/{sample}_val_1.fq.gz",
        r2="data/clean/{sample}_val_2.fq.gz",
        index=lambda wc: INDEX_FILES[wc.genome]
    output:
        bam=temp("results/{genome}/{sample}.unsorted.bam"),
        summary="results/qc/alignment/{genome}/{sample}.bowtie2.log"
    threads: config["threads"]["mapping"]
    params:
        index=lambda wc: REF[wc.genome + "_index"],
        bt_threads=lambda wc, threads: max(1, threads - 1),
        min_insert=config["mapping"]["min_insert"],
        max_insert=config["mapping"]["max_insert"],
        seed=config["mapping"]["seed"],
        extra=lambda wc: "--no-overlap --no-dovetail" if wc.genome == "spike" else ""
    log: "logs/mapping/{genome}/{sample}.samtools.log"
    shell:
        r"""
        bowtie2 --local --very-sensitive-local --no-mixed --no-discordant \
            --phred33 --seed {params.seed} {params.extra} \
            -I {params.min_insert} -X {params.max_insert} \
            -p {params.bt_threads} -x {params.index:q} \
            -1 {input.r1:q} -2 {input.r2:q} 2> {output.summary:q} \
          | samtools view -u -f 2 -F 2308 -o {output.bam:q} - 2> {log:q}
        """


rule sort_index_bam:
    input:
        bam="results/{genome}/{sample}.unsorted.bam"
    output:
        bam="results/{genome}/{sample}.bam",
        bai="results/{genome}/{sample}.bam.bai",
        stats="results/qc/alignment/{genome}/{sample}.flagstat.txt"
    threads: config["threads"]["samtools"]
    params:
        extra_threads=lambda wc, threads: max(0, threads - 1)
    log: "logs/sort_index/{genome}/{sample}.log"
    shell:
        r"""
        samtools sort -@ {params.extra_threads} -o {output.bam:q} {input.bam:q} 2> {log:q}
        samtools index -@ {params.extra_threads} {output.bam:q} {output.bai:q} 2>> {log:q}
        samtools flagstat -@ {params.extra_threads} {output.bam:q} > {output.stats:q} 2>> {log:q}
        """
```

`*.bowtie2.log` summarizes alignment relative to all input reads. `*.flagstat.txt` describes the **filtered BAM**, so its percentages must not be interpreted as overall alignment rates for the raw data. Snakemake removes temporary unsorted BAM files after their downstream jobs finish.

### 5.3 Record the normalization mode and spike-in fragment count

`-f 66` requires both proper pairing (2) and read 1 (64), counting each fragment once. BPM mode does not require a spike-in index: bedGraph retains raw fragment coverage, and BPM normalization applies only to bigWig. Section 7 describes the units of each output.

```python
rule normalization:
    input:
        spike=lambda wc: [f"results/spike/{wc.sample}.bam"] if MODE == "spikein" else []
    output:
        "results/normalization/{sample}.tsv"
    params:
        mode=MODE,
        constant=config["normalization"]["spike_constant"]
    run:
        import subprocess

        count = "NA"
        scale = 1.0
        if params.mode == "spikein":
            count = int(subprocess.check_output(
                ["samtools", "view", "-c", "-f", "66", "-F", "2308", input.spike[0]],
                text=True,
            ).strip())
            if count == 0:
                raise ValueError(f"{wildcards.sample}: zero spike-in fragments; review the sample and reference")
            scale = float(params.constant) / count
        Path(output[0]).write_text(
            f"mode\tspike_fragments\tscale_factor\n{params.mode}\t{count}\t{scale:.12g}\n",
            encoding="utf-8",
        )
```

A nonzero spike-in count allows calculation of a scaling factor but does not establish that calibration is reliable. Interpret the count alongside the spike-in addition procedure and amount, sequencing depth, and alignment results. Very few spike-in fragments can produce unstable scaling.

### 5.4 Generate paired-end fragment BED and bedGraph files

Sort the BAM by read name, then use BEDPE to represent each mate pair as one fragment. For each pair, take the smaller start coordinate and the larger end coordinate, retaining only positive-length intervals on the same chromosome. Bowtie2 already controls the permitted fragment length range. `genomecov -bg` reports only regions with nonzero signal, as required by SEACR. [BEDTools BEDPE documentation](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html), [SEACR bedGraph preparation](https://github.com/FredHutch/SEACR#preparing-input-bedgraph-files)

```python
rule fragment_bed:
    input:
        bam="results/target/{sample}.bam"
    output:
        bed="results/fragments/{sample}.bed"
    # The pipeline runs SAMtools, BEDTools, awk, and sort; budget for at least four processes.
    threads: config["threads"]["samtools"]
    log: "logs/fragments/{sample}.log"
    shell:
        r"""
        (
          samtools sort -n -@ 0 -T {output.bed:q}.sorttmp -O BAM {input.bam:q} \
            | bedtools bamtobed -bedpe -i stdin \
            | awk 'BEGIN{{OFS="\t"}} $1==$4 && $2>=0 && $5>=0 {{
                start=($2<$5 ? $2 : $5); end=($3>$6 ? $3 : $6);
                if (end>start) print $1,start,end
              }}' \
            | sort --parallel=1 -k1,1 -k2,2n -k3,3n > {output.bed:q}
          test -s {output.bed:q}
        ) 2> {log:q}
        """


rule fragment_bedgraph:
    input:
        bed="results/fragments/{sample}.bed",
        norm="results/normalization/{sample}.tsv",
        sizes=REF["chrom_sizes"]
    output:
        bg="results/bedgraph/{sample}.bedgraph"
    log: "logs/bedgraph/{sample}.log"
    shell:
        r"""
        scale=$(awk 'NR==2 {{print $3}}' {input.norm:q})
        bedtools genomecov -bg -scale "$scale" -i {input.bed:q} -g {input.sizes:q} \
            > {output.bg:q} 2> {log:q}
        """
```

### 5.5 SEACR peak calling and blacklist filtering

`non` skips SEACR control normalization, and `stringent` retains the original calling mode. `intersect -v` removes **entire peaks** that overlap the blacklist, without splitting intervals. All six SEACR columns are preserved so that the maximum-signal region can later be extracted from column 6. The rule avoids `grep chr`, which could incorrectly discard chromosomes named using conventions such as `1`.

```python
rule seacr:
    input:
        bg="results/bedgraph/{sample}.bedgraph",
        blacklist=REF["blacklist"],
        sh=SW["seacr"],
        r=str(Path(SW["seacr"]).with_suffix(".R"))
    output:
        raw="results/peaks/{sample}.stringent.bed",
        filtered="results/filtered_peaks/{sample}.bed"
    params:
        fraction=config["peaks"]["top_fraction"],
        prefix="results/peaks/{sample}"
    log: "logs/seacr/{sample}.log"
    shell:
        r"""
        bash {input.sh:q} {input.bg:q} {params.fraction} non stringent {params.prefix:q} > {log:q} 2>&1
        bedtools intersect -v -a {output.raw:q} -b {input.blacklist:q} > {output.filtered:q} 2>> {log:q}
        """
```

### 5.6 Generate bigWig files

`--extendReads --samFlagInclude 64` uses read 1 to represent each complete paired-end fragment, avoiding double counting of mates. BPM normalizes by the total signal across bins. Spike-in mode uses `None` and the recorded scaling factor without additional BPM normalization. Neither BPM nor simple scaling requires `--effectiveGenomeSize`, so the original fixed human value is omitted. [bamCoverage parameter documentation](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)

```python
rule bam_to_bigwig:
    input:
        bam="results/target/{sample}.bam",
        bai="results/target/{sample}.bam.bai",
        norm="results/normalization/{sample}.tsv",
        blacklist=REF["blacklist"]
    output:
        bw="results/bw/{sample}.bw"
    threads: config["threads"]["deeptools"]
    params:
        mode=MODE,
        bins=config["coverage"]["bin_size"]
    log: "logs/bigwig/{sample}.log"
    shell:
        r"""
        scale=$(awk 'NR==2 {{print $3}}' {input.norm:q})
        norm=None
        if [ "{params.mode}" = "BPM" ]; then norm=BPM; fi
        bamCoverage --bam {input.bam:q} --outFileName {output.bw:q} \
            --numberOfProcessors {threads} --binSize {params.bins} \
            --extendReads --samFlagInclude 64 \
            --normalizeUsing "$norm" --scaleFactor "$scale" \
            --blackListFileName {input.blacklist:q} > {log:q} 2>&1
        """
```

bigWig generation uses deepTools' chunk-based blacklist filtering, whereas peak filtering uses BEDTools to exclude entire overlapping intervals. Their boundary behavior differs: bigWig output should not be interpreted as precisely removing every fragment that spans a blacklisted region. [deepTools blacklist behavior](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html#optional-arguments)

### 5.7 Basic plots: fragment lengths, peak heatmaps, and average profiles

deepTools estimates fragment lengths by sampling the filtered BAM and produces a PDF, a statistics table, and a summary in the log. This is not an exact frequency table obtained by traversing every fragment. [bamPEFragmentSize documentation](https://deeptools.readthedocs.io/en/develop/content/tools/bamPEFragmentSize.html)

SEACR column 6 contains the **maximum-signal region**, not a single-base summit. This template plots 3 kb upstream and downstream of that region's center. `--skipZeros` removes rows with zero signal throughout the window. If no peaks remain after filtering, plotting stops with an explanatory log message rather than creating placeholder plots.

```python
rule fragment_length_plot:
    input:
        bam="results/target/{sample}.bam",
        bai="results/target/{sample}.bam.bai"
    output:
        pdf="results/qc/fragment_length/{sample}.pdf",
        table="results/qc/fragment_length/{sample}.tsv"
    threads: config["threads"]["deeptools"]
    log: "logs/fragment_length/{sample}.log"
    shell:
        r"""
        bamPEFragmentSize --bamfiles {input.bam:q} \
            --histogram {output.pdf:q} --table {output.table:q} \
            --samplesLabel {wildcards.sample:q} --numberOfProcessors {threads} \
            > {log:q} 2>&1
        """


rule peak_matrix:
    input:
        peaks="results/filtered_peaks/{sample}.bed",
        bw="results/bw/{sample}.bw"
    output:
        centers="results/plots/peaks/{sample}.max_signal_regions.bed",
        matrix="results/plots/peaks/{sample}.matrix.gz"
    threads: config["threads"]["deeptools"]
    params:
        before=config["plot"]["upstream"],
        after=config["plot"]["downstream"],
        bins=config["plot"]["bin_size"]
    log: "logs/peak_matrix/{sample}.log"
    shell:
        r"""
        if [ ! -s {input.peaks:q} ]; then
            echo "No filtered peaks; review QC and SEACR settings before plotting." > {log:q}
            exit 1
        fi
        awk 'BEGIN{{OFS="\t"}} {{split($6,a,":"); split(a[2],b,"-"); print a[1],b[1],b[2]}}' \
            {input.peaks:q} > {output.centers:q}
        computeMatrix reference-point --referencePoint center \
            -S {input.bw:q} -R {output.centers:q} \
            --beforeRegionStartLength {params.before} --afterRegionStartLength {params.after} \
            --binSize {params.bins} --skipZeros --numberOfProcessors {threads} \
            --outFileName {output.matrix:q} > {log:q} 2>&1
        """


rule peak_plots:
    input:
        matrix="results/plots/peaks/{sample}.matrix.gz"
    output:
        heatmap="results/plots/peaks/{sample}.heatmap.pdf",
        profile="results/plots/peaks/{sample}.profile.pdf"
    params:
        dpi=config["plot"]["dpi"]
    log: "logs/peak_plots/{sample}.log"
    shell:
        r"""
        plotHeatmap --matrixFile {input.matrix:q} --outFileName {output.heatmap:q} \
            --dpi {params.dpi} --regionsLabel Peaks --samplesLabel {wildcards.sample:q} \
            --boxAroundHeatmaps no --colorList white,red --missingDataColor white > {log:q} 2>&1
        plotProfile --matrixFile {input.matrix:q} --outFileName {output.profile:q} \
            --dpi {params.dpi} --regionsLabel Peaks --samplesLabel {wildcards.sample:q} \
            --refPointLabel "Peak maximum-signal center" >> {log:q} 2>&1
        """
```

These plots assess signal patterns within individual samples. Each sample uses its own peaks and an automatic color scale, so the plots do not establish quantitative differences between groups. Customized group plots around TSSs or gene bodies are outside the scope of this template.

### 5.8 Aggregate basic quality control reports

MultiQC waits for all FASTQ quality control results, trimming reports, and alignment logs before aggregating them. Fragment length and peak plots are saved separately and do not depend on whether MultiQC recognizes their formats.

```python
rule multiqc:
    input:
        raw=expand("results/qc/fastqc_raw/{sample}_R{mate}_fastqc.zip",
                   sample=SAMPLES, mate=[1, 2]),
        clean=expand("results/qc/fastqc_clean/{sample}_val_{mate}_fastqc.zip",
                     sample=SAMPLES, mate=[1, 2]),
        trim=expand("data/clean/{sample}_R{mate}.fq.gz_trimming_report.txt",
                    sample=SAMPLES, mate=[1, 2]),
        mapping=expand("results/qc/alignment/{genome}/{sample}.bowtie2.log",
                       genome=GENOMES, sample=SAMPLES),
        stats=expand("results/qc/alignment/{genome}/{sample}.flagstat.txt",
                     genome=GENOMES, sample=SAMPLES)
    output:
        html="results/qc/multiqc_report.html",
        data=directory("results/qc/multiqc_data")
    log: "logs/multiqc.log"
    shell:
        r"""
        multiqc results/qc/fastqc_raw results/qc/fastqc_clean \
            results/qc/alignment data/clean --dirs --dirs-depth 1 \
            --outdir results/qc --filename multiqc_report.html --force > {log:q} 2>&1
        """
```

`--dirs --dirs-depth 1` adds the containing directory name to each sample label, preventing reports with identical names, such as target and spike-in alignment reports, from overwriting each other during aggregation.

## 6. Execution and reproducibility records

### 6.1 Activate the fixed environment and check software availability

```bash
cd "[path_to_your_dir]/CUT_TAG_project"
# If using Conda, activate the validated analysis environment first.
conda activate "[path_to_your_software_environment]"
export PATH="[path_to_your_software]/bin:$PATH"

for tool in snakemake fastqc multiqc cutadapt bowtie2 samtools bedtools \
            bamCoverage bamPEFragmentSize computeMatrix plotHeatmap plotProfile \
            Rscript java perl pigz; do
    command -v "$tool" || exit 1
done
test -x "[path_to_your_software]/TrimGalore/trim_galore"
test -f "[path_to_your_software]/SEACR/SEACR_1.3.sh"
test -f "[path_to_your_software]/SEACR/SEACR_1.3.R"
```

Software paths in the run commands and YAML configuration must refer to the same environment. If you do not use Conda, omit `conda activate` and the Conda export command below, and retain an installation manifest or container image digest for the actual environment.

### 6.2 Check dependencies before running the workflow

```bash
# -n checks and displays the plan without running analysis; -p prints the shell commands.
snakemake --snakefile workflow/Snakefile --cores 16 --dry-run --printshellcmds

# Run the workflow. --cores limits concurrent core usage; match it to your allocation.
# --rerun-incomplete reruns interrupted jobs whose outputs were not completed.
snakemake --snakefile workflow/Snakefile --cores 16 \
    --rerun-incomplete --printshellcmds
```

For background execution, use the following command **instead of** the foreground execution command above. Do not launch both workflow instances at once:

```bash
nohup snakemake --snakefile workflow/Snakefile --cores 16 \
    --rerun-incomplete --printshellcmds > logs/snakemake.run.log 2>&1 &
```

The core limit does not constrain memory usage. Reference indexes, BAM sorting, and matrix calculations also require sufficient memory and temporary disk space. On shared servers, run within your allocated compute node. This template assumes a budget of at least eight cores; setting a very small `--cores` value may not adequately account for external subprocesses.

### 6.3 Preserve the run's provenance

Replace the placeholders before archiving the configuration and rules. If you change parameters, create a separate project directory for the new run and retain the previous results and records.

```bash
mkdir -p results/provenance
cp config/config.yaml results/provenance/config.used.yaml
cp workflow/Snakefile workflow/all_rules.smk results/provenance/

# Record exact installed versions and builds; requires the complete analysis Conda environment.
conda list --explicit > results/provenance/conda-explicit.txt
snakemake --version > results/provenance/snakemake.version.txt
"[path_to_your_software]/TrimGalore/trim_galore" --version \
    > results/provenance/trim_galore.version.txt 2>&1
sha256sum "[path_to_your_software]/SEACR/SEACR_1.3.sh" \
    "[path_to_your_software]/SEACR/SEACR_1.3.R" \
    > results/provenance/seacr.sha256

# Record SHA-256 checksums for the configured FASTQ, reference index, blacklist, and chrom.sizes files.
# Reading large files takes time but does not modify inputs; BPM mode excludes unused spike indexes.
python - <<'PY'
import hashlib
from pathlib import Path
import yaml

c = yaml.safe_load(Path("config/config.yaml").read_text())
paths = [Path(f"data/raw/{s}/{s}_R{mate}.fq.gz")
         for s in c["samples"] for mate in (1, 2)]
paths += [Path(c["ref"][key]) for key in ("chrom_sizes", "blacklist")]
genomes = ["target"] + (["spike"] if c["normalization"]["mode"] == "spikein" else [])
for genome in genomes:
    prefix = c["ref"][genome + "_index"]
    ext = c["ref"][genome + "_index_ext"]
    paths += [Path(f"{prefix}.{part}.{ext}") for part in ("1", "2", "3", "4", "rev.1", "rev.2")]
with open("results/provenance/inputs.sha256", "w") as out:
    for path in paths:
        digest = hashlib.sha256()
        with path.open("rb") as source:
            for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
                digest.update(block)
        out.write(f"{digest.hexdigest()}  {path}\n")
PY
```

If you change the FASTQ naming pattern, update the input manifest code above as well. Before rerunning, verify inputs with `sha256sum -c results/provenance/inputs.sha256`. If directories change, update the paths in the checksum manifest accordingly. For an external Trim Galore installation, also retain the program directory for the version used.

## 7. Output files and completion checks

In the table below, `{sample}` represents a sample ID, and `{genome}` is `target` or `spike`. Files associated with `spike/` are generated only in spike-in mode.

| Output | Contents / checks |
| --- | --- |
| `results/qc/multiqc_report.html` | Summary of raw and trimmed FASTQ quality, adapters, and alignment; inspect for sample-specific anomalies |
| `data/clean/{sample}_val_1.fq.gz`, `_val_2.fq.gz` | Trimmed paired-end reads |
| `results/target/{sample}.bam`, `.bam.bai` | Coordinate-sorted BAM of properly paired primary alignments and its index |
| `results/spike/{sample}.bam`, `.bam.bai` | Spike-in alignment BAM and index |
| `results/qc/alignment/{genome}/{sample}.bowtie2.log` | Original Bowtie2 alignment summary |
| `results/qc/alignment/{genome}/{sample}.flagstat.txt` | Statistics for the filtered BAM |
| `results/normalization/{sample}.tsv` | Actual mode, spike-in fragment count, and scaling factor; the count is NA in BPM mode |
| `results/fragments/{sample}.bed` | Fragment intervals, with one row per read pair |
| `results/bedgraph/{sample}.bedgraph` | Raw fragment coverage in BPM mode or calibrated coverage in spike-in mode; input to SEACR |
| `results/peaks/{sample}.stringent.bed` | Original six-column SEACR peaks |
| `results/filtered_peaks/{sample}.bed` | Peaks remaining after removing blacklist-overlapping intervals |
| `results/bw/{sample}.bw` | BPM-normalized or spike-in-scaled signal in 25 bp bins, viewable in a genome browser |
| `results/qc/fragment_length/{sample}.pdf`, `.tsv` | Sampled fragment length plot and summary statistics |
| `results/plots/peaks/{sample}.max_signal_regions.bed` | Maximum-signal regions within the filtered peaks |
| `results/plots/peaks/{sample}.matrix.gz` | Signal matrix around peaks |
| `results/plots/peaks/{sample}.heatmap.pdf`, `.profile.pdf` | Per-sample peak heatmap and average signal profile |
| `results/provenance/`, `logs/`, `.snakemake/` | Configuration, input checksums, environment records, execution logs, and workflow state |

After completion, repeat the dry run. It should report that no jobs need to run. This confirms that the required outputs and dependencies are satisfied; biological quality control still requires review.

```bash
snakemake --snakefile workflow/Snakefile --cores 16 --dry-run
samtools quickcheck -v results/target/*.bam
```

Manually inspect whether trimming improves adapter contamination and read quality, whether target alignment rates and valid fragment counts are reasonable, whether spike-in counts support calibration, whether fragment length distributions are abnormal, whether peaks are empty, and whether the bigWig tracks and plots show the expected enrichment. Choose thresholds according to the experimental design and sample type; a single percentage is not a universal pass criterion for all CUT&Tag experiments.

If peaks are empty, the peak plotting rule stops. Review `logs/seacr/`, the alignment logs, and blacklist compatibility to distinguish a sample with little signal from a parameter or reference error. Do not create empty output files to bypass the checks.

> This is a reusable workflow template. Before applying it to an experiment, replace the placeholders, freeze the actual software environment, and complete a dry run, a representative sample run, and output inspection using real FASTQ and reference files.
