# WES Snakemake Upstream Analysis Pipeline Template

This template processes **paired-end whole-exome sequencing (WES)** data in either **tumor-only** or **matched tumor-normal** mode using **fastp → BWA-MEM → duplicate removal → BQSR → Mutect2 → FilterMutectCalls**. It produces analysis-ready BAM files, filtered candidate variant calls, and basic quality control reports. Run all commands in **Linux with Bash**, from the project root directory.

The workflow includes adapter trimming, read quality filtering, alignment, duplicate handling, base quality score recalibration (BQSR), capture coverage QC, contamination estimation, read-orientation artifact modeling, variant filtering, and PASS-only VCF export. Basic plots are available in the fastp and MultiQC HTML reports. Functional annotation, MAF conversion, mutation burden calculations, group comparisons, clinical interpretation, CNV/SV analysis, and other project-specific analyses are outside its scope.

## 1. Assumptions and analysis boundaries

- Replace every `[path_to_your_...]` placeholder with a real absolute path and replace the example sample IDs. Use only letters, digits, underscores, and hyphens in sample IDs. Avoid spaces in project and software paths.
- Provide one pair of **gzip-compressed, Phred+33 FASTQ files per tumor or normal sample**, representing one non-UMI library. If multiple lanes belong to the same library, merge R1 and R2 in the same lane order and preserve that record. Multiple libraries or samples with different capture kits require corresponding read-group or interval handling before using this template.
- Choose the calling mode through the `tumors` mapping in Section 3. A `null` normal selects **tumor-only** calling; a normal sample ID selects **matched tumor-normal** calling. A matched normal must be an appropriate non-tumor sample from the same individual. An unrelated normal, a PoN, and a known-sites VCF are not substitutes for that individual's matched normal.
- Both modes retain the existing population germline resource and PoN. PASS records remain **candidate somatic SNVs/indels** requiring QC; residual germline variation is a particular limitation of tumor-only calling. This is not a germline variant calling workflow. [Mutect2 documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360045800552-Mutect2)
- All reference resources must use the same assembly, contig names, and compatible sequence dictionaries. A label such as `GRCh38` alone does not guarantee identical reference sequences or contig sets. The BWA index must have been built from the exact FASTA used by GATK.
- Supply bait and target intervals for the actual capture kit. The template calls variants in **targets plus the configured padding**, while capture QC and contamination estimation use unpadded targets. The resulting PASS VCF therefore includes the padded calling territory.
- The original duplicate policy is retained: `MarkDuplicates --REMOVE_DUPLICATES true` **removes duplicates from its output BAM**. Duplicate metrics describe the input before removal. This is not UMI consensus processing. [MarkDuplicates documentation](https://gatk.broadinstitute.org/hc/en-us/articles/21905036102043-MarkDuplicates-Picard)
- No project-specific allele fraction, depth, gene, or chromosome filter is added after FilterMutectCalls. The PASS export selects the VCF FILTER field with bcftools and preserves headers and FORMAT definitions.

### Reproducibility improvements included here

The template tracks actual BWA index components, BAM/VCF indexes, and Mutect2 statistics as dependencies. It uses capture intervals, adds capture QC and read-orientation artifact filtering, records explicit preprocessing parameters, and lets Snakemake manage intermediate BAM files. It also supports matched normals through the same preprocessing rules. It does not create dummy reference files or remove inputs from inside downstream shell commands. The original toolchain and tumor-only option are retained, but results are not expected to be numerically identical to the original unbounded, less fully specified run.

## 2. Project and software setup

### 2.1 Directory structure

```bash
mkdir -p "[path_to_your_dir]/WES_project"
cd "[path_to_your_dir]/WES_project"
mkdir -p config workflow data/raw data/clean results logs
```

Prepare this input layout. Save the code in Sections 3–5 to the indicated files, excluding the Markdown fences.

```text
WES_project/
├── config/
│   └── config.yaml
├── workflow/
│   ├── Snakefile
│   └── all_rules.smk
├── data/raw/
│   ├── tumor_01_1.fq.gz
│   ├── tumor_01_2.fq.gz
│   ├── tumor_02_1.fq.gz
│   └── tumor_02_2.fq.gz
├── data/clean/
├── results/
└── logs/
```

If your input filenames differ, update the raw input pattern in `rule fastp` and the checksum manifest in Section 6, or create symbolic links that follow this naming convention. Keep the original FASTQ files unchanged.

### 2.2 Required software

| Software | Purpose |
| --- | --- |
| Snakemake 8.x, Python, PyYAML | Workflow scheduling, configuration, and input manifests |
| fastp | Adapter trimming, read filtering, and before/after QC plots |
| Classic BWA with the `mem` subcommand | Alignment using `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa` index files |
| SAMtools | Sorting, BAM indexing, and alignment QC |
| GATK 4.1.1 or later and the Java version required by that release | Duplicate removal, recalibration, capture metrics, and variant calling/filtering |
| bcftools | PASS selection, VCF indexing, and variant summary statistics |
| MultiQC | Combined HTML QC report and basic plots |
| Bash and coreutils | Shell execution and file checksums |

Use a laboratory-validated environment and freeze exact versions after validation. `software.bin` must contain the command-line tools and a compatible Java executable. `software.gatk` points to the GATK launcher. This template uses one existing environment and does not declare per-rule `conda:` environments. A BWA-MEM2 index cannot be substituted for the classic BWA index.

### 2.3 Reference preparation

Use existing, verified reference files. Prepare the following before the dry run:

- FASTA, its `.fai`, and the corresponding sequence dictionary.
- Five classic BWA index components sharing the prefix in `ref.bwa_index`.
- BGZF-compressed known-sites, germline, PoN, and common biallelic SNP VCFs with matching `.tbi` indexes.
- Capture bait and target files in Picard `.interval_list` format, with the reference dictionary in their headers.

The common SNP resource used for contamination estimation must contain suitable biallelic SNPs with population allele frequencies in `INFO/AF`. The germline resource also needs population allele frequencies; a variant database without the required annotations is not an interchangeable substitute.

If the capture provider supplies BED files, convert them with the dictionary matching your FASTA. BED coordinates are zero-based and half-open; let the conversion tool handle the coordinate convention.

```bash
mkdir -p reference/capture
"[path_to_your_software]/gatk/gatk" BedToIntervalList \
    --INPUT "[path_to_your_reference]/capture/targets.bed" \
    --SEQUENCE_DICTIONARY "[path_to_your_reference]/genome.dict" \
    --OUTPUT reference/capture/targets.interval_list
"[path_to_your_software]/gatk/gatk" BedToIntervalList \
    --INPUT "[path_to_your_reference]/capture/baits.bed" \
    --SEQUENCE_DICTIONARY "[path_to_your_reference]/genome.dict" \
    --OUTPUT reference/capture/baits.interval_list
```

Then set the two capture paths in the configuration to these outputs. Bait intervals describe capture probes, whereas target intervals describe the regions intended for coverage. Do not silently treat them as interchangeable. [CollectHsMetrics input requirements](https://gatk.broadinstitute.org/hc/en-us/articles/360036725371-CollectHsMetrics-Picard)

## 3. Configuration file: `config/config.yaml`

### 3.1 Default configuration: tumor-only calling

```yaml
# Replace all paths and sample IDs, then freeze the configuration for the run.
# List every sample that needs preprocessing, including normals when used.
samples:
  - tumor_01
  - tumor_02

# Keys are tumor IDs; values are matched normal IDs or YAML null for tumor-only.
# Only the keys in this mapping receive variant callsets.
tumors:
  tumor_01: null
  tumor_02: null

software:
  bin: "[path_to_your_software]/bin"
  gatk: "[path_to_your_software]/gatk/gatk"

ref:
  fasta: "[path_to_your_reference]/genome.fa"
  dict: "[path_to_your_reference]/genome.dict"
  # Prefix only; the five actual index files are declared in the Snakefile.
  bwa_index: "[path_to_your_reference]/bwa/genome"
  known_sites:
    - "[path_to_your_reference]/known_sites/known_snps.vcf.gz"
    - "[path_to_your_reference]/known_sites/high_confidence_snps.vcf.gz"
    - "[path_to_your_reference]/known_sites/known_indels.vcf.gz"
  germline: "[path_to_your_reference]/resources/af_only_germline.vcf.gz"
  pon: "[path_to_your_reference]/resources/panel_of_normals.vcf.gz"
  common_snps: "[path_to_your_reference]/resources/common_biallelic_snps.vcf.gz"
  targets: "[path_to_your_reference]/capture/targets.interval_list"
  baits: "[path_to_your_reference]/capture/baits.interval_list"

threads:
  fastp: 4
  mapping: 8    # BWA receives n-1 threads; reserve one for SAMtools view.
  samtools: 4  # -@ receives n-1 additional threads.
  gatk: 2     # Also passed to Java's ActiveProcessorCount.
  mutect2: 4  # PairHMM receives n-1 threads, leaving room for the main process.

memory:
  fastp_mb: 4000
  mapping_mb: 12000
  sort_mb: 8000
  gatk_mb: 12000   # Scheduler reservation, including non-heap overhead.
  java_heap_mb: 8000  # Java -Xmx; keep lower than gatk_mb.
  multiqc_mb: 4000

fastp:
  qualified_quality_phred: 15  # Bases below Q15 count as unqualified.
  unqualified_percent_limit: 40  # Maximum allowed percentage of unqualified bases.
  n_base_limit: 5            # Maximum number of N bases per read.
  length_required: 15       # Minimum read length after adapter trimming.
  # PE adapter detection is enabled explicitly in the rule.
  # No sliding-window quality trimming or FASTQ-level deduplication is enabled.

analysis:
  interval_padding: 100  # bp added around targets for BQSR training and Mutect2.
  # ApplyBQSR transforms the entire BAM; QC/pileup targets are not padded.

capture_qc:
  min_mapping_quality: 20  # Minimum MAPQ for target coverage metrics.
  min_base_quality: 20     # Minimum recalibrated base quality for coverage metrics.
  coverage_cap: 1000       # Depth cap for theoretical sensitivity calculations.
  # These are QC settings, not extra Mutect2 or PASS-selection filters.
```

The fastp quality parameters filter reads based on their base qualities; they do not enable sliding-window removal of low-quality ends. Adapter trimming is enabled, and paired-end adapter autodetection is requested explicitly. Review fastp's report for the actual adapters and retained reads. [fastp options](https://github.com/OpenGene/fastp)

### 3.2 Configuration with known matched normal samples

For two matched pairs, replace only the `samples` and `tumors` blocks above with the following. Keep the remaining configuration unchanged. Do not append duplicate YAML keys.

```yaml
samples:
  - tumor_01
  - normal_01
  - tumor_02
  - normal_02

tumors:
  tumor_01: normal_01
  tumor_02: normal_02
```

Add `normal_01_1.fq.gz`, `normal_01_2.fq.gz`, `normal_02_1.fq.gz`, and `normal_02_2.fq.gz` to `data/raw/`. Each normal follows the same trimming, alignment, duplicate removal, recalibration, and QC rules. Variant calling remains one job per tumor, with its matched normal supplied to that job. The normal does not receive a separate somatic callset.

The mapping identifies **sample IDs**, not BAM paths or VCF resources. Mutect2's `-normal` argument must match the normal BAM's read-group `SM` value. The BWA rule sets `SM` from the configured sample ID, so `tumor_01: normal_01` produces a call with the tumor BAM, the normal BAM, and `-normal normal_01`. [Mutect2 normal-sample argument](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)

Use unquoted YAML `null` when no normal exists; the string `"null"` would be treated as a sample name. Pairings are explicit: the workflow does not infer them from filename similarities. One normal can be reused for multiple tumors from the same individual when that relationship is appropriate and recorded. Confirm sample identity and normal suitability before assigning pairs.

## 4. Main workflow: `workflow/Snakefile`

```python
import os
import re
import shlex

configfile: "config/config.yaml"

SAMPLES = config["samples"]
PAIRS = config["tumors"]
TUMORS = list(PAIRS)
SW = config["software"]
REF = config["ref"]
MEM = config["memory"]

if not SAMPLES or len(SAMPLES) != len(set(SAMPLES)):
    raise ValueError("Provide a non-empty list of unique sample IDs")
if any(not isinstance(s, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", s) for s in SAMPLES):
    raise ValueError("Sample IDs may contain only letters, digits, underscores, and hyphens")
if not TUMORS or any(t not in SAMPLES for t in TUMORS):
    raise ValueError("Each tumor must be listed in samples; provide at least one tumor")
for tumor, normal in PAIRS.items():
    if normal is not None and normal not in SAMPLES:
        raise ValueError(f"{tumor}: matched normal must be listed in samples, or use YAML null")
    if normal == tumor:
        raise ValueError(f"{tumor}: tumor and matched normal must be different samples")
if not REF["known_sites"]:
    raise ValueError("At least one compatible known-sites VCF is required for BQSR")
if config["analysis"]["interval_padding"] < 0:
    raise ValueError("interval_padding must be non-negative")
if not 0 < MEM["java_heap_mb"] < MEM["gatk_mb"]:
    raise ValueError("java_heap_mb must be positive and smaller than gatk_mb")

os.environ["PATH"] = SW["bin"] + os.pathsep + os.environ["PATH"]
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; export LC_ALL=C; ")

# Declare real reference files rather than treating the BWA prefix as a file.
BWA_INDEX = [REF["bwa_index"] + suffix for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")]
REF_SUPPORT = [REF["fasta"] + ".fai", REF["dict"]]
KNOWN_INDEXES = [p + ".tbi" for p in REF["known_sites"]]


def java_options(wildcards, threads):
    # Respect the job's CPU allocation and use a project-local scratch directory.
    return (f"-Xmx{MEM['java_heap_mb']}m -XX:ActiveProcessorCount={threads} "
            f"-Djava.io.tmpdir=results/tmp/{wildcards.sample}")


wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLES)

rule all:
    input:
        "results/qc/multiqc_report.html",
        expand("results/bam/recal/{sample}.bam", sample=SAMPLES),
        expand("results/bam/recal/{sample}.bam.bai", sample=SAMPLES),
        expand("results/qc/capture/{sample}.per_target.tsv", sample=SAMPLES),
        expand("results/recal/{sample}.recal.table", sample=SAMPLES),
        expand("results/vcf/filtered/{sample}.vcf.gz", sample=TUMORS),
        expand("results/vcf/filtered/{sample}.vcf.gz.tbi", sample=TUMORS),
        expand("results/vcf/filtered/{sample}.filtering.stats", sample=TUMORS),
        expand("results/vcf/pass/{sample}.pass.vcf.gz", sample=TUMORS),
        expand("results/vcf/pass/{sample}.pass.vcf.gz.tbi", sample=TUMORS),
        expand("results/qc/variants/{sample}.bcftools.stats.txt", sample=TUMORS)

include: "all_rules.smk"
```

`rule all` defines the final deliverables. Input and output paths are relative to the project working directory; the include path is relative to the Snakefile. The execution commands below therefore specify `--snakefile workflow/Snakefile` explicitly.

## 5. Rules: `workflow/all_rules.smk`

Concatenate the Python blocks in Sections 5.1–5.10, in order, into this single file. Snakemake creates parent directories for declared outputs and logs. Scratch directories used by Java are created explicitly.

### 5.1 Adapter trimming, read filtering, and FASTQ QC

`-i/-I` specify R1/R2 and `-o/-O` specify the paired outputs. `--detect_adapter_for_pe` enables paired-end adapter autodetection. The HTML report includes before/after quality plots; JSON is retained for MultiQC. Reads failing pair retention criteria are not passed to alignment.

```python
rule fastp:
    input:
        r1="data/raw/{sample}_1.fq.gz",
        r2="data/raw/{sample}_2.fq.gz"
    output:
        r1="data/clean/{sample}_1.fq.gz",
        r2="data/clean/{sample}_2.fq.gz",
        html="results/qc/fastp/{sample}.html",
        json="results/qc/fastp/{sample}.json"
    threads: config["threads"]["fastp"]
    resources:
        mem_mb=MEM["fastp_mb"]
    params:
        quality=config["fastp"]["qualified_quality_phred"],
        unqualified=config["fastp"]["unqualified_percent_limit"],
        n_limit=config["fastp"]["n_base_limit"],
        length=config["fastp"]["length_required"]
    log: "logs/fastp/{sample}.log"
    shell:
        r"""
        fastp -i {input.r1:q} -I {input.r2:q} -o {output.r1:q} -O {output.r2:q} \
            --html {output.html:q} --json {output.json:q} --thread {threads} \
            --detect_adapter_for_pe --qualified_quality_phred {params.quality} \
            --unqualified_percent_limit {params.unqualified} \
            --n_base_limit {params.n_limit} --length_required {params.length} > {log:q} 2>&1
        """
```

### 5.2 BWA-MEM alignment and coordinate sorting

The read group records sample (`SM`), read-group ID (`ID`), library (`LB`), and sequencing platform (`PL`). The generated IDs represent the assumed single library per sample. `-M` marks shorter split alignments as secondary, retaining the original compatibility setting. BWA's fixed `-K` batch size helps keep its input batching consistent when thread counts change. This rule retains alignment records for downstream GATK read filtering rather than imposing a separate MAPQ cutoff. [BWA manual](https://bio-bwa.sourceforge.net/bwa.shtml)

```python
rule bwa_mem:
    input:
        r1="data/clean/{sample}_1.fq.gz",
        r2="data/clean/{sample}_2.fq.gz",
        index=BWA_INDEX
    output:
        bam=temp("results/bam/unsorted/{sample}.bam")
    threads: config["threads"]["mapping"]
    resources:
        mem_mb=MEM["mapping_mb"]
    params:
        index=REF["bwa_index"],
        bwa_threads=lambda wc, threads: max(1, threads - 1)
    log:
        bwa="logs/bwa/{sample}.log",
        samtools="logs/bwa/{sample}.samtools.log"
    shell:
        r"""
        bwa mem -M -K 100000000 -t {params.bwa_threads} \
            -R '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tLB:{wildcards.sample}.lib1\tPL:ILLUMINA' \
            {params.index:q} {input.r1:q} {input.r2:q} 2> {log.bwa:q} \
          | samtools view -u -o {output.bam:q} - 2> {log.samtools:q}
        """


rule sort_bam:
    input:
        bam="results/bam/unsorted/{sample}.bam"
    output:
        bam=temp("results/bam/aligned/{sample}.bam"),
        flagstat="results/qc/alignment/aligned/{sample}.flagstat.txt"
    threads: config["threads"]["samtools"]
    resources:
        mem_mb=MEM["sort_mb"]
    params:
        extra=lambda wc, threads: max(0, threads - 1)
    log: "logs/sort/{sample}.log"
    shell:
        r"""
        samtools sort -@ {params.extra} -m 1G -o {output.bam:q} {input.bam:q} 2> {log:q}
        samtools flagstat -@ {params.extra} {output.bam:q} > {output.flagstat:q} 2>> {log:q}
        """
```

`samtools sort -m 1G` is an approximate memory budget **per sorting thread**. Review the total sorting reservation if increasing the thread count. The aligned flagstat report is retained before duplicate removal.

### 5.3 Remove duplicates and index the resulting BAM

`--REMOVE_DUPLICATES true` physically excludes duplicate records, as in the source workflow. `--CREATE_INDEX false` prevents implicit index naming; SAMtools then creates the explicitly declared `.bam.bai`. The duplication metrics are retained even when intermediate BAM files are cleaned up.

```python
rule remove_duplicates:
    input:
        bam="results/bam/aligned/{sample}.bam"
    output:
        bam=temp("results/bam/dedup/{sample}.bam"),
        bai=temp("results/bam/dedup/{sample}.bam.bai"),
        metrics="results/qc/duplicates/{sample}.duplicate_metrics.txt"
    threads: config["threads"]["gatk"]
    resources:
        mem_mb=MEM["gatk_mb"]
    params:
        gatk=SW["gatk"],
        java=java_options
    log: "logs/duplicates/{sample}.log"
    shell:
        r"""
        mkdir -p results/tmp/{wildcards.sample}
        {params.gatk:q} --java-options {params.java:q} MarkDuplicates \
            --INPUT {input.bam:q} --OUTPUT {output.bam:q} --METRICS_FILE {output.metrics:q} \
            --REMOVE_DUPLICATES true --CREATE_INDEX false > {log:q} 2>&1
        samtools index {output.bam:q} {output.bai:q} 2>> {log:q}
        """
```

### 5.4 Base quality score recalibration

BaseRecalibrator models systematic quality-score errors using known variant sites. `--known-sites` masks established polymorphisms during model training; it does not label those positions as somatic. Training uses targets plus the configured padding. ApplyBQSR applies the model to the whole deduplicated BAM and writes an explicitly indexed final BAM.

```python
rule base_recalibrator:
    input:
        bam="results/bam/dedup/{sample}.bam",
        bai="results/bam/dedup/{sample}.bam.bai",
        ref=REF["fasta"],
        support=REF_SUPPORT,
        known=REF["known_sites"],
        known_indexes=KNOWN_INDEXES,
        targets=REF["targets"]
    output:
        table="results/recal/{sample}.recal.table"
    threads: config["threads"]["gatk"]
    resources:
        mem_mb=MEM["gatk_mb"]
    params:
        gatk=SW["gatk"],
        java=java_options,
        known_args=" ".join("--known-sites " + shlex.quote(p) for p in REF["known_sites"]),
        padding=config["analysis"]["interval_padding"]
    log: "logs/bqsr/{sample}.base_recalibrator.log"
    shell:
        r"""
        mkdir -p results/tmp/{wildcards.sample}
        {params.gatk:q} --java-options {params.java:q} BaseRecalibrator \
            -R {input.ref:q} -I {input.bam:q} {params.known_args} \
            -L {input.targets:q} --interval-padding {params.padding} \
            -O {output.table:q} > {log:q} 2>&1
        """


rule apply_bqsr:
    input:
        bam="results/bam/dedup/{sample}.bam",
        bai="results/bam/dedup/{sample}.bam.bai",
        ref=REF["fasta"],
        support=REF_SUPPORT,
        table="results/recal/{sample}.recal.table"
    output:
        bam="results/bam/recal/{sample}.bam",
        bai="results/bam/recal/{sample}.bam.bai"
    threads: config["threads"]["gatk"]
    resources:
        mem_mb=MEM["gatk_mb"]
    params:
        gatk=SW["gatk"],
        java=java_options
    log: "logs/bqsr/{sample}.apply_bqsr.log"
    shell:
        r"""
        mkdir -p results/tmp/{wildcards.sample}
        {params.gatk:q} --java-options {params.java:q} ApplyBQSR \
            -R {input.ref:q} -I {input.bam:q} --bqsr-recal-file {input.table:q} \
            --create-output-bam-index false -O {output.bam:q} > {log:q} 2>&1
        samtools index {output.bam:q} {output.bai:q} 2>> {log:q}
        """
```

### 5.5 Alignment and capture coverage QC

SAMtools reports describe the final, duplicate-removed BAM. Their counts differ from the aligned BAM report because they represent different processing stages. CollectHsMetrics uses the actual capture targets and baits, with explicit MAPQ/base-quality thresholds. Inspect mean target coverage, coverage breadth, and coverage uniformity alongside the duplicate metrics. Metrics labeled `PCT_*` in the raw Picard report are generally fractions. Coverage is assessed after duplicate removal; use the separate duplicate report to assess the original duplication rate. [CollectHsMetrics documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360036725371-CollectHsMetrics-Picard)

```python
rule alignment_qc:
    input:
        bam="results/bam/recal/{sample}.bam",
        bai="results/bam/recal/{sample}.bam.bai",
        ref=REF["fasta"],
        support=REF_SUPPORT
    output:
        flagstat="results/qc/alignment/recal/{sample}.flagstat.txt",
        stats="results/qc/alignment/recal/{sample}.samtools.stats.txt"
    threads: config["threads"]["samtools"]
    params:
        extra=lambda wc, threads: max(0, threads - 1)
    log: "logs/alignment_qc/{sample}.log"
    shell:
        r"""
        samtools flagstat -@ {params.extra} {input.bam:q} > {output.flagstat:q} 2> {log:q}
        samtools stats -@ {params.extra} -r {input.ref:q} {input.bam:q} > {output.stats:q} 2>> {log:q}
        """


rule capture_qc:
    input:
        bam="results/bam/recal/{sample}.bam",
        bai="results/bam/recal/{sample}.bam.bai",
        ref=REF["fasta"],
        support=REF_SUPPORT,
        targets=REF["targets"],
        baits=REF["baits"]
    output:
        metrics="results/qc/capture/{sample}.hs_metrics.txt",
        per_target="results/qc/capture/{sample}.per_target.tsv"
    threads: config["threads"]["gatk"]
    resources:
        mem_mb=MEM["gatk_mb"]
    params:
        gatk=SW["gatk"],
        java=java_options,
        mapq=config["capture_qc"]["min_mapping_quality"],
        baseq=config["capture_qc"]["min_base_quality"],
        cap=config["capture_qc"]["coverage_cap"]
    log: "logs/capture_qc/{sample}.log"
    shell:
        r"""
        mkdir -p results/tmp/{wildcards.sample}
        {params.gatk:q} --java-options {params.java:q} CollectHsMetrics \
            --INPUT {input.bam:q} --OUTPUT {output.metrics:q} --REFERENCE_SEQUENCE {input.ref:q} \
            --BAIT_INTERVALS {input.baits:q} --TARGET_INTERVALS {input.targets:q} \
            --PER_TARGET_COVERAGE {output.per_target:q} \
            --MINIMUM_MAPPING_QUALITY {params.mapq} --MINIMUM_BASE_QUALITY {params.baseq} \
            --COVERAGE_CAP {params.cap} > {log:q} 2>&1
        """
```

### 5.6 Mutect2 in tumor-only or matched tumor-normal mode

Mutect2 uses the germline population resource and PoN in both modes. When `PAIRS[sample]` names a normal, the rule also provides its recalibrated BAM and explicitly labels it with `-normal`. Omitting `-normal` while supplying two BAMs would not correctly identify the normal sample. With a `null` pairing, only the tumor BAM is supplied.

`--f1r2-tar-gz` records orientation counts; LearnReadOrientationModel learns artifact priors for filtering. This is relevant to orientation-dependent damage, including artifacts that may occur in FFPE samples. The `.vcf.gz.stats` file is tracked explicitly and passed to FilterMutectCalls. [GATK somatic calling workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)

```python
rule mutect2:
    input:
        bam="results/bam/recal/{sample}.bam",
        bai="results/bam/recal/{sample}.bam.bai",
        normal=lambda wc: [f"results/bam/recal/{PAIRS[wc.sample]}.bam"] if PAIRS[wc.sample] else [],
        normal_index=lambda wc: [f"results/bam/recal/{PAIRS[wc.sample]}.bam.bai"] if PAIRS[wc.sample] else [],
        ref=REF["fasta"],
        support=REF_SUPPORT,
        germline=REF["germline"],
        germline_index=REF["germline"] + ".tbi",
        pon=REF["pon"],
        pon_index=REF["pon"] + ".tbi",
        targets=REF["targets"]
    output:
        vcf="results/vcf/unfiltered/{sample}.vcf.gz",
        tbi="results/vcf/unfiltered/{sample}.vcf.gz.tbi",
        stats="results/vcf/unfiltered/{sample}.vcf.gz.stats",
        f1r2="results/orientation/{sample}.f1r2.tar.gz"
    threads: config["threads"]["mutect2"]
    resources:
        mem_mb=MEM["gatk_mb"]
    params:
        gatk=SW["gatk"],
        java=java_options,
        pairhmm=lambda wc, threads: max(1, threads - 1),
        normal_args=lambda wc, input: (
            "-I " + shlex.quote(input.normal[0]) + " -normal " + shlex.quote(PAIRS[wc.sample])
            if input.normal else ""
        ),
        padding=config["analysis"]["interval_padding"]
    log: "logs/mutect2/{sample}.log"
    shell:
        r"""
        mkdir -p results/tmp/{wildcards.sample}
        {params.gatk:q} --java-options {params.java:q} Mutect2 \
            -R {input.ref:q} -I {input.bam:q} {params.normal_args} \
            --germline-resource {input.germline:q} --panel-of-normals {input.pon:q} \
            -L {input.targets:q} --interval-padding {params.padding} \
            --native-pair-hmm-threads {params.pairhmm} --f1r2-tar-gz {output.f1r2:q} \
            --create-output-variant-index true -O {output.vcf:q} > {log:q} 2>&1
        """


rule learn_orientation:
    input:
        f1r2="results/orientation/{sample}.f1r2.tar.gz"
    output:
        priors="results/orientation/{sample}.artifact_priors.tar.gz"
    threads: config["threads"]["gatk"]
    resources:
        mem_mb=MEM["gatk_mb"]
    params:
        gatk=SW["gatk"],
        java=java_options
    log: "logs/orientation/{sample}.log"
    shell:
        r"""
        mkdir -p results/tmp/{wildcards.sample}
        {params.gatk:q} --java-options {params.java:q} LearnReadOrientationModel \
            -I {input.f1r2:q} -O {output.priors:q} > {log:q} 2>&1
        """
```

### 5.7 Pileup summaries and contamination estimation

GetPileupSummaries counts alleles at common SNPs overlapping capture targets. With multiple `-L` inputs, `--interval-set-rule INTERSECTION` is essential here: the default union would include sites outside the intended intersection. The same rule runs on the matched normal when one is configured. [GetPileupSummaries interval arguments](https://gatk.broadinstitute.org/hc/en-us/articles/360042913771-GetPileupSummaries)

CalculateContamination uses the tumor pileups and, in paired mode, the normal pileups supplied with `--matched-normal`. Its outputs are the tumor contamination estimate and minor-allele-fraction segments used by FilterMutectCalls. This segmentation is not a standalone copy-number analysis. [CalculateContamination matched-normal mode](https://gatk.broadinstitute.org/hc/en-us/articles/4418054253211-CalculateContamination)

```python
rule pileup_summaries:
    input:
        bam="results/bam/recal/{sample}.bam",
        bai="results/bam/recal/{sample}.bam.bai",
        ref=REF["fasta"],
        support=REF_SUPPORT,
        common=REF["common_snps"],
        common_index=REF["common_snps"] + ".tbi",
        targets=REF["targets"]
    output:
        table="results/contamination/{sample}.pileups.table"
    threads: config["threads"]["gatk"]
    resources:
        mem_mb=MEM["gatk_mb"]
    params:
        gatk=SW["gatk"],
        java=java_options
    log: "logs/contamination/{sample}.pileups.log"
    shell:
        r"""
        mkdir -p results/tmp/{wildcards.sample}
        {params.gatk:q} --java-options {params.java:q} GetPileupSummaries \
            -R {input.ref:q} -I {input.bam:q} -V {input.common:q} \
            -L {input.common:q} -L {input.targets:q} --interval-set-rule INTERSECTION \
            -O {output.table:q} > {log:q} 2>&1
        """


rule calculate_contamination:
    input:
        table="results/contamination/{sample}.pileups.table",
        normal=lambda wc: [f"results/contamination/{PAIRS[wc.sample]}.pileups.table"] if PAIRS[wc.sample] else []
    output:
        table="results/contamination/{sample}.contamination.table",
        segments="results/contamination/{sample}.segments.table"
    threads: config["threads"]["gatk"]
    resources:
        mem_mb=MEM["gatk_mb"]
    params:
        gatk=SW["gatk"],
        java=java_options,
        matched_args=lambda wc, input: "--matched-normal " + shlex.quote(input.normal[0]) if input.normal else ""
    log: "logs/contamination/{sample}.estimate.log"
    shell:
        r"""
        mkdir -p results/tmp/{wildcards.sample}
        {params.gatk:q} --java-options {params.java:q} CalculateContamination \
            -I {input.table:q} {params.matched_args} -O {output.table:q} \
            --tumor-segmentation {output.segments:q} > {log:q} 2>&1
        """
```

Inspect warnings and the number of informative SNPs. An empty or uninformative pileup table does not establish zero contamination; review interval compatibility, coverage, and the SNP resource instead of inserting a zero-valued result.

### 5.8 Filter Mutect2 calls

FilterMutectCalls combines the original calls, Mutect2 statistics, contamination estimate, segmentation, and orientation priors. It labels records in the FILTER field and retains filtered records in the output VCF for review. The next rule exports PASS records separately. [FilterMutectCalls documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360042476952-FilterMutectCalls)

```python
rule filter_mutect_calls:
    input:
        vcf="results/vcf/unfiltered/{sample}.vcf.gz",
        tbi="results/vcf/unfiltered/{sample}.vcf.gz.tbi",
        stats="results/vcf/unfiltered/{sample}.vcf.gz.stats",
        contamination="results/contamination/{sample}.contamination.table",
        segments="results/contamination/{sample}.segments.table",
        priors="results/orientation/{sample}.artifact_priors.tar.gz",
        ref=REF["fasta"],
        support=REF_SUPPORT
    output:
        vcf="results/vcf/filtered/{sample}.vcf.gz",
        tbi="results/vcf/filtered/{sample}.vcf.gz.tbi",
        stats="results/vcf/filtered/{sample}.filtering.stats"
    threads: config["threads"]["gatk"]
    resources:
        mem_mb=MEM["gatk_mb"]
    params:
        gatk=SW["gatk"],
        java=java_options
    log: "logs/filter_mutect/{sample}.log"
    shell:
        r"""
        mkdir -p results/tmp/{wildcards.sample}
        {params.gatk:q} --java-options {params.java:q} FilterMutectCalls \
            -R {input.ref:q} -V {input.vcf:q} --stats {input.stats:q} \
            --contamination-table {input.contamination:q} --tumor-segmentation {input.segments:q} \
            --ob-priors {input.priors:q} --filtering-stats {output.stats:q} \
            --create-output-variant-index true -O {output.vcf:q} > {log:q} 2>&1
        """
```

### 5.9 Export PASS records and summarize variant calls

`bcftools view --apply-filters PASS` retains records whose FILTER value is PASS, while preserving the full VCF header. `-Oz` writes BGZF-compressed VCF and `index --tbi` generates its index. No fixed positions within the FORMAT column are assumed. A valid VCF with no PASS variants may contain only headers and should not be replaced by a fabricated call. [bcftools documentation](https://www.htslib.org/doc/bcftools.html)

```python
rule export_pass:
    input:
        vcf="results/vcf/filtered/{sample}.vcf.gz",
        tbi="results/vcf/filtered/{sample}.vcf.gz.tbi"
    output:
        vcf="results/vcf/pass/{sample}.pass.vcf.gz",
        tbi="results/vcf/pass/{sample}.pass.vcf.gz.tbi",
        stats="results/qc/variants/{sample}.bcftools.stats.txt"
    log: "logs/pass/{sample}.log"
    shell:
        r"""
        bcftools view --apply-filters PASS -Oz -o {output.vcf:q} {input.vcf:q} 2> {log:q}
        bcftools index --tbi --force {output.vcf:q} 2>> {log:q}
        bcftools stats {output.vcf:q} > {output.stats:q} 2>> {log:q}
        """
```

### 5.10 Aggregate quality reports and basic plots

MultiQC waits for fastp, both alignment stages, duplicate metrics, capture metrics, and variant summaries. Directory prefixes distinguish aligned and recalibrated reports. The HTML report provides supported plots such as base-quality profiles, insert-size distributions, and target coverage summaries. Per-target coverage tables and detailed GATK filtering logs remain separate files. [MultiQC Picard module](https://docs.seqera.io/multiqc/modules/picard), [SAMtools module](https://docs.seqera.io/multiqc/modules/samtools)

```python
rule multiqc:
    input:
        fastp=expand("results/qc/fastp/{sample}.json", sample=SAMPLES),
        aligned=expand("results/qc/alignment/aligned/{sample}.flagstat.txt", sample=SAMPLES),
        recal=expand("results/qc/alignment/recal/{sample}.flagstat.txt", sample=SAMPLES),
        stats=expand("results/qc/alignment/recal/{sample}.samtools.stats.txt", sample=SAMPLES),
        duplicates=expand("results/qc/duplicates/{sample}.duplicate_metrics.txt", sample=SAMPLES),
        capture=expand("results/qc/capture/{sample}.hs_metrics.txt", sample=SAMPLES),
        variants=expand("results/qc/variants/{sample}.bcftools.stats.txt", sample=TUMORS)
    output:
        html="results/qc/multiqc_report.html",
        data=directory("results/qc/multiqc_data")
    resources:
        mem_mb=MEM["multiqc_mb"]
    log: "logs/multiqc.log"
    shell:
        r"""
        multiqc results/qc/fastp results/qc/alignment results/qc/duplicates \
            results/qc/capture results/qc/variants --dirs --dirs-depth 1 \
            --outdir results/qc --filename multiqc_report.html --force > {log:q} 2>&1
        """
```

## 6. Run the workflow and record provenance

### 6.1 Activate the environment and check tools

```bash
cd "[path_to_your_dir]/WES_project"
conda activate "[path_to_your_software_environment]"
export PATH="[path_to_your_software]/bin:$PATH"

for tool in snakemake python fastp bwa samtools bcftools multiqc java; do
    command -v "$tool" || exit 1
done
test -x "[path_to_your_software]/gatk/gatk"
"[path_to_your_software]/gatk/gatk" --version
```

The shell environment and YAML software paths must agree. If not using Conda, omit the activation and Conda export commands and retain an equivalent environment manifest or container digest.

### 6.2 Dry run and execution

```bash
# Check inputs and dependency construction; print commands without executing analysis.
snakemake --snakefile workflow/Snakefile --cores 16 --resources mem_mb=48000 \
    --dry-run --printshellcmds

# Run within the allocated CPU and memory budgets; rerun incomplete jobs after interruption.
snakemake --snakefile workflow/Snakefile --cores 16 --resources mem_mb=48000 \
    --rerun-incomplete --printshellcmds
```

For background execution, use this command **instead of** the foreground run:

```bash
nohup snakemake --snakefile workflow/Snakefile --cores 16 --resources mem_mb=48000 \
    --rerun-incomplete --printshellcmds > logs/snakemake.run.log 2>&1 &
```

`--cores` and `--resources mem_mb` limit scheduled concurrency based on rule requests; they are not operating-system enforcement of every subprocess's actual usage. `-Xmx` limits Java heap, not total process memory. Keep non-heap and operating-system headroom, and adjust requests to observed usage. On a cluster, run inside an allocated compute job or use a separately validated executor profile. The example assumes at least eight available cores and sufficient temporary disk space.

Snakemake automatically removes BAMs marked `temp()` after their consumers finish. Final recalibrated BAMs, clean FASTQ files, VCFs, metrics, and logs remain. Add `--notemp` during validation if you need to inspect intermediate BAMs. Do not add manual `rm` commands to rules.

### 6.3 Archive configuration, software versions, and input checksums

Replace placeholders before taking this snapshot. For a run with changed parameters, create a separate project directory so that earlier results and records remain available.

```bash
mkdir -p results/provenance
cp config/config.yaml results/provenance/config.used.yaml
cp workflow/Snakefile workflow/all_rules.smk results/provenance/

# Export exact versions/builds from the complete analysis environment.
conda list --explicit > results/provenance/conda-explicit.txt
snakemake --version > results/provenance/snakemake.version.txt
"[path_to_your_software]/gatk/gatk" --version > results/provenance/gatk.version.txt 2>&1
java -version > results/provenance/java.version.txt 2>&1

# Hash large input files without modifying them. This may take time.
python - <<'PY'
from pathlib import Path
import hashlib
import yaml

c = yaml.safe_load(Path("config/config.yaml").read_text())
r = c["ref"]
paths = [Path(f"data/raw/{sample}_{mate}.fq.gz")
         for sample in c["samples"] for mate in (1, 2)]
paths += [Path(r["fasta"]), Path(r["fasta"] + ".fai"), Path(r["dict"])]
paths += [Path(r["bwa_index"] + suffix) for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")]
vcfs = r["known_sites"] + [r["germline"], r["pon"], r["common_snps"]]
paths += [Path(name + suffix) for name in vcfs for suffix in ("", ".tbi")]
paths += [Path(r["targets"]), Path(r["baits"])]
with open("results/provenance/inputs.sha256", "w") as output:
    for path in dict.fromkeys(paths):
        digest = hashlib.sha256()
        with path.open("rb") as source:
            for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
                digest.update(block)
        output.write(f"{digest.hexdigest()}  {path}\n")
PY
```

Retain the matching GATK installation if it is outside the exported Conda environment. For every run, also record the reference/resource release identifiers and capture-kit version. Before rerunning, check `sha256sum -c results/provenance/inputs.sha256`. If files move, update their paths in the manifest while preserving the recorded checksums.

## 7. Outputs and completion checks

Preprocessing and QC outputs are generated for every ID in `samples`. Variant callsets, orientation models, contamination estimates, and filtering statistics are generated for tumor IDs only. Pileup tables are additionally generated for the matched normals. In paired mode, a tumor-named VCF contains both tumor and normal sample columns; identify them by the VCF header, not by fixed column positions.

| Output | Purpose |
| --- | --- |
| `data/clean/{sample}_1.fq.gz`, `_2.fq.gz` | Paired reads after fastp processing |
| `results/qc/fastp/{sample}.html`, `.json` | Per-sample before/after FASTQ QC |
| `results/qc/alignment/aligned/{sample}.flagstat.txt` | Alignment counts before duplicate removal |
| `results/qc/duplicates/{sample}.duplicate_metrics.txt` | Duplicate statistics for the original aligned library |
| `results/recal/{sample}.recal.table` | BQSR model trained on padded target regions |
| `results/bam/recal/{sample}.bam`, `.bam.bai` | Final coordinate-sorted, duplicate-removed, recalibrated BAM and index |
| `results/qc/alignment/recal/{sample}.*.txt` | Alignment statistics for the final BAM |
| `results/qc/capture/{sample}.hs_metrics.txt`, `.per_target.tsv` | Capture metrics and per-target coverage after duplicate removal |
| `results/vcf/unfiltered/{sample}.vcf.gz`, `.vcf.gz.tbi`, `.vcf.gz.stats` | Raw Mutect2 calls, index, and calling statistics |
| `results/orientation/{sample}.f1r2.tar.gz`, `.artifact_priors.tar.gz` | Read-orientation counts and learned artifact priors |
| `results/contamination/{sample}.pileups.table` | Allele counts at common SNPs within unpadded targets |
| `results/contamination/{sample}.contamination.table`, `.segments.table` | Contamination estimate and segmentation used for filtering |
| `results/vcf/filtered/{sample}.vcf.gz`, `.vcf.gz.tbi` | Full VCF with FILTER labels and its index |
| `results/vcf/filtered/{sample}.filtering.stats` | FilterMutectCalls statistics |
| `results/vcf/pass/{sample}.pass.vcf.gz`, `.pass.vcf.gz.tbi` | PASS-only candidate calls and index, including padded calling territory |
| `results/qc/variants/{sample}.bcftools.stats.txt` | Basic summary of PASS records |
| `results/qc/multiqc_report.html`, `multiqc_data/` | Combined QC report, plots, and underlying report data |
| `results/provenance/`, `logs/`, `.snakemake/` | Frozen inputs/configuration, software records, logs, and workflow state |

After completion, confirm that another dry run has no pending jobs and perform basic output checks:

```bash
snakemake --snakefile workflow/Snakefile --cores 16 --resources mem_mb=48000 --dry-run
samtools quickcheck -v results/bam/recal/*.bam

# Example checks for one sample; repeat for every configured sample.
bcftools view -h results/vcf/pass/tumor_01.pass.vcf.gz > /dev/null
bcftools index --nrecords results/vcf/pass/tumor_01.pass.vcf.gz

# In paired mode, expect both tumor_01 and normal_01; order is not assumed.
bcftools query --list-samples results/vcf/pass/tumor_01.pass.vcf.gz
```

Review fastp retention and quality plots, alignment rates, duplicate levels, insert sizes, target depth/breadth, contamination warnings, FilterMutectCalls statistics, and the number of PASS variants. Unexpectedly low coverage or zero calls warrants review of both the data and resource compatibility. An empty PASS callset alone is not proof that the sample has no somatic mutations. Set acceptance criteria for the experiment rather than treating a single depth or allele fraction threshold as universally sufficient.

> This is a reusable upstream template. Complete placeholder replacement, environment validation, a dry run, and a representative real-data run before applying it to a full cohort. A successful dependency check does not validate reference compatibility, variant accuracy, or biological sample quality.
