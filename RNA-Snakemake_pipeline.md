# RNA-seq Snakemake Upstream Analysis Pipeline Template

This template processes **paired-end bulk RNA-seq** data using **FastQC → fastp → STAR → featureCounts**, with alignment QC, library-strand checks, bigWig tracks, and a combined MultiQC report. Run the commands in **Linux with Bash**, from the project root. All paths and sample names below are examples to replace.

The workflow ends at unnormalized gene counts and basic sequencing/alignment QC. Differential expression, expression-based PCA and clustering, TPM calculations, gene-set enrichment, and project-specific biological interpretation are outside its scope.

## 1. Assumptions and analysis choices

- Provide one pair of gzip-compressed, Phred+33 FASTQ files for each non-UMI library. This template does not cover single-end, single-cell, long-read, or UMI processing. If one library spans several lanes, combine R1 and R2 in identical lane order and record the operation before starting.
- Replace `[path_to_your_...]` placeholders and example sample IDs. Use letters, digits, underscores, or hyphens for sample IDs; avoid spaces in project and software paths.
- Select the actual species and genome assembly. The FASTA, GTF, STAR index, and RSeQC BED12 annotation must have matching sequences/contigs and compatible annotation versions. Do not infer species from an example reference name.
- Set library strandedness explicitly. The supplied placeholder permits the `strand_check` target, but prevents a full counting run until replaced by the integer `0`, `1`, or `2`. All samples in one run must use the same orientation setting.
- Reads are not deduplicated. Highly expressed transcripts can produce repeated reads, while excessive duplication can also indicate limited library complexity. Interpret duplication alongside the other QC metrics.
- The source toolchain is retained, with explicit preprocessing and resource settings. Counting uses exons grouped by `gene_id`, replacing the source's gene-span and gene-symbol settings. Separate sorting, strand checks, and coverage tracks are also added. These choices can change results relative to the original notes.
- Gene ID version suffixes are preserved. The final count matrix contains integer fragment counts, including zero-count genes; it is not normalized expression.

## 2. Project, software, and reference setup

### 2.1 Project layout

```bash
mkdir -p "[path_to_your_dir]/RNA_project"
cd "[path_to_your_dir]/RNA_project"
mkdir -p config workflow/scripts data/raw data/clean results logs
```

Save the following sections into the indicated files, without the Markdown fences.

```text
RNA_project/
|-- config/config.yaml
|-- workflow/Snakefile
|-- workflow/all_rules.smk
|-- workflow/scripts/export_counts.py
|-- data/raw/
|   |-- sample_01.R1.fq.gz
|   |-- sample_01.R2.fq.gz
|   |-- sample_02.R1.fq.gz
|   `-- sample_02.R2.fq.gz
|-- data/clean/
|-- results/
`-- logs/
```

Use symbolic links with these names if necessary, preserving the source files. This naming convention also appears in the checksum command in Section 7.

### 2.2 Software requirements

| Software | Purpose |
| --- | --- |
| Snakemake 8.x, Python, PyYAML | Workflow scheduling, configuration, count export, and manifests |
| FastQC and compatible Java | Raw-read QC |
| fastp | Paired adapter trimming, read filtering, and before/after QC plots |
| STAR 2.7.x | Splice-aware alignment, transcriptome BAM, and supplementary gene counts |
| SAMtools | Coordinate sorting, BAM indexing, and alignment statistics |
| Subread featureCounts with `--countReadPairs` support | Gene-level fragment counting |
| RSeQC | Library-strand diagnostics using `infer_experiment.py` |
| deepTools | bigWig coverage generation using `bamCoverage` |
| MultiQC | Combined interactive QC plots |
| Bash, gzip, coreutils | Shell commands, decompression, and checksums |

Use one validated environment containing these tools and freeze the exact versions. The version families above are compatibility guidance, not an environment lock. Set `software.bin` to its executable directory; required supporting libraries and Java must also be available. This template uses an existing environment and does not declare per-rule Conda environments.

### 2.3 Prepare a compatible STAR index

Prefer an existing index with documented FASTA, GTF, STAR version, and indexing parameters. When building an index, use a new directory and save its command and logs with the reference. Replace the read-length placeholder by an integer: for 150-base reads, use `149`.

```bash
mkdir -p "[path_to_your_reference]/STAR_index"
"[path_to_your_software]/bin/STAR" \
    --runMode genomeGenerate \
    --runThreadN 8 \
    --genomeDir "[path_to_your_reference]/STAR_index" \
    --genomeFastaFiles "[path_to_your_reference]/genome.fa" \
    --sjdbGTFfile "[path_to_your_reference]/annotation.gtf" \
    --sjdbOverhang "[maximum_read_length_minus_1]" \
    --outFileNamePrefix "[path_to_your_reference]/STAR_index/build."
```

`--sjdbGTFfile` supplies annotated splice junctions; `--sjdbOverhang` should normally be the maximum input read length minus one. Small genomes may need a lower `--genomeSAindexNbases`; follow the STAR manual for the actual genome. Use an index compatible with the alignment executable. The workflow explicitly tracks the genome and annotation index components used by STAR. [STAR parameter definitions](https://github.com/alexdobin/STAR/blob/master/source/parametersDefault), [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

Also provide a **transcript-model BED12** file for RSeQC, generated from the same annotation release and using the same contig names. A simple three-column BED file is insufficient. Retain the conversion command or reference-provider metadata. [RSeQC documentation](https://rseqc.sourceforge.net/)

## 3. Configuration: `config/config.yaml`

```yaml
# Replace paths and sample IDs before running.
samples:
  - sample_01
  - sample_02

software:
  bin: "[path_to_your_software]/bin"

ref:
  fasta: "[path_to_your_reference]/genome.fa"
  gtf: "[path_to_your_reference]/annotation.gtf"
  star_index: "[path_to_your_reference]/STAR_index"
  bed12: "[path_to_your_reference]/annotation.bed12"

threads:
  fastqc: 2
  fastp: 4
  star: 8
  sort: 4
  samtools: 2
  featurecounts: 8
  bamcoverage: 4

# Total memory requested per job, in MB; adjust for the actual reference/data.
# STAR loads a separate index copy for every concurrently running STAR job.
memory_mb:
  fastqc: 2048
  fastp: 4096
  star: 48000
  sort: 8000
  samtools: 2048
  rseqc: 4096
  featurecounts: 8000
  bamcoverage: 8000
  multiqc: 4096

fastp:
  qualified_quality_phred: 15  # Bases below Q15 are considered unqualified.
  unqualified_percent_limit: 40  # Reject reads exceeding this percentage.
  n_base_limit: 5  # Reject reads with more than this many N bases.
  length_required: 15  # Minimum retained read length; validate for your assay.

featurecounts:
  # Replace with an integer: 0 = unstranded, 1 = forward, 2 = reverse.
  # If unknown, run strand_check first; see Section 6.
  strandedness: "[set_strandedness_to_0_1_or_2]"

coverage:
  bin_size: 10  # Width in bases of each coverage bin.
```

Memory resources constrain scheduling when supplied through `--resources`; they do not enforce operating-system memory limits. The example STAR request is a starting point for a large genome. Measure actual usage and adjust it before scaling up.

## 4. Workflow entry point: `workflow/Snakefile`

```python
import os
import re
from pathlib import Path

configfile: "config/config.yaml"

shell.executable("/bin/bash")
os.environ["PATH"] = config["software"]["bin"] + os.pathsep + os.environ["PATH"]

SAMPLES = config["samples"]
if not SAMPLES or len(SAMPLES) != len(set(SAMPLES)):
    raise ValueError("Provide a nonempty list of unique sample IDs.")
if any(not isinstance(s, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", s) for s in SAMPLES):
    raise ValueError("Sample IDs may contain only letters, digits, underscores, and hyphens.")

STAR_INDEX_NAMES = [
    "Genome", "SA", "SAindex", "genomeParameters.txt",
    "chrName.txt", "chrLength.txt", "chrNameLength.txt", "chrStart.txt",
    "geneInfo.tab", "transcriptInfo.tab", "exonInfo.tab", "exonGeTrInfo.tab",
    "sjdbInfo.txt", "sjdbList.out.tab",
]
STAR_INDEX = [str(Path(config["ref"]["star_index"]) / name) for name in STAR_INDEX_NAMES]
BAMS = expand("results/bam/{sample}.bam", sample=SAMPLES)
BAIS = expand("results/bam/{sample}.bam.bai", sample=SAMPLES)


def strand_setting(wildcards):
    # Deferred validation lets strand_check run before orientation is known.
    value = config["featurecounts"]["strandedness"]
    if type(value) is not int or value not in (0, 1, 2):
        raise ValueError(
            "Set featurecounts.strandedness to integer 0, 1, or 2. "
            "If unknown, run the strand_check target and inspect its reports first."
        )
    return value


wildcard_constraints:
    sample="[A-Za-z0-9_-]+",
    mate="[12]"


rule all:
    input:
        BAMS,
        BAIS,
        expand("results/star/{sample}/Aligned.toTranscriptome.out.bam", sample=SAMPLES),
        expand("results/star/{sample}/ReadsPerGene.out.tab", sample=SAMPLES),
        expand("results/star/{sample}/SJ.out.tab", sample=SAMPLES),
        expand("results/bigwig/{sample}.CPM.bw", sample=SAMPLES),
        "results/counts/gene_counts.tsv",
        "results/counts/gene_annotation.tsv",
        "results/multiqc/multiqc_report.html",
        "results/multiqc/multiqc_data"


rule strand_check:
    input:
        expand("results/qc/strandedness/{sample}.txt", sample=SAMPLES)


include: "all_rules.smk"
```

## 5. Analysis rules: `workflow/all_rules.smk`

Save **all** rules in this section to one file. Snakemake creates output dependencies and runs independent samples concurrently within the available resources.

### 5.1 Raw FASTQ QC

FastQC reports sequence quality, adapters, GC distribution, and duplication. Each mate has an explicit report. Processed-read QC is supplied by fastp in the next step.

```python
rule fastqc_raw:
    input:
        "data/raw/{sample}.R{mate}.fq.gz"
    output:
        html="results/qc/fastqc_raw/{sample}.R{mate}_fastqc.html",
        zip="results/qc/fastqc_raw/{sample}.R{mate}_fastqc.zip"
    threads: config["threads"]["fastqc"]
    resources:
        mem_mb=config["memory_mb"]["fastqc"]
    log:
        "logs/fastqc/{sample}.R{mate}.log"
    shell:
        r"""
        mkdir -p results/qc/fastqc_raw logs/fastqc
        fastqc --threads {threads} --outdir results/qc/fastqc_raw \
            {input:q} > {log:q} 2>&1
        """
```

### 5.2 Paired adapter trimming and quality filtering

`--detect_adapter_for_pe` enables paired-end adapter detection. The configured quality, N-base, and length limits determine read rejection. Only retained pairs enter STAR. JSON and HTML reports include before/after metrics and basic plots. No read merging, UMI extraction, or deduplication is enabled. [fastp documentation](https://github.com/OpenGene/fastp)

```python
rule fastp:
    input:
        r1="data/raw/{sample}.R1.fq.gz",
        r2="data/raw/{sample}.R2.fq.gz"
    output:
        r1="data/clean/{sample}.R1.fq.gz",
        r2="data/clean/{sample}.R2.fq.gz",
        html="results/qc/fastp/{sample}.html",
        json="results/qc/fastp/{sample}.json"
    threads: config["threads"]["fastp"]
    resources:
        mem_mb=config["memory_mb"]["fastp"]
    params:
        quality=config["fastp"]["qualified_quality_phred"],
        percent=config["fastp"]["unqualified_percent_limit"],
        ns=config["fastp"]["n_base_limit"],
        length=config["fastp"]["length_required"]
    log:
        "logs/fastp/{sample}.log"
    shell:
        r"""
        mkdir -p data/clean results/qc/fastp logs/fastp
        fastp --in1 {input.r1:q} --in2 {input.r2:q} \
            --out1 {output.r1:q} --out2 {output.r2:q} \
            --thread {threads} --detect_adapter_for_pe \
            --qualified_quality_phred {params.quality} \
            --unqualified_percent_limit {params.percent} \
            --n_base_limit {params.ns} --length_required {params.length} \
            --html {output.html:q} --json {output.json:q} > {log:q} 2>&1
        """
```

### 5.3 STAR splice-aware alignment

One-pass alignment uses the annotated index. `--quantMode` retains transcriptome alignments and STAR gene-count diagnostics. Unique alignments receive MAPQ 255; reads with more than ten mapping loci have no reported alignment. `--outSAMunmapped Within` retains unmapped records for QC. `NH` supports mapping-ambiguity checks; `NM` and `MD` describe alignment differences. The genome BAM is sorted separately. [STAR parameter definitions](https://github.com/alexdobin/STAR/blob/master/source/parametersDefault)

```python
rule star:
    input:
        r1="data/clean/{sample}.R1.fq.gz",
        r2="data/clean/{sample}.R2.fq.gz",
        index=STAR_INDEX,
        fasta=config["ref"]["fasta"],
        gtf=config["ref"]["gtf"]
    output:
        bam=temp("results/star/{sample}/Aligned.out.bam"),
        transcriptome="results/star/{sample}/Aligned.toTranscriptome.out.bam",
        counts="results/star/{sample}/ReadsPerGene.out.tab",
        junctions="results/star/{sample}/SJ.out.tab",
        final_log="results/star/{sample}/Log.final.out",
        detailed_log="results/star/{sample}/Log.out",
        progress_log="results/star/{sample}/Log.progress.out"
    threads: config["threads"]["star"]
    resources:
        mem_mb=config["memory_mb"]["star"]
    params:
        index=config["ref"]["star_index"],
        prefix="results/star/{sample}/"
    log:
        "logs/star/{sample}.log"
    shell:
        r"""
        mkdir -p {params.prefix:q} logs/star
        STAR --runThreadN {threads} --genomeDir {params.index:q} \
            --genomeLoad NoSharedMemory --runRNGseed 777 \
            --readFilesIn {input.r1:q} {input.r2:q} --readFilesCommand gunzip -c \
            --outFileNamePrefix {params.prefix:q} --twopassMode None \
            --outSAMtype BAM Unsorted --outSAMunmapped Within \
            --outSAMattributes NH HI AS nM NM MD \
            --outSAMmapqUnique 255 --outFilterMultimapNmax 10 \
            --quantMode TranscriptomeSAM GeneCounts > {log:q} 2>&1
        """
```

The FASTA and GTF are tracked for provenance, but STAR reads their prebuilt index during alignment. Changing them requires rebuilding the index in a new location and updating the configuration; dependency tracking cannot prove index/reference consistency. The transcriptome BAM is retained in STAR's output order and is not a genome-coordinate BAM.

### 5.4 Coordinate sorting, indexing, and alignment QC

The SAMtools rules reserve one thread for the main process when setting `-@`. Sorting allocates approximately 75% of the job's memory across its reserved threads, leaving overhead. `quickcheck` checks basic BAM integrity before indexing; `flagstat` and `stats` provide QC summaries. [SAMtools sorting documentation](https://www.htslib.org/doc/samtools-sort.html)

```python
rule sort_bam:
    input:
        "results/star/{sample}/Aligned.out.bam"
    output:
        "results/bam/{sample}.bam"
    threads: config["threads"]["sort"]
    resources:
        mem_mb=config["memory_mb"]["sort"]
    params:
        extra=lambda wildcards, threads: max(0, threads - 1),
        per_thread=lambda wildcards, threads, resources: max(1, int(resources.mem_mb * 0.75 / threads))
    log:
        "logs/samtools/{sample}.sort.log"
    shell:
        r"""
        mkdir -p results/bam logs/samtools
        samtools sort -@ {params.extra} -m {params.per_thread}M \
            -T results/bam/{wildcards.sample}.sort_tmp \
            -o {output:q} {input:q} > {log:q} 2>&1
        """


rule index_bam:
    input:
        "results/bam/{sample}.bam"
    output:
        "results/bam/{sample}.bam.bai"
    threads: config["threads"]["samtools"]
    resources:
        mem_mb=config["memory_mb"]["samtools"]
    params:
        extra=lambda wildcards, threads: max(0, threads - 1)
    log:
        "logs/samtools/{sample}.index.log"
    shell:
        r"""
        mkdir -p logs/samtools
        samtools quickcheck -v {input:q} > {log:q} 2>&1
        samtools index -@ {params.extra} {input:q} {output:q} >> {log:q} 2>&1
        """


rule alignment_qc:
    input:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bam.bai"
    output:
        flagstat="results/qc/samtools/{sample}.flagstat.txt",
        stats="results/qc/samtools/{sample}.stats.txt"
    threads: config["threads"]["samtools"]
    resources:
        mem_mb=config["memory_mb"]["samtools"]
    params:
        extra=lambda wildcards, threads: max(0, threads - 1)
    log:
        "logs/samtools/{sample}.qc.log"
    shell:
        r"""
        mkdir -p results/qc/samtools logs/samtools
        samtools flagstat -@ {params.extra} {input.bam:q} > {output.flagstat:q} 2> {log:q}
        samtools stats -@ {params.extra} {input.bam:q} > {output.stats:q} 2>> {log:q}
        """
```

BAI indexing assumes reference contigs shorter than its coordinate limit (approximately 512 Mb). A reference with larger contigs needs CSI indexing and corresponding dependency paths throughout this template. [SAMtools index documentation](https://www.htslib.org/doc/samtools-index.html)

### 5.5 Library-strand diagnostics

RSeQC samples up to 200,000 alignments and compares their orientations with transcript models. `-q 255` selects STAR's unique-MAPQ alignments under this workflow's settings. The reports inform the counting choice; they do not automatically change it. [RSeQC documentation](https://rseqc.sourceforge.net/)

```python
rule infer_strand:
    input:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bam.bai",
        bed=config["ref"]["bed12"]
    output:
        "results/qc/strandedness/{sample}.txt"
    resources:
        mem_mb=config["memory_mb"]["rseqc"]
    log:
        "logs/rseqc/{sample}.log"
    shell:
        r"""
        mkdir -p results/qc/strandedness logs/rseqc
        infer_experiment.py -i {input.bam:q} -r {input.bed:q} \
            -s 200000 -q 255 > {output:q} 2> {log:q}
        """
```

### 5.6 Gene-level fragment counting

`-p --countReadPairs` counts pairs. `-t exon -g gene_id` groups exons by stable gene ID. `-B -C` requires both mates to align and excludes chimeric pairs. `-s` controls strand interpretation. Multimapping and multi-gene overlaps are excluded by default; duplicate removal and fractional assignment are not enabled. BAM order follows the configured sample list. [Subread user guide](https://subread.sourceforge.net/SubreadUsersGuide.pdf)

```python
rule featurecounts:
    input:
        bams=BAMS,
        gtf=config["ref"]["gtf"]
    output:
        table="results/counts/featurecounts.txt",
        summary="results/counts/featurecounts.txt.summary"
    threads: config["threads"]["featurecounts"]
    resources:
        mem_mb=config["memory_mb"]["featurecounts"]
    params:
        strand=strand_setting
    log:
        "logs/featurecounts/counts.log"
    shell:
        r"""
        mkdir -p results/counts logs/featurecounts
        featureCounts -T {threads} -p --countReadPairs -B -C \
            -s {params.strand} -t exon -g gene_id \
            -a {input.gtf:q} -o {output.table:q} {input.bams:q} > {log:q} 2>&1
        """


rule export_counts:
    input:
        "results/counts/featurecounts.txt"
    output:
        counts="results/counts/gene_counts.tsv",
        annotation="results/counts/gene_annotation.tsv"
    params:
        samples=SAMPLES,
        bams=BAMS
    resources:
        mem_mb=1024
    script:
        "scripts/export_counts.py"
```

Save this helper as **`workflow/scripts/export_counts.py`**. It verifies column identity before assigning sample names, preserves annotation fields, and rejects duplicate IDs or non-integer counts.

```python
import csv
from pathlib import Path


def export_counts(source, counts_path, annotation_path, samples, bams):
    expected_annotation = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    with open(source, encoding="utf-8") as source_handle:
        reader = csv.reader((line for line in source_handle if not line.startswith("#")), delimiter="\t")
        header = next(reader)
        if header[:6] != expected_annotation:
            raise ValueError("Unexpected featureCounts annotation columns.")
        if [str(Path(p).resolve()) for p in header[6:]] != [str(Path(p).resolve()) for p in bams]:
            raise ValueError("featureCounts BAM columns do not match the configured sample order.")
        if len(samples) != len(bams) or len(set(samples)) != len(samples):
            raise ValueError("Sample names must uniquely match the BAM columns.")
        with open(counts_path, "w", newline="", encoding="utf-8") as count_handle, \
                open(annotation_path, "w", newline="", encoding="utf-8") as annotation_handle:
            counts = csv.writer(count_handle, delimiter="\t", lineterminator="\n")
            annotation = csv.writer(annotation_handle, delimiter="\t", lineterminator="\n")
            counts.writerow(["gene_id"] + list(samples))
            annotation.writerow(["gene_id"] + header[1:6])
            seen = set()
            for row in reader:
                if len(row) != len(header) or not row[0] or row[0] in seen:
                    raise ValueError("Malformed row or duplicate gene ID in featureCounts output.")
                values = row[6:]
                if any(not value.isascii() or not value.isdigit() for value in values):
                    raise ValueError("Expected nonnegative integer fragment counts.")
                seen.add(row[0])
                counts.writerow([row[0]] + values)
                annotation.writerow(row[:6])


if "snakemake" in globals():
    export_counts(
        snakemake.input[0], snakemake.output.counts, snakemake.output.annotation,
        snakemake.params.samples, snakemake.params.bams,
    )
```

### 5.7 bigWig coverage tracks

The tracks contain **combined-strand aligned-read coverage**, normalized to CPM. MAPQ 255 selects unique STAR alignments; flag mask 2308 excludes unmapped, secondary, and supplementary records. Reads are not extended across introns. Both mates can contribute coverage, so track values are not the fragment counts in the gene matrix. These tracks support genome-browser inspection and are not differential-expression inputs. [deepTools bamCoverage documentation](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)

```python
rule bamcoverage:
    input:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bam.bai"
    output:
        "results/bigwig/{sample}.CPM.bw"
    threads: config["threads"]["bamcoverage"]
    resources:
        mem_mb=config["memory_mb"]["bamcoverage"]
    params:
        bins=config["coverage"]["bin_size"]
    log:
        "logs/bigwig/{sample}.log"
    shell:
        r"""
        mkdir -p results/bigwig logs/bigwig
        bamCoverage --bam {input.bam:q} --outFileName {output:q} \
            --numberOfProcessors {threads} --binSize {params.bins} \
            --normalizeUsing CPM --minMappingQuality 255 \
            --samFlagExclude 2308 --exactScaling > {log:q} 2>&1
        """
```

`--exactScaling` computes the scaling factor from all reads instead of sampling, at additional runtime cost. CPM normalization is fixed in this rule to match the output filename.

### 5.8 Combined QC report and basic plots

MultiQC runs after all declared QC outputs exist. Its HTML contains available raw quality, trimming, alignment, counting-assignment, and strand-diagnostic summaries. Individual reports remain available if a pinned MultiQC version does not parse a particular module. [MultiQC featureCounts module](https://docs.seqera.io/multiqc/modules/featurecounts), [MultiQC RSeQC module](https://docs.seqera.io/multiqc/modules/rseqc)

```python
rule multiqc:
    input:
        fastqc=expand("results/qc/fastqc_raw/{sample}.R{mate}_fastqc.zip", sample=SAMPLES, mate=[1, 2]),
        fastp=expand("results/qc/fastp/{sample}.json", sample=SAMPLES),
        star=expand("results/star/{sample}/Log.final.out", sample=SAMPLES),
        flagstat=expand("results/qc/samtools/{sample}.flagstat.txt", sample=SAMPLES),
        stats=expand("results/qc/samtools/{sample}.stats.txt", sample=SAMPLES),
        strand=expand("results/qc/strandedness/{sample}.txt", sample=SAMPLES),
        counts="results/counts/featurecounts.txt.summary"
    output:
        html="results/multiqc/multiqc_report.html",
        data=directory("results/multiqc/multiqc_data")
    resources:
        mem_mb=config["memory_mb"]["multiqc"]
    log:
        "logs/multiqc.log"
    shell:
        r"""
        mkdir -p results/multiqc logs
        multiqc {input:q} --force --outdir results/multiqc \
            --filename multiqc_report.html > {log:q} 2>&1
        """
```

## 6. Determine library strandedness

If the protocol specifies orientation, set the integer in `config.yaml` and confirm it against the reports. If orientation is unknown, leave the placeholder temporarily and run:

```bash
export PATH="[path_to_your_software]/bin:$PATH"

# Preview only the preprocessing, mapping, and strand-diagnostic jobs.
snakemake --snakefile workflow/Snakefile --cores 16 \
    --resources mem_mb=96000 --dry-run --printshellcmds strand_check

# Produce strand reports without requiring a counting orientation.
snakemake --snakefile workflow/Snakefile --cores 16 \
    --resources mem_mb=96000 --printshellcmds --rerun-incomplete strand_check
```

Read `results/qc/strandedness/*.txt` and compare the two paired-end categories below with the library protocol.

| Observed informative alignments | Interpretation | featureCounts setting |
| --- | --- | --- |
| Similar fractions in both categories | Consistent with an unstranded library | `0` |
| Dominant `1++,1--,2+-,2-+` | Read 1 follows transcript strand | `1` |
| Dominant `1+-,1-+,2++,2--` | Read 1 opposes transcript strand | `2` |

RSeQC also reports alignments whose orientation cannot be determined. Low informative counts, mismatched annotation, contamination, or intermediate fractions need investigation; this template sets no automatic acceptance threshold. Confirm all libraries agree before choosing the shared setting. Record the protocol and diagnostic evidence in the run notes. [RSeQC strand interpretation](https://rseqc.sourceforge.net/)

## 7. Run and record the workflow

### 7.1 Record software and input identities

Keep a frozen environment, the workflow files, configuration, reference-build records, and checksums. Use a new run directory when changing inputs or major analysis settings. The following commands capture the current run's provenance; keep these local records separate from the generic public template.

```bash
mkdir -p results/provenance
export PATH="[path_to_your_software]/bin:$PATH"

{
    date -u
    snakemake --version
    python --version
    fastqc --version
    fastp --version
    STAR --version
    samtools --version
    featureCounts -v
    bamCoverage --version
    multiqc --version
    java -version
    python -m pip freeze
} > results/provenance/software_versions.txt 2>&1

# If this is a Conda environment, preserve its package builds as well.
if command -v conda >/dev/null 2>&1; then
    conda list --explicit > results/provenance/conda-explicit.txt
fi

cp config/config.yaml results/provenance/config.used.yaml
sha256sum config/config.yaml workflow/Snakefile workflow/all_rules.smk \
    workflow/scripts/export_counts.py > results/provenance/workflow.sha256

# Checksumming large FASTQ and index files can take substantial time.
python - <<'PY'
import hashlib
from pathlib import Path
import yaml

with open("config/config.yaml", encoding="utf-8") as handle:
    cfg = yaml.safe_load(handle)
paths = [Path(f"data/raw/{sample}.R{mate}.fq.gz")
         for sample in cfg["samples"] for mate in (1, 2)]
paths += [Path(cfg["ref"][key]) for key in ("fasta", "gtf", "bed12")]
paths += sorted(p for p in Path(cfg["ref"]["star_index"]).rglob("*") if p.is_file())
with open("results/provenance/inputs.sha256", "w", encoding="utf-8") as manifest:
    for path in paths:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        manifest.write(f"{digest.hexdigest()}  {path}\n")
PY
```

Python package versions include RSeQC and deepTools when installed in this environment. Review the version log for missing commands; version recording is not an installation check. Preserve relevant non-Python package builds in the environment lock as well.

### 7.2 Dry run, execution, and HTML workflow report

After replacing the strand placeholder with an integer, run:

```bash
# Resolve inputs and preview commands without executing analysis tools.
snakemake --snakefile workflow/Snakefile --cores 16 \
    --resources mem_mb=96000 --dry-run --printshellcmds

# Execute locally; pipefail preserves a failed workflow's exit status through tee.
set -o pipefail
snakemake --snakefile workflow/Snakefile --cores 16 \
    --resources mem_mb=96000 --printshellcmds --rerun-incomplete \
    2>&1 | tee results/provenance/snakemake.log

# Save the execution report after a successful run.
snakemake --snakefile workflow/Snakefile \
    --report results/provenance/snakemake_report.html

sha256sum results/counts/gene_counts.tsv results/counts/gene_annotation.tsv \
    results/counts/featurecounts.txt results/counts/featurecounts.txt.summary \
    > results/provenance/count_outputs.sha256
```

Snakemake removes the temporary unsorted genome BAM only after its consumers finish. Clean FASTQ, coordinate-sorted BAMs, transcriptome BAMs, and final reports remain. Reserve enough disk space for these retained files and STAR/sorting temporary files. When resuming, use the same command and configuration; inspect failed-job logs before rerunning.

## 8. Outputs and basic review

| Output | Purpose |
| --- | --- |
| `results/qc/fastqc_raw/` | Per-mate raw-read QC |
| `results/qc/fastp/` | Before/after filtering metrics and plots |
| `data/clean/*.fq.gz` | Retained paired reads used for alignment |
| `results/bam/*.bam` and `*.bam.bai` | Sorted genome alignments and indexes |
| `results/star/*/Log.final.out` | STAR alignment summary |
| `results/star/*/SJ.out.tab` | Detected splice junctions |
| `results/star/*/Aligned.toTranscriptome.out.bam` | Transcriptome-coordinate alignments |
| `results/star/*/ReadsPerGene.out.tab` | Supplementary STAR counting diagnostics |
| `results/qc/samtools/` | Alignment counts and statistics |
| `results/qc/strandedness/` | Evidence for library orientation |
| `results/counts/featurecounts.txt` and `.summary` | Original counting table and assignment QC |
| `results/counts/gene_counts.tsv` | Gene-by-sample integer fragment matrix |
| `results/counts/gene_annotation.tsv` | Gene identifiers and featureCounts annotation fields |
| `results/bigwig/*.CPM.bw` | Combined-strand coverage tracks for browser inspection |
| `results/multiqc/multiqc_report.html` | Combined basic QC plots |
| `results/provenance/` and `logs/` | Versions, parameters, checksums, reports, and execution logs |

Before handing off the count matrix:

1. Confirm all expected samples and both mates are present; compare raw input identities with the manifest.
2. Inspect quality, adapters, retained read pairs, and read-length changes in FastQC/fastp.
3. Review STAR unique/multiple mapping, unmapped categories, and splice-junction metrics. Investigate outliers using the actual species, library preparation, and annotation.
4. Check orientation evidence against the protocol and the saved counting setting.
5. Review featureCounts assignment categories and verify matrix columns match `samples` in `config.yaml`. Count allocation depends on this template's exon, pairing, and ambiguity policies; STAR counts need not match exactly.
6. Inspect a few representative loci with the BAM and bigWig in a browser configured for the same assembly. Coverage gaps over introns are expected for spliced alignments.
7. Preserve unfiltered integer counts and all QC records. Define project-specific acceptance criteria before biological analysis; this template does not impose universal QC cutoffs.

## 9. Validation boundary

This document is a reusable workflow template, not evidence that a sequencing run has passed QC. Replace placeholders, freeze the environment, verify reference compatibility, and run a small representative library before applying it to a complete dataset. A Snakemake dry run checks dependency construction and command formatting; it cannot verify FASTQ contents, annotation validity, installed tool behavior, or biological quality.
