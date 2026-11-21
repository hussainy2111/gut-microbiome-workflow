# Gut Microbiome Analysis Workflow

## Project overview

This project implements a reproducible analysis workflow for paired-end 16S rRNA gene sequencing data from a longitudinal mouse gut microbiome study. The workflow starts with raw FASTQ files and performs sequence assembly, quality filtering, reference alignment, chimera removal and taxonomic classification. It then generates genus-level abundance profiles and compares microbial diversity between early and late sampling periods. Mothur is used for sequence processing and taxonomic classification, Python for downstream analysis, and Nextflow to coordinate the workflow. Docker provides the software environment.

## 1. Dataset

The sequencing data are publicly available from [Zenodo](https://zenodo.org/records/800651).

The dataset contains ten biological samples collected from one mouse:

| Sampling period | Samples | Sampling days |
|---|---|---|
| Early | F3D0, F3D1, F3D2, F3D3, F3D5 | 0–5 |
| Late | F3D141, F3D142, F3D143, F3D144, F3D145 | 141–145 |

A mock community was also included in the original dataset. It was processed during the workflow but excluded from the final biological comparison. The input data consist of paired-end FASTQ files. The reference files include a SILVA alignment reference and the Mothur training-set FASTA and taxonomy files. Sample information is stored in `data/metadata.csv`, which links each sample identifier to its sampling day and sampling period.

## 2. Repository structure

```text
gut-microbiome-workflow/
├── data/
│   └── metadata.csv
├── scripts/
│   ├── download_data.sh
│   ├── stability.batch
│   └── diversity.py
├── results/
│   ├── mothur/
│   │   ├── genus.shared
│   │   ├── genus.taxonomy
│   │   └── qc/
│   └── diversity/
│       ├── genus_abundance.tsv
│       ├── genus_relative_abundance.tsv
│       ├── genus_composition.png
│       ├── phylotype_abundance.tsv
│       ├── shannon_diversity.tsv
│       └── shannon_comparison.png
├── pipeline_info/
├── Dockerfile
├── environment.yml
├── main.nf
├── nextflow.config
└── README.md
```

The `data/raw/` and `data/reference/` directories are created when the input data are downloaded. These directories are excluded from Git because the sequencing reads and reference databases are large files.

## 3. Reproducing the analysis

The analysis is organized so that another user can obtain the repository, retrieve the input data, build the software environment and execute the complete workflow.

### Input data retrieval

The script `scripts/download_data.sh` downloads the paired-end FASTQ files for all ten biological samples and the mock community from the Zenodo dataset. It also retrieves the reference files required by Mothur:

- `silva.v4.fasta`
- `trainset9_032012.pds.fasta`
- `trainset9_032012.pds.tax`

The script checks whether each file already exists before downloading it. This allows the download procedure to be repeated without downloading existing files again.

### Workflow execution

The main workflow is defined in `main.nf`. It contains two connected Nextflow processes.

**Process 1: Mothur sequence processing**

The first process receives the raw reads, reference files and Mothur batch script. It copies the FASTQ files into its task working directory and executes the commands in `scripts/stability.batch`. The process then excludes unwanted lineages and the mock community before generating the phylotype abundance table and taxonomy assignments. Mothur logs and quality-control summaries are retained as outputs.

**Process 2: Diversity and composition analysis**

The second process receives the Mothur abundance and taxonomy outputs, the sample metadata and `scripts/diversity.py`.

### Re-running analysis

The following commands describe how to reproduce the analysis from the repository.

To obtain the repository and enter its directory:

```bash
git clone https://github.com/hussainy2111/gut-microbiome-workflow.git
cd gut-microbiome-workflow
```

In order to retrieve the sequencing reads and reference databases:

```bash
bash scripts/download_data.sh
```

To build the Docker image using the supplied environment definition:

```bash
docker build -t gut-microbiome:1.0 .
```

To run the complete Nextflow workflow:

```bash
nextflow run main.nf
```

The workflow runs the Mothur and Python processes in sequence and publishes the results under `results/`. Nextflow execution information is written to `pipeline_info/`.


If the Mothur outputs have already been generated, the Python analysis can also be run independently inside the Docker image:

```bash
docker run --rm \
  -v "$PWD:/workspace" \
  -w /workspace \
  gut-microbiome:1.0 \
  python scripts/diversity.py \
    --shared results/mothur/genus.shared \
    --taxonomy results/mothur/genus.taxonomy \
    --metadata data/metadata.csv \
    --outdir results/diversity
```


## 4. Quality-control and processing decisions

The Mothur commands are maintained in `scripts/stability.batch`. The workflow follows the general Mothur MiSeq SOP approach.

| Processing step | Method | Reason |
|---|---|---|
| Paired-read assembly | `make.contigs` | Joins forward and reverse reads into contigs. |
| Initial sequence filtering | `screen.seqs` | Removes sequences with ambiguous bases, excessive length or long homopolymer runs. |
| Dereplication | `unique.seqs` | Groups identical sequences to reduce redundancy. |
| Reference alignment | `align.seqs` | Aligns sequences against the supplied SILVA reference. |
| Alignment screening | `screen.seqs` | Retains sequences within the selected alignment region. |
| Alignment filtering | `filter.seqs` | Removes uninformative alignment columns. |
| Further dereplication | `unique.seqs` | Reduces redundant sequences after alignment filtering. |
| Pre-clustering | `pre.cluster` | Groups closely related sequences to reduce variation associated with likely sequencing errors. |
| Chimera detection | `chimera.vsearch` | Identifies and removes likely chimeric sequences. |
| Taxonomic classification | `classify.seqs` | Assigns taxonomy using the supplied training reference and a confidence cutoff of 80. |
| Lineage filtering | `remove.lineage` | Excludes chloroplast, mitochondrial, unknown, archaeal and eukaryotic lineages. |
| Mock community exclusion | `remove.groups` | Excludes the mock community from the biological comparison. |

The workflow retains Mothur logs, contig reports, sequence summaries and taxonomic summaries in `results/mothur/qc/`. For the diversity comparison, I used phylotype-level counts rather than genus-collapsed counts. All ten biological samples were rarefied to 2,404 reads, the lowest sequencing depth among the samples, to compare Shannon diversity at a common read depth.

## 5. Output files

### Mothur outputs: `results/mothur/`

| File | Description |
|---|---|
| `genus.shared` | Sample-by-phylotype abundance table. |
| `genus.taxonomy` | Taxonomic assignments for the phylotypes. |
| `qc/run.logfile` | Mothur execution log. |
| `qc/` sequence summaries | Intermediate sequence and taxonomy summaries. |
| `qc/` contig report | Information about paired-read assembly. |

The filenames `genus.shared` and `genus.taxonomy` are the names used for the published workflow outputs. The shared file contains phylotype-level counts.

### Diversity outputs: `results/diversity/`

| File | Description |
|---|---|
| `phylotype_abundance.tsv` | Phylotype-level abundance counts. |
| `genus_abundance.tsv` | Genus-level read counts. |
| `genus_relative_abundance.tsv` | Genus-level relative abundances for each sample. |
| `genus_composition.png` | Plot of genus-level composition across samples. |
| `shannon_diversity.tsv` | Sample-level Shannon estimates, sequencing depth and metadata. |
| `shannon_comparison.png` | Comparison of Shannon diversity between early and late sampling periods. |

### Nextflow execution information: `pipeline_info/`

- `timeline.html`: process execution timeline.
- `report.html`: workflow execution report.
- `trace.txt`: process-level execution information.
- `flowchart.html`: workflow diagram.

## 6. Findings and interpretation

### Genus-level composition

The workflow generated genus-level abundance and relative abundance profiles for all ten biological samples. Approximately 69–79% of reads per sample could not be assigned to a genus. These reads were retained under an `Unclassified` category rather than being assigned to taxa without sufficient classification support. This limits the taxonomic detail that can be interpreted from the genus-level profiles. The unclassified fraction should be considered when interpreting the relative abundance of the identified genera.

### Phylotype-level Shannon diversity

Shannon diversity was calculated using the phylotype abundance table. The samples were rarefied to 2,404 reads before calculating the diversity estimates.

| Sampling period | Number of samples | Mean Shannon diversity | Standard deviation |
|---|---:|---:|---:|
| Early (days 0–5) | 5 | 1.844 | 0.087 |
| Late (days 141–145) | 5 | 1.686 | 0.100 |

The late-period mean was approximately 8.6% lower than the early-period mean. There was some overlap between individual samples. For example, F3D141 had a higher Shannon estimate than F3D2 and F3D3, although its estimate was lower than those of F3D0, F3D1 and F3D5. The results therefore show a lower average Shannon diversity in the late-period samples, rather than complete separation between the two sampling periods. 

## 7. Limitations and what I would not conclude

All ten samples came from the same mouse, collected at different times. They are repeated samples from one animal, not samples from ten different mice. The early and late groups therefore do not represent two separate groups of animals. Many reads could not be identified at the genus level, so the results do not give a complete picture of the bacteria present. Rarefaction allowed me to compare samples at the same sequencing depth, but it meant removing some reads from samples with more reads.

Based on these results, I would not say that gut microbial diversity always decreases with age, that the same pattern would be seen in other mice, or that these changes affect the mouse's health. The results only show that the average Shannon diversity was lower in the later samples from this mouse. Samples from more mice would be needed to see whether the same pattern occurs in other animals.


