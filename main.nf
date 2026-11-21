nextflow.enable.dsl = 2

params.raw_dir = "${projectDir}/data/raw"
params.reference_dir = "${projectDir}/data/reference"
params.metadata = "${projectDir}/data/metadata.csv"
params.outdir = "${projectDir}/results"

process MOTHUR {
    tag "16S read processing"

    container "gut-microbiome:1.0"
    cpus 2

    publishDir "${params.outdir}/mothur", mode: "copy"

    input:
    path raw_dir
    path reference_dir
    path batch_file

    output:
    path "genus.shared", emit: shared
    path "genus.taxonomy", emit: taxonomy
    path "qc", emit: qc

    script:
    """
    set -euo pipefail

    # Copy raw reads into the isolated Nextflow task directory
    mkdir -p reads
    cp "${raw_dir}"/*.fastq reads/

    export DATA="\$(realpath reads)"
    export REFERENCE_LOCATION="\$(realpath "${reference_dir}")"
    export TYPE="fastq"
    export PROC="${task.cpus}"
    export ALIGNREF="silva.v4.fasta"
    export TAXONREF_FASTA="trainset9_032012.pds.fasta"
    export TAXONREF_TAX="trainset9_032012.pds.tax"
    export CONTAMINANTS="Chloroplast-Mitochondria-unknown-Archaea-Eukaryota"

    BATCH="\$(realpath "${batch_file}")"

    test -n "\$(find "\$DATA" -maxdepth 1 -name '*_R1.fastq' -print -quit)"
    test -s "\$REFERENCE_LOCATION/\$ALIGNREF"
    test -s "\$REFERENCE_LOCATION/\$TAXONREF_FASTA"
    test -s "\$REFERENCE_LOCATION/\$TAXONREF_TAX"

    echo "Starting mothur analysis..."

    cd "\$DATA"
    mothur "\$BATCH"

    cp final.tx.shared ../genus.shared
    cp final.tx.1.cons.taxonomy ../genus.taxonomy

    cd ..

    mkdir -p qc
    cp reads/*.logfile qc/ 2>/dev/null || true
    cp reads/*.summary qc/ 2>/dev/null || true
    cp reads/*.contigs_report qc/ 2>/dev/null || true
    cp reads/*.tax.summary qc/ 2>/dev/null || true
    cp reads/stability.files qc/ 2>/dev/null || true

    echo "Mothur analysis completed."
    """
}

process DIVERSITY {
    tag "Genus-level diversity"

    container "gut-microbiome:1.0"

    publishDir "${params.outdir}/diversity", mode: "copy"

    input:
    path shared
    path taxonomy
    path metadata
    path analysis_script

    output:
    path "*.tsv", emit: tables
    path "*.png", emit: figures

    script:
    """
    set -euo pipefail

    python "${analysis_script}" \
        --shared "${shared}" \
        --taxonomy "${taxonomy}" \
        --metadata "${metadata}" \
        --outdir .

    echo "Diversity analysis completed."
    """
}

workflow {
    raw_ch = Channel.value(file(params.raw_dir, checkIfExists: true))
    reference_ch = Channel.value(file(params.reference_dir, checkIfExists: true))
    batch_ch = Channel.value(file("${projectDir}/scripts/stability.batch", checkIfExists: true))
    metadata_ch = Channel.value(file(params.metadata, checkIfExists: true))
    script_ch = Channel.value(file("${projectDir}/scripts/diversity.py", checkIfExists: true))

    MOTHUR(raw_ch, reference_ch, batch_ch)

    DIVERSITY(
        MOTHUR.out.shared,
        MOTHUR.out.taxonomy,
        metadata_ch,
        script_ch
    )
}
