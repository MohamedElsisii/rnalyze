#!/bin/bash

# Exit immediately if any command fails, unset variables are used, or a pipeline fails
set -euo pipefail

# Logging function
log() {
    local level=$1
    local message=$2
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[$timestamp] [$level] $message"
}

# Initialize default values
TRIM="no"
IDENTIFIER="gene_id"
MAX_THREADS=$(( $(nproc) - 1 ))
EXECUTION_TIME=$(date +"%Y-%m-%d-%H-%M-%S")
LOCATION="$PWD/RNAlyze-$EXECUTION_TIME"

# Function to display usage instructions
usage() {
    cat <<EOF
Usage: rnalyze [OPTIONS]

Options:
  -p, --pipeline <Full|Alignment>   Specify the pipeline type.
  -t, --tool <BWA|Bowtie2|HISAT2>                          Specify the alignment tool.
  -d, --data <Download|Directory>                          Specify data handling method.
  -s, --srr-path <PATH>                                    Path to SRR_Acc_List.txt (required if -d Download).
  -D, --data-dir <PATH>                                    Directory containing data (required if -d Directory).
  -y, --type <SE|PE>                                       Specify RNA-Seq data type (Single-End or Paired-End).
  -r, --ref <UnindexedURL|IndexedURL|IndexedPath|UnindexedPath>
                                                            Specify reference genome handling.
  -u, --ref-url <URL>                                      URL for downloading reference genome.
  -R, --ref-path <PATH>                                    Path to reference genome files.
  -T, --trim <yes|no>                                      Enable or disable trimming (only for Full pipeline).
  -l, --leading <INT>                                      Trimmomatic LEADING threshold (default: 3).
  -L, --trailing <INT>                                     Trimmomatic TRAILING threshold (default: 3).
  -w, --sliding-window <INT:INT>                           Trimmomatic SLIDINGWINDOW parameters (default: 4:25).
  -m, --minlen <INT>                                       Trimmomatic MINLEN threshold (default: 36).
  -a, --adapter-se <TruSeq2-SE|TruSeq3-SE>                 Adapter for Single-End data (required if --trim yes and --type SE).
  -A, --adapter-pe <NexteraPE-PE|TruSeq2-PE|TruSeq3-PE|TruSeq3-PE-2>
                                                            Adapter for Paired-End data (required if --trim yes and --type PE).
  -g, --gtf <Download|Path>                                Specify GTF file handling.
  -G, --gtf-url <URL>                                      URL for downloading GTF file.
  -P, --gtf-path <PATH>                                    Path to GTF file.
  -i, --identifier <STRING>                                Gene identifier attribute in GTF file (default: gene_id).
  -h, --help                                               Display this help message.

Examples:
  rnalyze -p Full -t HISAT2 -d Download -s /path/to/SRR_Acc_List.txt -y SE -r UnindexedURL -u http://example.com/ref.fna -T yes -a TruSeq3-SE -l 5 -g Download -G http://example.com/annotation.gtf -i gene_id
  rnalyze -p Alignment -t BWA -d Directory -D /path/to/data -y PE -r IndexedPath -R /path/to/indexed_genome
EOF
}

# Parse command-line arguments
while getopts ":p:t:d:s:D:y:r:u:R:T:l:L:w:m:a:A:g:G:P:i:h" opt; do
    case $opt in
        p) PIPELINE="$OPTARG" ;;
        t) TOOL="$OPTARG" ;;
        d) DATA_METHOD="$OPTARG" ;;
        s) SRR_DIR="$OPTARG" ;;
        D) DATA_DIR="$OPTARG" ;;
        y) DATA_TYPE="$OPTARG" ;;
        r) REF_METHOD="$OPTARG" ;;
        u) REF_URL="$OPTARG" ;;
        R) REF_PATH="$OPTARG" ;;
        T) TRIM="$OPTARG" ;;
        l) LEADING="$OPTARG" ;;
        L) TRAILING="$OPTARG" ;;
        w) SLIDING_WINDOW="$OPTARG" ;;
        m) MINLEN="$OPTARG" ;;
        a) ADAPTER_SE="$OPTARG" ;;
        A) ADAPTER_PE="$OPTARG" ;;
        g) GTF_METHOD="$OPTARG" ;;
        G) GTF_URL="$OPTARG" ;;
        P) GTF_PATH="$OPTARG" ;;
        i) IDENTIFIER="$OPTARG" ;;
        h) usage; exit 0 ;;
        \?) log "ERROR" "Invalid option: -$OPTARG"; usage; exit 1 ;;
        :) log "ERROR" "Option -$OPTARG requires an argument."; usage; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${PIPELINE:-}" || -z "${TOOL:-}" || -z "${DATA_METHOD:-}" || -z "${DATA_TYPE:-}" || -z "${REF_METHOD:-}" ]]; then
    log "ERROR" "Missing required arguments."
    usage
    exit 1
fi

# Validate pipeline-specific arguments
case "$PIPELINE" in
    Full)
        if [[ -z "${GTF_METHOD:-}" ]]; then
            log "ERROR" "For Full pipeline, --gtf argument is required."
            usage
            exit 1
        fi
        ;;
    Alignment)
        if [[ -n "${GTF_METHOD:-}" ]]; then
            log "ERROR" "For Alignment pipeline, --gtf argument should not be provided."
            usage
            exit 1
        fi
        ;;
esac

# Validate data handling method
if [[ "$DATA_METHOD" == "Download" && -z "${SRR_DIR:-}" ]]; then
    log "ERROR" "--srr-path is required when --data is set to Download."
    usage
    exit 1
elif [[ "$DATA_METHOD" == "Download" ]]; then
    SRR_PATH="$SRR_DIR/SRR_Acc_List.txt"
    if [[ ! -f "$SRR_PATH" ]]; then
        log "ERROR" "SRR_Acc_List.txt not found in the provided directory: $SRR_DIR"
        usage
        exit 1
    fi
    log "INFO" "✅ Found SRR_Acc_List.txt at $SRR_PATH"
elif [[ "$DATA_METHOD" == "Directory" && -z "${DATA_DIR:-}" ]]; then
    log "ERROR" "--data-dir is required when --data is set to Directory."
    usage
    exit 1
fi

# Validate reference genome handling method
if [[ "$REF_METHOD" == "UnindexedURL" && -z "${REF_URL:-}" ]]; then
    log "ERROR" "--ref-url is required when --ref is set to UnindexedURL."
    usage
    exit 1
elif [[ "$REF_METHOD" == "IndexedURL" && -z "${REF_URL:-}" ]]; then
    log "ERROR" "--ref-url is required when --ref is set to IndexedURL."
    usage
    exit 1
elif [[ "$REF_METHOD" == "IndexedPath" && -z "${REF_PATH:-}" ]]; then
    log "ERROR" "--ref-path is required when --ref is set to IndexedPath."
    usage
    exit 1
elif [[ "$REF_METHOD" == "UnindexedPath" && -z "${REF_PATH:-}" ]]; then
    log "ERROR" "--ref-path is required when --ref is set to UnindexedPath."
    usage
    exit 1
fi

# Validate GTF handling method
if [[ "${GTF_METHOD:-}" == "Download" && -z "${GTF_URL:-}" ]]; then
    log "ERROR" "--gtf-url is required when --gtf is set to Download."
    usage
    exit 1
elif [[ "${GTF_METHOD:-}" == "Path" && -z "${GTF_PATH:-}" ]]; then
    log "ERROR" "--gtf-path is required when --gtf is set to Path."
    usage
    exit 1
fi

# Validate trimming and adapter options (only for Full pipeline)
if [[ "$PIPELINE" == "Full" && "$TRIM" == "yes" ]]; then
    if [[ "$DATA_TYPE" == "SE" && -z "${ADAPTER_SE:-}" ]]; then
        log "ERROR" "--adapter-se is required when --trim is set to yes and --type is SE."
        usage
        exit 1
    elif [[ "$DATA_TYPE" == "PE" && -z "${ADAPTER_PE:-}" ]]; then
        log "ERROR" "--adapter-pe is required when --trim is set to yes and --type is PE."
        usage
        exit 1
    fi
fi

# Start the timer to measure the total execution time of the pipeline
start_time=$(date +%s)

# Create the working directory
mkdir -p "$LOCATION"
cd "$LOCATION"

# Create necessary directories for storing pipeline outputs
log "INFO" "📂 Creating directories..."
mkdir -p "Data" "Mapping" "Ref_Genome"

# Conditionally create trimming-related directories only if TRIM=yes
if [[ "$TRIM" == "yes" ]]; then
    mkdir -p "Fastqc_Reports/Before_Trimming" "Fastqc_Reports/After_Trimming" \
             "MultiQC_Reports/Before_Trimming" "MultiQC_Reports/After_Trimming"
else
    mkdir -p "Fastqc_Reports" "MultiQC_Reports"
fi

# Create Featurecount and Annotation directories only for Full pipeline
if [[ "$PIPELINE" == "Full" ]]; then
    mkdir -p "Featurecount" "Annotation"
fi

# Handle data based on the data method
if [[ "$DATA_METHOD" == "Download" ]]; then
    log "INFO" "📥 Downloading data using SRR_Acc_List.txt..."
    while read -r srr_id; do
        if [[ -n "$srr_id" ]]; then
            log "INFO" "⏳ Downloading $srr_id..."
            fastq-dump --gzip --split-files --outdir "$LOCATION/Data" "$srr_id" || { log "ERROR" "Failed to download $srr_id. Exiting..."; exit 1; }
            log "INFO" "✅ Downloaded $srr_id."
        fi
    done < "$SRR_PATH"

    log "INFO" "✅ Data download completed."
elif [[ "$DATA_METHOD" == "Directory" ]]; then
    log "INFO" "🚚 Moving .fastq.gz files from $DATA_DIR to $LOCATION/Data..."
    find "$DATA_DIR" -type f -iname "*.fastq.gz" -exec cp -v {} "$LOCATION/Data/" \; || { log "ERROR" "Failed to move .fastq.gz files. Exiting..."; exit 1; }
    log "INFO" "✅ .fastq.gz files moved successfully."
fi

# Rename files with _raw extension if trimming is disabled
if [[ "$TRIM" == "no" ]]; then
    log "INFO" "✂️ No trimming requested. Renaming files with _raw extension..."
    
    if [[ "$DATA_TYPE" == "SE" ]]; then
        # Single-End (SE) Data: Rename all .fastq.gz files with _raw
        for file in "$LOCATION/Data/"*.fastq.gz; do
            filename=$(basename "$file" .fastq.gz)
            mv "$file" "$LOCATION/Data/${filename}_raw.fastq.gz"
        done
    elif [[ "$DATA_TYPE" == "PE" ]]; then
        # Paired-End (PE) Data: Rename both _1.fastq.gz and _2.fastq.gz files with _1_raw and _2_raw
        for r1 in "$LOCATION/Data/"*_1.fastq.gz; do
            r2="${r1/_1.fastq.gz/_2.fastq.gz}"
            
            if [[ -f "$r1" && -f "$r2" ]]; then
                # Extract base name without _1 or _2
                sample=$(basename "$r1" _1.fastq.gz)
                
                # Rename R1 and R2 files
                mv "$r1" "$LOCATION/Data/${sample}_1_raw.fastq.gz"
                mv "$r2" "$LOCATION/Data/${sample}_2_raw.fastq.gz"
            else
                log "ERROR" "Missing pair for $r1 or $r2. Skipping renaming."
            fi
        done
    fi
fi

# Reference Genome Handling
log "INFO" "📚 Handling reference genome..."
case "$REF_METHOD" in
    "UnindexedURL")
        log "INFO" "Downloading reference genome from: $REF_URL"
        wget -q --show-progress -O "$LOCATION/Ref_Genome/ref_genome.fna.gz" "$REF_URL" || { log "ERROR" "Download failed. Exiting..."; exit 1; }
        if [[ "$REF_URL" == *.gz ]]; then
            log "INFO" "Uncompressing reference genome..."
            gunzip "$LOCATION/Ref_Genome/ref_genome.fna.gz" || { log "ERROR" "Failed to uncompress file. Exiting..."; exit 1; }
        fi
        case $TOOL in
            "BWA")
                log "INFO" "📌 Indexing reference genome with BWA..."
                bwa index "$LOCATION/Ref_Genome/ref_genome.fna" || { log "ERROR" "Indexing failed. Exiting..."; exit 1; }
                ;;
            "Bowtie2")
                log "INFO" "📌 Indexing reference genome with Bowtie2..."
                bowtie2-build "$LOCATION/Ref_Genome/ref_genome.fna" "$LOCATION/Ref_Genome/ref_genome" || { log "ERROR" "Indexing failed. Exiting..."; exit 1; }
                ;;
            "HISAT2")
                log "INFO" "📌 Indexing reference genome with HISAT2..."
                hisat2-build "$LOCATION/Ref_Genome/ref_genome.fna" "$LOCATION/Ref_Genome/ref_genome" || { log "ERROR" "Indexing failed. Exiting..."; exit 1; }
                ;;
        esac
        ;;
    "IndexedURL")
        log "INFO" "Downloading indexed reference genome..."
        wget -q --show-progress -O "$LOCATION/Ref_Genome/indexed_ref_genome.tar.gz" "$REF_URL" || { log "ERROR" "Download failed. Exiting..."; exit 1; }
        log "INFO" "Extracting indexed reference genome..."
        tar -xzf "$LOCATION/Ref_Genome/indexed_ref_genome.tar.gz" -C "$LOCATION/Ref_Genome" || { log "ERROR" "Extraction failed. Exiting..."; exit 1; }
        rm "$LOCATION/Ref_Genome/indexed_ref_genome.tar.gz"
        ;;
    "UnindexedPath")
        log "INFO" "Moving unindexed reference genome..."
        cp "$REF_PATH" "$LOCATION/Ref_Genome/ref_genome.fna" || { log "ERROR" "Failed to move file. Exiting..."; exit 1; }
        case $TOOL in
            "BWA")
                log "INFO" "📌 Indexing reference genome with BWA..."
                bwa index "$LOCATION/Ref_Genome/ref_genome.fna" || { log "ERROR" "Indexing failed. Exiting..."; exit 1; }
                ;;
            "Bowtie2")
                log "INFO" "📌 Indexing reference genome with Bowtie2..."
                bowtie2-build "$LOCATION/Ref_Genome/ref_genome.fna" "$LOCATION/Ref_Genome/ref_genome" || { log "ERROR" "Indexing failed. Exiting..."; exit 1; }
                ;;
            "HISAT2")
                log "INFO" "📌 Indexing reference genome with HISAT2..."
                hisat2-build "$LOCATION/Ref_Genome/ref_genome.fna" "$LOCATION/Ref_Genome/ref_genome" || { log "ERROR" "Indexing failed. Exiting..."; exit 1; }
                ;;
        esac
        ;;
    "IndexedPath")
        log "INFO" "Moving indexed reference genome..."
        cp -r "$REF_PATH"/* "$LOCATION/Ref_Genome/" || { log "ERROR" "Failed to move files. Exiting..."; exit 1; }
        ;;
esac

# Quality Control & Trimming
if [[ "$PIPELINE" == "Full" || "$PIPELINE" == "Alignment" ]]; then
    log "INFO" "🔍 Running FastQC on raw data..."
    for file in "$LOCATION/Data/"*fastq.gz; do
        log "INFO" "⚙️ Processing $file..."
        filename=$(basename "$file" .fastq.gz)
        if [[ "$TRIM" == "yes" ]]; then
            fastqc "$file" -o "$LOCATION/Fastqc_Reports/Before_Trimming"
        else
            fastqc "$file" -o "$LOCATION/Fastqc_Reports"
        fi
    done

    # Generate MultiQC report for raw data
    if [[ "$TRIM" == "yes" ]]; then
        multiqc "$LOCATION/Fastqc_Reports/Before_Trimming" -o "$LOCATION/MultiQC_Reports/Before_Trimming"
    else
        multiqc "$LOCATION/Fastqc_Reports" -o "$LOCATION/MultiQC_Reports"
    fi

    if [[ "$DATA_TYPE" == "SE" && "$TRIM" == "yes" ]]; then
        log "INFO" "✂️ Running Trimmomatic for trimming..."
        # Single-End (SE) Trimming
        log "INFO" "⚙️ Trimming Single-End (SE) data..."
        for file in "$LOCATION/Data"/*.fastq.gz; do
            filename=$(basename "$file" .fastq.gz)
            log "INFO" "⚙️ Trimming $filename..."
            trimmomatic SE -phred33 \
                "$file" "$LOCATION/Data/${filename}_trimmed.fastq.gz" \
                ILLUMINACLIP:"$CONDA_PREFIX/share/trimmomatic/adapters/${ADAPTER_SE:-TruSeq3-SE}.fa:2:30:10" \
                LEADING:"${LEADING:-3}" TRAILING:"${TRAILING:-3}" SLIDINGWINDOW:"${SLIDING_WINDOW:-4:25}" MINLEN:"${MINLEN:-36}"
            log "INFO" "✅ Trimmed $filename."
            log "INFO" "✅ Trimming completed."
        done
    elif [[ "$DATA_TYPE" == "PE" && "$TRIM" == "yes" ]]; then
        # Paired-End (PE) Trimming
        log "INFO" "✂️ Running Trimmomatic for trimming..."
        log "INFO" "⚙️ Trimming Paired-End (PE) data..."
        for r1 in "$LOCATION/Data/"*_1.fastq.gz; do
            r2="${r1/_1.fastq.gz/_2.fastq.gz}"
            filename=$(basename "$r1" _1.fastq.gz)
            log "INFO" "⚙️ Trimming $filename..."
            trimmomatic PE -phred33 \
                "$r1" "$r2" \
                "$LOCATION/Data/${filename}_1_trimmed.fastq.gz" "$LOCATION/Data/${filename}_1_unpaired.fastq.gz" \
                "$LOCATION/Data/${filename}_2_trimmed.fastq.gz" "$LOCATION/Data/${filename}_2_unpaired.fastq.gz" \
                ILLUMINACLIP:"$CONDA_PREFIX/share/trimmomatic/adapters/${ADAPTER_PE:-TruSeq3-PE}.fa:2:30:10" \
                LEADING:"${LEADING:-3}" TRAILING:"${TRAILING:-3}" SLIDINGWINDOW:"${SLIDING_WINDOW:-4:25}" MINLEN:"${MINLEN:-36}"
            log "INFO" "✅ Trimmed $filename."
            log "INFO" "✅ Trimming completed."
        done
    fi

    # Run FastQC on trimmed files (only if TRIM=yes)
    if [[ "$TRIM" == "yes" ]]; then
        log "INFO" "🔍 Running FastQC on trimmed files..."
        for file in "$LOCATION/Data"/*_trimmed.fastq.gz; do
            log "INFO" "⚙️ Processing $file..."
            filename=$(basename "$file" .fastq.gz)
            fastqc "$file" -o "$LOCATION/Fastqc_Reports/After_Trimming"
        done

        # Generate MultiQC report for trimmed data
        multiqc "$LOCATION/Fastqc_Reports/After_Trimming" -o "$LOCATION/MultiQC_Reports/After_Trimming"
    fi
fi

# Alignment, Sorting, and Indexing
log "INFO" "🎯 Performing sequence alignment..."

# Determine the file extension based on TRIM value
if [[ "$TRIM" == "yes" ]]; then
    EXTENSION="trimmed"
else
    EXTENSION="raw"
fi

if [[ "$DATA_TYPE" == "SE" ]]; then
    # Single-End (SE) Alignment
    for file in "$LOCATION/Data"/*_${EXTENSION}.fastq.gz; do
        filename=$(basename "$file" .fastq.gz)
        log "INFO" "⚙️ Aligning ${filename} with $TOOL..."
        case $TOOL in
            "BWA")
                bwa mem -t "$MAX_THREADS" "$LOCATION/Ref_Genome/ref_genome.fna" "$file" > "$LOCATION/Mapping/${filename}.sam"
                ;;
            "Bowtie2")
                bowtie2 -p "$MAX_THREADS" -x "$LOCATION/Ref_Genome/ref_genome" -U "$file" -S "$LOCATION/Mapping/${filename}.sam"
                ;;
            "HISAT2")
                hisat2 -p "$MAX_THREADS" -x "$LOCATION/Ref_Genome/ref_genome" -U "$file" -S "$LOCATION/Mapping/${filename}.sam"
                ;;
            *)
                log "ERROR" "Invalid aligner. Exiting..."
                exit 1
                ;;
        esac
        # Convert SAM to BAM, sort, and index
        samtools view -bS "$LOCATION/Mapping/${filename}.sam" > "$LOCATION/Mapping/${filename}.bam"
        samtools sort --threads "$((MAX_THREADS - 1))" -o "$LOCATION/Mapping/${filename}.sorted.bam" "$LOCATION/Mapping/${filename}.bam"
        samtools index "$LOCATION/Mapping/${filename}.sorted.bam"
        samtools flagstat "$LOCATION/Mapping/${filename}.sorted.bam" > "$LOCATION/Mapping/${filename}_Flagstated.txt"
        log "INFO" "✅ Alignment completed for ${filename}."
    done
elif [[ "$DATA_TYPE" == "PE" ]]; then
    # Paired-End (PE) Alignment
    for r1 in "$LOCATION/Data"/*_1_${EXTENSION}.fastq.gz; do
        # Extract base name without the _1_trimmed or _1_raw suffix
        if [[ "$r1" == *"_1_${EXTENSION}.fastq.gz" ]]; then
            sample=$(basename "$r1" _1_${EXTENSION}.fastq.gz)

            # Determine the corresponding R2 file
            r2="$LOCATION/Data/${sample}_2_${EXTENSION}.fastq.gz"
            if [[ ! -f "$r2" ]]; then
                log "ERROR" "No matching R2 file found for $r1. Skipping alignment."
                continue
            fi

            # Perform alignment
            log "INFO" "⚙️ Aligning $r1 (R1) and $r2 (R2) with $TOOL..."
            case $TOOL in
                "BWA")
                    bwa mem -t "$MAX_THREADS" "$LOCATION/Ref_Genome/ref_genome.fna" "$r1" "$r2" > "$LOCATION/Mapping/${sample}.sam"
                    ;;
                "Bowtie2")
                    bowtie2 -p "$MAX_THREADS" -x "$LOCATION/Ref_Genome/ref_genome" -1 "$r1" -2 "$r2" -S "$LOCATION/Mapping/${sample}.sam"
                    ;;
                "HISAT2")
                    hisat2 -p "$MAX_THREADS" -x "$LOCATION/Ref_Genome/ref_genome" -1 "$r1" -2 "$r2" -S "$LOCATION/Mapping/${sample}.sam"
                    ;;
                *)
                    log "ERROR" "Invalid aligner. Exiting..."
                    exit 1
                    ;;
            esac
            # Convert SAM to BAM, sort, and index
            samtools view -bS "$LOCATION/Mapping/${sample}.sam" > "$LOCATION/Mapping/${sample}.bam"
            samtools sort --threads "$((MAX_THREADS - 1))" -o "$LOCATION/Mapping/${sample}.sorted.bam" "$LOCATION/Mapping/${sample}.bam"
            samtools index "$LOCATION/Mapping/${sample}.sorted.bam"
            # Generate flagstat report
            samtools flagstat "$LOCATION/Mapping/${sample}.sorted.bam" > "$LOCATION/Mapping/${sample}_Flagstated.txt"
            log "INFO" "✅ Alignment completed for $sample."
        fi
    done
fi

# Handle GTF file (only for Full pipeline)
if [[ "$PIPELINE" == "Full" ]]; then
    log "INFO" "📚 Handling GTF file..."
    case "${GTF_METHOD:-}" in
        "Download")
            log "INFO" "Downloading GTF file..."
            wget -q -O "$LOCATION/Annotation/annotation.gtf.gz" "${GTF_URL:-}" || { log "ERROR" "Download failed. Exiting..."; exit 1; }
            if [[ "${GTF_URL:-}" == *.gz ]]; then
                log "INFO" "Uncompressing GTF file..."
                gunzip "$LOCATION/Annotation/annotation.gtf.gz" || { log "ERROR" "Failed to uncompress file. Exiting..."; exit 1; }
            fi
            ;;
        "Path")
            log "INFO" "Moving GTF file..."
            cp "${GTF_PATH:-}" "$LOCATION/Annotation/annotation.gtf" || { log "ERROR" "Failed to move file. Exiting..."; exit 1; }
            ;;
    esac
fi

# Feature Counting (only for Full pipeline)
if [[ "$PIPELINE" == "Full" ]]; then
    log "INFO" "📊 Running FeatureCounts"
    for file in "$LOCATION/Mapping"/*.sorted.bam; do
        filename=$(basename "$file" .sorted.bam)

        if [[ "$DATA_TYPE" == "SE" ]]; then
            # Single-End (SE) Pipeline
            featureCounts -T "$MAX_THREADS" -g "$IDENTIFIER" -a "$LOCATION/Annotation/annotation.gtf" \
                -o "$LOCATION/Featurecount/$filename.counts.txt" "$file" || { log "ERROR" "FeatureCounts failed for $file. Exiting..."; exit 1; }
        elif [[ "$DATA_TYPE" == "PE" ]]; then
            # Paired-End (PE) Pipeline
            featureCounts -T "$MAX_THREADS" -g "$IDENTIFIER" -a "$LOCATION/Annotation/annotation.gtf" \
                -p -o "$LOCATION/Featurecount/$filename.counts.txt" "$file" || { log "ERROR" "FeatureCounts failed for $file. Exiting..."; exit 1; }
        fi
    done
    log "INFO" "✅ FeatureCounts completed."

    # Automatically merge and preprocess the count files
    log "INFO" "⚙️ Merging and preprocessing count files..."

    # Check if any .counts.txt files exist
    shopt -s nullglob  # Ensure globbing returns an empty list if no files match
    count_files=("$LOCATION/Featurecount"/*.counts.txt)
    if [[ ${#count_files[@]} -eq 0 ]]; then
        log "ERROR" "No count files found in $LOCATION/Featurecount. Skipping merging..."
    else
        # Generate the merged count matrix using a single awk command
        merged_matrix="$LOCATION/Featurecount/merged_counts.tsv"
        preprocessed_matrix="$LOCATION/Featurecount/preprocessed_counts.tsv"

        log "INFO" "📊 Generating merged count matrix and consolidating duplicate genes..."

        # Use awk to process all files in one pass
        awk -v identifier="$IDENTIFIER" '
        BEGIN {
            FS = OFS = "\t"  # Set input and output field separators to tab
        }
        FNR > 2 {  # Skip header lines in each file
            gene = $1  # Extract the identifier (first column)
            count = $NF  # Extract the count (last column)
            sample = FILENAME  # Get the file name

            # Remove the directory path from the filename
            sub(".*/", "", sample)

            # Remove the .counts.txt suffix
            sub(/\.counts\.txt$/, "", sample)

            # Use a concatenated key for gene-sample pairs
            key = gene "," sample
            counts[key] += count  # Accumulate counts for each gene-sample pair

            # Track unique genes and samples
            genes[gene]
            samples[sample]
        }
        END {
            # Print the header row
            printf "%s", identifier
            for (sample in samples) {
                printf "\t%s", sample
            }
            print ""

            # Print the gene counts
            for (gene in genes) {
                printf "%s", gene
                for (sample in samples) {
                    key = gene "," sample
                    printf "\t%d", counts[key]
                }
                print ""
            }
        }
        ' "${count_files[@]}" | sort > "$preprocessed_matrix"

        log "INFO" "✅ Merged and preprocessed count matrix created at $preprocessed_matrix."
    fi
fi

# End the timer
end_time=$(date +%s)
processing_time=$((end_time - start_time))
log "INFO" "⏰ End Time: $(date)"
log "INFO" "⏱️ Total Processing Time: $((processing_time / 60)) minutes and $((processing_time % 60)) seconds."
log "INFO" "✅ Pipeline execution completed!"
