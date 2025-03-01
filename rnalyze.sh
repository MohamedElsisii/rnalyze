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

# Function to get the number of CPU cores in a platform-independent way
get_cpu_cores() {
    case "$(uname -s)" in
        Linux*)  nproc ;;
        Darwin*) sysctl -n hw.ncpu ;;
        CYGWIN*|MINGW32*|MSYS*|MINGW*) echo 1 ;;  
        *)       echo 1 ;;  
    esac
}

# Function to validate all required arguments
validate_required_args() {
    log "INFO" "🔍 Validating required arguments..."

    # Check if essential pipeline arguments are missing
    if [[ -z "${PIPELINE:-}" || -z "${TOOL:-}" || -z "${DATA_METHOD:-}" || -z "${DATA_TYPE:-}" || -z "${REF_METHOD:-}" ]]; then
        log "ERROR" "Missing required arguments. Please specify: --pipeline, --tool, --data, --type, and --ref."
        usage
        exit 1
    fi

    # Validate pipeline type
    if [[ "$PIPELINE" != "Full" && "$PIPELINE" != "Alignment" ]]; then
        log "ERROR" "Invalid pipeline type: $PIPELINE. Allowed values: Full, Alignment."
        usage
        exit 1
    fi

    # Validate alignment tool
    if [[ "$TOOL" != "BWA" && "$TOOL" != "Bowtie2" && "$TOOL" != "HISAT2" ]]; then
        log "ERROR" "Invalid tool specified: $TOOL. Allowed values: BWA, Bowtie2, HISAT2."
        usage
        exit 1
    fi

    # Validate data method and ensure required files/directories exist
    if [[ "$DATA_METHOD" == "Download" && -z "${SRR_PATH:-}" ]]; then
        log "ERROR" "Missing SRR list (--srr-path). Required when using --data Download."
        usage
        exit 1
    elif [[ "$DATA_METHOD" == "Directory" && -z "${DATA_DIR:-}" ]]; then
        log "ERROR" "Missing data directory (--data-dir). Required when using --data Directory."
        usage
        exit 1
    elif [[ "$DATA_METHOD" == "Directory" && ! -d "$DATA_DIR" ]]; then
        log "ERROR" "Specified data directory does not exist: $DATA_DIR."
        exit 1
    fi

    # Validate RNA-Seq data type
    if [[ "$DATA_TYPE" != "SE" && "$DATA_TYPE" != "PE" ]]; then
        log "ERROR" "Invalid data type: $DATA_TYPE. Allowed values: SE (Single-End), PE (Paired-End)."
        usage
        exit 1
    fi

    # Validate reference genome method and ensure required paths/URLs exist
    case "$REF_METHOD" in
        UnindexedURL|IndexedURL)
            if [[ -z "${REF_URL:-}" ]]; then
                log "ERROR" "Missing reference genome URL (--ref-url). Required for --ref $REF_METHOD."
                usage
                exit 1
            fi
            ;;
        UnindexedPath|IndexedPath)
            if [[ -z "${REF_PATH:-}" ]]; then
                log "ERROR" "Missing reference genome path (--ref-path). Required for --ref $REF_METHOD."
                usage
                exit 1
            elif [[ ! -e "$REF_PATH" ]]; then
                log "ERROR" "Reference genome file/directory does not exist: $REF_PATH."
                exit 1
            fi
            ;;
        *)
            log "ERROR" "Invalid reference method: $REF_METHOD. Allowed values: UnindexedURL, IndexedURL, UnindexedPath, IndexedPath."
            usage
            exit 1
            ;;
    esac

    # Validate trimming tool and parameters if trimming is enabled
    if [[ "$TRIM" == "yes" ]]; then
        if [[ -z "${TRIM_TOOL:-}" ]]; then
            log "ERROR" "--trim-tool is required when --trim is set to yes."
            usage
            exit 1
        fi

        if [[ "$TRIM_TOOL" == "Trimmomatic" ]]; then
            if [[ "$DATA_TYPE" == "SE" && -z "${ADAPTER_SE:-}" ]]; then
                log "ERROR" "--adapter-se is required when --trim-tool is Trimmomatic and --type is SE."
                usage
                exit 1
            elif [[ "$DATA_TYPE" == "PE" && -z "${ADAPTER_PE:-}" ]]; then
                log "ERROR" "--adapter-pe is required when --trim-tool is Trimmomatic and --type is PE."
                usage
                exit 1
            fi
        elif [[ "$TRIM_TOOL" == "Cutadapt" ]]; then
            if [[ "$DATA_TYPE" == "SE" && -z "${ADAPTER_SE:-}" && -z "${FRONT_SE:-}" ]]; then
                log "ERROR" "At least one of --adapter-se or --front-se is required when --trim-tool is Cutadapt and --type is SE."
                usage
                exit 1
            elif [[ "$DATA_TYPE" == "PE" ]]; then
                if [[ (-z "${ADAPTER_SE:-}" || -z "${ADAPTER_PE:-}") && (-z "${FRONT_SE:-}" || -z "${FRONT_PE:-}") ]]; then
                    log "ERROR" "For paired-end data, you must specify either both --adapter-se and --adapter-pe, or both --front-se and --front-pe, or all four options."
                    usage
                    exit 1
                fi
            fi
        else
            log "ERROR" "Invalid trimming tool specified. Use Trimmomatic or Cutadapt."
            usage
            exit 1
        fi
    fi

    # Validate GTF file if Full pipeline is selected
    if [[ "$PIPELINE" == "Full" ]]; then
        if [[ "$GTF_METHOD" == "Download" && -z "${GTF_URL:-}" ]]; then
            log "ERROR" "Missing GTF file URL (--gtf-url). Required when using --gtf Download."
            usage
            exit 1
        elif [[ "$GTF_METHOD" == "Path" && -z "${GTF_PATH:-}" ]]; then
            log "ERROR" "Missing GTF file path (--gtf-path). Required when using --gtf Path."
            usage
            exit 1
        elif [[ "$GTF_METHOD" == "Path" && ! -f "$GTF_PATH" ]]; then
            log "ERROR" "Specified GTF file does not exist: $GTF_PATH."
            exit 1
        fi
    fi

    # Validate CPU core count
    if [[ "$MAX_THREADS" -lt 1 ]]; then
        log "ERROR" "Invalid CPU core count detected: $MAX_THREADS. Must be at least 1."
        exit 1
    fi

    log "INFO" "✅ All required arguments validated successfully."
}



# Function to create directories
create_directories() {
    log "INFO" "📂 Creating directories..."
    mkdir -p "Data" "Mapping" "Ref_Genome"

    if [[ "$TRIM" == "yes" ]]; then
        mkdir -p "Fastqc_Reports/Before_Trimming" "Fastqc_Reports/After_Trimming" \
                 "MultiQC_Reports/Before_Trimming" "MultiQC_Reports/After_Trimming"
    else
        mkdir -p "Fastqc_Reports" "MultiQC_Reports"
    fi

    if [[ "$PIPELINE" == "Full" ]]; then
        mkdir -p "Featurecount" "Annotation"
    fi
}

# Function to handle data download or move
handle_data() {
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
        
        # Use find with -maxdepth 1 to only process files in the top-level directory
        find "$DATA_DIR" -maxdepth 1 -type f -iname "*.fastq.gz" -exec cp -v {} "$LOCATION/Data/" \; || { log "ERROR" "Failed to move .fastq.gz files. Exiting..."; exit 1; }
        
        log "INFO" "✅ .fastq.gz files moved successfully."
    fi
}

# Function to rename files with _raw extension if trimming is disabled
rename_files() {
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
                    log "ERROR" "Missing pair for $r1 or $2. Skipping renaming."
                fi
            done
        fi
    fi
}

# Function to handle reference genome
handle_reference_genome() {
    log "INFO" "📚 Handling reference genome..."
    case "$REF_METHOD" in
        "UnindexedURL")
            log "INFO" "Downloading reference genome from: $REF_URL"
            wget -q --show-progress -O "$LOCATION/Ref_Genome/ref_genome.fna.gz" "$REF_URL" || { log "ERROR" "Download failed. Exiting..."; exit 1; }
            if [[ "$REF_URL" == *.gz ]]; then
                log "INFO" "Uncompressing reference genome..."
                gunzip "$LOCATION/Ref_Genome/ref_genome.fna.gz" || { log "ERROR" "Failed to uncompress file. Exiting..."; exit 1; }
            fi
            # Ensure the file has a valid extension
            ref_file="$LOCATION/Ref_Genome/ref_genome.fna"
            if [[ ! -f "$ref_file" ]]; then
                log "ERROR" "Reference genome file not found or has an invalid extension. Expected .fna, .fa, or .fasta."
                exit 1
            fi
            index_reference_genome "$ref_file"
            ;;
        "IndexedURL")
            log "INFO" "Downloading indexed reference genome..."
            # Determine the expected index file name based on the selected tool
            case $TOOL in
                "BWA") expected_index="bwa_index.tar.gz" ;;
                "Bowtie2") expected_index="bowtie_index.tar.gz" ;;
                "HISAT2") expected_index="hisat2_index.tar.gz" ;;
                *) log "ERROR" "Invalid alignment tool specified for indexed reference genome."; exit 1 ;;
            esac
            # Download the indexed reference genome
            wget -q --show-progress -O "$LOCATION/Ref_Genome/$expected_index" "$REF_URL" || { log "ERROR" "Download failed. Exiting..."; exit 1; }
            log "INFO" "Extracting indexed reference genome..."
            tar -xzf "$LOCATION/Ref_Genome/$expected_index" -C "$LOCATION/Ref_Genome" || { log "ERROR" "Extraction failed. Exiting..."; exit 1; }
            rm "$LOCATION/Ref_Genome/$expected_index"
            log "INFO" "✅ Indexed reference genome extracted successfully."
            ;;
        "UnindexedPath")
            log "INFO" "Moving unindexed reference genome..."
            # Automatically detect the file name and validate its extension
            ref_file=$(basename "$REF_PATH")
            ref_ext="${ref_file##*.}"
            if [[ "$ref_ext" != "fna" && "$ref_ext" != "fa" && "$ref_ext" != "fasta" ]]; then
                log "ERROR" "Invalid reference genome extension: .$ref_ext. Expected .fna, .fa, or .fasta."
                exit 1
            fi
            # Copy the file and ensure it has a .fna extension
            cp "$REF_PATH" "$LOCATION/Ref_Genome/ref_genome.fna" || { log "ERROR" "Failed to move file. Exiting..."; exit 1; }
            index_reference_genome "$LOCATION/Ref_Genome/ref_genome.fna"
            ;;
        "IndexedPath")
            log "INFO" "Moving indexed reference genome..."
            # Validate that all required index files are present based on the selected tool
            case $TOOL in
                "BWA") required_files=("$REF_PATH"/*.{amb,ann,bwt,pac,sa}) ;;
                "Bowtie2") required_files=("$REF_PATH"/*.{1.bt2,2.bt2,3.bt2,4.bt2,rev.1.bt2,rev.2.bt2}) ;;
                "HISAT2") required_files=("$REF_PATH"/*.{1.ht2,2.ht2,3.ht2,4.ht2,5.ht2,6.ht2,7.ht2,8.ht2}) ;;
                *) log "ERROR" "Invalid alignment tool specified for indexed reference genome."; exit 1 ;;
            esac
            missing_files=()
            for file in "${required_files[@]}"; do
                if [[ ! -f "$file" ]]; then
                    missing_files+=("$(basename "$file")")
                fi
            done
            if [[ ${#missing_files[@]} -gt 0 ]]; then
                log "ERROR" "Missing required index files: ${missing_files[*]}"
                exit 1
            fi
            log "INFO" "✅ All required index files found."
            # Move all files to the Ref_Genome directory
            cp -r "$REF_PATH"/* "$LOCATION/Ref_Genome/" || { log "ERROR" "Failed to move files. Exiting..."; exit 1; }
            log "INFO" "✅ Indexed reference genome files moved successfully."
            ;;
    esac
}

# Function to index reference genome
index_reference_genome() {
    local ref_file=$1
    case $TOOL in
        "BWA")
            log "INFO" "📌 Indexing reference genome with BWA..."
            bwa index "$ref_file" || { log "ERROR" "Indexing failed. Exiting..."; exit 1; }
            ;;
        "Bowtie2")
            log "INFO" "📌 Indexing reference genome with Bowtie2..."
            bowtie2-build "$ref_file" "$LOCATION/Ref_Genome/ref_genome" || { log "ERROR" "Indexing failed. Exiting..."; exit 1; }
            ;;
        "HISAT2")
            log "INFO" "📌 Indexing reference genome with HISAT2..."
            hisat2-build "$ref_file" "$LOCATION/Ref_Genome/ref_genome" || { log "ERROR" "Indexing failed. Exiting..."; exit 1; }
            ;;
    esac
}

# Function to perform quality control and trimming
perform_quality_control_and_trimming() {
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
            log "INFO" "✂️ Running $TRIM_TOOL for trimming..."
            if [[ "$TRIM_TOOL" == "Trimmomatic" ]]; then
                # Single-End (SE) Trimming with Trimmomatic
                log "INFO" "⚙️ Trimming Single-End (SE) data with Trimmomatic..."
                for file in "$LOCATION/Data"/*.fastq.gz; do
                    filename=$(basename "$file" .fastq.gz)
                    log "INFO" "⚙️ Trimming $filename..."
                    trimmomatic SE -phred33 \
                        "$file" "$LOCATION/Data/${filename}_trimmed.fastq.gz" \
                        ILLUMINACLIP:"$CONDA_PREFIX/share/trimmomatic/adapters/${ADAPTER_SE:-TruSeq3-SE}.fa:2:30:10" \
                        LEADING:"${LEADING:-3}" TRAILING:"${TRAILING:-3}" SLIDINGWINDOW:"${SLIDING_WINDOW:-4:25}" MINLEN:"${MINLEN:-36}"
                    log "INFO" "✅ Trimmed $filename."
                done
            elif [[ "$TRIM_TOOL" == "Cutadapt" ]]; then
                # Single-End (SE) Trimming with Cutadapt
                log "INFO" "⚙️ Trimming Single-End (SE) data with Cutadapt..."
                for file in "$LOCATION/Data"/*.fastq.gz; do
                    filename=$(basename "$file" .fastq.gz)
                    log "INFO" "⚙️ Trimming $filename..."
                    cutadapt -a "${ADAPTER_SE}" -g "${FRONT_SE}" -q "${QUALITY_CUTOFF:-20}" -m "${MIN_READ_LENGTH:-30}" \
                        -o "$LOCATION/Data/${filename}_trimmed.fastq.gz" "$file"
                    log "INFO" "✅ Trimmed $filename."
                done
            fi
            log "INFO" "✅ Trimming completed."
        elif [[ "$DATA_TYPE" == "PE" && "$TRIM" == "yes" ]]; then
            log "INFO" "✂️ Running $TRIM_TOOL for trimming..."
            if [[ "$TRIM_TOOL" == "Trimmomatic" ]]; then
                # Paired-End (PE) Trimming with Trimmomatic
                log "INFO" "⚙️ Trimming Paired-End (PE) data with Trimmomatic..."
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
                done
            elif [[ "$TRIM_TOOL" == "Cutadapt" ]]; then
                # Paired-End (PE) Trimming with Cutadapt
                log "INFO" "⚙️ Trimming Paired-End (PE) data with Cutadapt..."
                for r1 in "$LOCATION/Data/"*_1.fastq.gz; do
                    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
                    filename=$(basename "$r1" _1.fastq.gz)
                    log "INFO" "⚙️ Trimming $filename..."
                    cutadapt -a "${ADAPTER_SE}" -A "${ADAPTER_PE}" -g "${FRONT_SE}" -G "${FRONT_PE}" -q "${QUALITY_CUTOFF:-20}" -m "${MIN_READ_LENGTH:-30}" \
                        -o "$LOCATION/Data/${filename}_1_trimmed.fastq.gz" -p "$LOCATION/Data/${filename}_2_trimmed.fastq.gz" "$r1" "$r2"
                    log "INFO" "✅ Trimmed $filename."
                done
            fi
            log "INFO" "✅ Trimming completed."
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
}

# Function to perform alignment, sorting, and indexing
perform_alignment() {
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
}

# Function to handle GTF file
handle_gtf() {
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
                # Automatically detect the file name and ensure it has a .gtf extension
                gtf_file=$(basename "$GTF_PATH")
                gtf_ext="${gtf_file##*.}"
                if [[ "$gtf_ext" != "gtf" ]]; then
                    cp "$GTF_PATH" "$LOCATION/Annotation/annotation.gtf" || { log "ERROR" "Failed to move file. Exiting..."; exit 1; }
                else
                    cp "$GTF_PATH" "$LOCATION/Annotation/annotation.gtf" || { log "ERROR" "Failed to move file. Exiting..."; exit 1; }
                fi
                ;;
        esac
    fi
}

# Function to perform feature counting
perform_feature_counting() {
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

        # Automatically merge the count files
        log "INFO" "⚙️ Merging count files..."

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

            log "INFO" "✅ Merged count matrix created at $preprocessed_matrix."
        fi
    fi
}

# Main script execution
start_time=$(date +%s)

# Initialize default values
TRIM="no"
IDENTIFIER="gene_id"
MAX_THREADS=$(( $(get_cpu_cores) - 1 ))
EXECUTION_TIME=$(date +"%Y-%m-%d-%H-%M-%S")
LOCATION="$PWD/RNAlyze-$EXECUTION_TIME"

# Function to display usage instructions
usage() {
    cat <<EOF

Usage: rnalyze [OPTIONS]

RNA-Seq Analysis Pipeline: A comprehensive tool for RNA-Seq data processing, alignment, and analysis.

Version: 1.0.1 [1/3/2025]

Options:
┌───────────────────┬──────────────────────────────────────────────────────────────────────────────┐
│ Option            │ Description                                                                  │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -p, --pipeline    │ Specify the pipeline type:                                                   │
│                   │   - Full: Perform full analysis (trimming, alignment, counting).             │
│                   │   - Alignment: Perform only alignment (no  counting).                        │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -t, --tool        │ Specify the alignment tool:                                                  │
│                   │   - BWA: Burrows-Wheeler Aligner.                                            │
│                   │   - Bowtie2: Bowtie2 aligner.                                                │
│                   │   - HISAT2: Hierarchical Indexing for Spliced Transcripts.                   │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -d, --data        │ Specify how to handle input data:                                            │
│                   │   - Download: Download data using SRR accession numbers.                     │
│                   │   - Directory: Use data from a local directory.                              │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -s, --srr-path    │ Path to SRR_Acc_List.txt (required if --data Download).                      │
│                   │ This file should contain a list of SRR accession numbers.                    │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -D, --data-dir    │ Directory containing input data (required if --data Directory).              │
│                   │ The directory should contain .fastq.gz files.                                │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -y, --type        │ Specify the RNA-Seq data type:                                               │
│                   │   - SE: Single-End data.                                                     │
│                   │   - PE: Paired-End data.                                                     │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -r, --ref         │ Specify how to handle the reference genome:                                  │
│                   │   - UnindexedURL: Download and index an unindexed reference genome from URL. │
│                   │   - IndexedURL: Download a pre-indexed reference genome from URL.            │
│                   │   - IndexedPath: Use a pre-indexed reference genome from a local path.       │
│                   │   - UnindexedPath: Use an unindexed reference genome from a local path.      │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -u, --ref-url     │ URL for downloading the reference genome (required if --ref UnindexedURL or  │
│                   │ IndexedURL).                                                                 │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -R, --ref-path    │ Path to the reference genome files (required if --ref IndexedPath or         │
│                   │ UnindexedPath).                                                              │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -T, --trim        │ Enable or disable trimming:                                                  │
│                   │   - yes: Enable trimming.                                                    │
│                   │   - no: Disable trimming.                                                    │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -x, --trim-tool   │ Specify the trimming tool (required if --trim yes):                          │
│                   │   - Trimmomatic: Use Trimmomatic for trimming.                               │
│                   │   - Cutadapt: Use Cutadapt for trimming.                                     │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -f, --gtf         │ Specify how to handle the GTF file:                                          │
│                   │   - Download: Download the GTF file from a URL.                              │
│                   │   - Path: Use a GTF file from a local path.                                  │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -F, --gtf-url     │ URL for downloading the GTF file (required if --gtf Download).               │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -P, --gtf-path    │ Path to the GTF file (required if --gtf Path).                               │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -i, --identifier  │ Gene identifier attribute in the GTF file (default: gene_id).                │
│                   │ This is used for feature counting.                                           │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -h, --help        │ Display this help message and exit.                                          │
└───────────────────┴──────────────────────────────────────────────────────────────────────────────┘

Trimmomatic Parameters (used when --trim-tool Trimmomatic):
┌───────────────────┬──────────────────────────────────────────────────────────────────────────────┐
│ Option            │ Description                                                                  │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -l, --leading     │ Remove bases with quality below this threshold from the start of the read.   │
│                   │ Default: 3.                                                                  │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -L, --trailing    │ Remove bases with quality below this threshold from the end of the read.     │
│                   │ Default: 3.                                                                  │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -w, --sliding-win │ Perform sliding window trimming with window size and average quality.        │
│                   │ Format: WINDOW_SIZE:QUALITY_THRESHOLD (e.g., 4:25).                          │
│                   │ Default: 4:25.                                                               │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -m, --minlen      │ Discard reads shorter than this length after trimming.                       │
│                   │ Default: 36.                                                                 │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -a, --adapter-se  │ Adapter sequence for Single-End data (required if --type SE):                │
│                   │   - TruSeq2-SE: TruSeq2 Single-End adapter.                                  │
│                   │   - TruSeq3-SE: TruSeq3 Single-End adapter.                                  │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -A, --adapter-pe  │ Adapter sequence for Paired-End data (required if --type PE):                │
│                   │   - NexteraPE-PE: Nextera Paired-End adapter.                                │
│                   │   - TruSeq2-PE: TruSeq2 Paired-End adapter.                                  │
│                   │   - TruSeq3-PE: TruSeq3 Paired-End adapter.                                  │
│                   │   - TruSeq3-PE-2: TruSeq3 Paired-End adapter.                                │
└───────────────────┴──────────────────────────────────────────────────────────────────────────────┘

Cutadapt Parameters (used when --trim-tool Cutadapt):
┌───────────────────┬──────────────────────────────────────────────────────────────────────────────┐
│ Option            │ Description                                                                  │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -a, --adapter-se  │ 3' adapter sequence to trim from single-end reads. (required if --type SE).  │
│                   │ Example: AGATCGGAAGAGC.                                                      │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -g, --front-se    │ 5' adapter sequence to trim from single-end reads. (optional if --type SE).  │
│                   │ Example: AGATCGGAAGAGC.                                                      │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -A, --adapter-pe  │ 3' adapter sequence for the second read in paired-end data.                  │
│                   │ (required if --type PE and --adapter-se is specified).                       │
│                   │ Example: AGATCGGAAGAGC.                                                      │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -G, --front-pe    │ 5' adapter sequence for the second read in paired-end data.                  │
│                   │ (required if --type PE and --front-se is specified).                         │
│                   │ Example: AGATCGGAAGAGC.                                                      │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -q, --quality-cut │ Trim bases with quality below this threshold.                                │
│                   │ Default: 20.                                                                 │
├───────────────────┼──────────────────────────────────────────────────────────────────────────────┤
│ -M, --min-read-len│ Discard reads shorter than this length after trimming.                       │
│                   │ Default: 30.                                                                 │
└───────────────────┴──────────────────────────────────────────────────────────────────────────────┘

Examples:
 
1. Full Pipeline with Trimmomatic (Single-End Data)
   
   - Objective: Analyze single-end RNA-Seq data using the full pipeline (trimming, alignment, and counting).
  
   - Command:
     
     rnalyze -p Full -t HISAT2 -d Download -s /path/to/SRR_Acc_List.txt -y SE -r UnindexedURL -u http://example.com/ref.fna -T yes -x Trimmomatic -a TruSeq3-SE -l 5 -f Download -F http://example.com/annotation.gtf -i gene_id
   
   - Explanation:
     - Download single-end RNA-Seq data using SRR accession numbers.
     - Use HISAT2 for alignment.
     - Trim reads using Trimmomatic with the TruSeq3-SE adapter.
     - Download and index the reference genome from a URL.
     - Download the GTF file for annotation.
     - Perform feature counting using the gene_id attribute.


2. Alignment Pipeline with Cutadapt (Paired-End Data)
  
   - Objective: Analyze paired-end RNA-Seq data using the alignment pipeline (trimming and alignment only).
  
   - Command:
    
     rnalyze -p Alignment -t BWA -d Directory -D /path/to/data -y PE -r IndexedPath -R /path/to/indexed_genome -T yes -x Cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20 -M 30
   
   - Explanation:
     - Use paired-end RNA-Seq data from a local directory.
     - Use BWA for alignment.
     - Trim reads using Cutadapt with 3' adapters for both reads.
     - Use a pre-indexed reference genome from a local path.
     - Set a quality cutoff of 20 and a minimum read length of 30 after trimming.


3. Full Pipeline with Cutadapt (Paired-End Data)
   
   - Objective: Analyze paired-end RNA-Seq data using the full pipeline (trimming, alignment, and counting).
   
   - Command:
    
     rnalyze -p Full -t HISAT2 -d Directory -D /path/to/data -y PE -r UnindexedPath -R /path/to/ref_genome.fna -T yes -x Cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g AGATCGGAAGAGC -G AGATCGGAAGAGC -q 20 -M 30 -f Path -P /path/to/annotation.gtf -i gene_id
  
   - Explanation:
     - Use paired-end RNA-Seq data from a local directory.
     - Use HISAT2 for alignment.
     - Trim reads using Cutadapt with both 3' and 5' adapters for both reads.
     - Use an unindexed reference genome from a local path and index it.
     - Use a GTF file from a local path for annotation.
     - Perform feature counting using the gene_id attribute.


4. Alignment Pipeline with Trimmomatic (Single-End Data)
   
   - Objective: Analyze single-end RNA-Seq data using the alignment pipeline (trimming and alignment only).
   
   - Command:
     
     rnalyze -p Alignment -t Bowtie2 -d Download -s /path/to/SRR_Acc_List.txt -y SE -r IndexedURL -u http://example.com/indexed_genome.tar.gz -T yes -x Trimmomatic -a TruSeq3-SE -l 5
   
   - Explanation:
     - Download single-end RNA-Seq data using SRR accession numbers.
     - Use Bowtie2 for alignment.
     - Trim reads using Trimmomatic with the TruSeq3-SE adapter.
     - Use a pre-indexed reference genome downloaded from a URL.
EOF
}


# Parse command-line arguments
while getopts ":p:t:d:s:D:y:r:u:R:T:x:l:L:w:m:a:A:g:G:q:M:f:F:P:i:h" opt; do
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
        x) TRIM_TOOL="$OPTARG" ;;
        l) LEADING="$OPTARG" ;;
        L) TRAILING="$OPTARG" ;;
        w) SLIDING_WINDOW="$OPTARG" ;;
        m) MINLEN="$OPTARG" ;;
        a) ADAPTER_SE="$OPTARG" ;;
        A) ADAPTER_PE="$OPTARG" ;;
        g) FRONT_SE="$OPTARG" ;;
        G) FRONT_PE="$OPTARG" ;;
        q) QUALITY_CUTOFF="$OPTARG" ;;
        M) MIN_READ_LENGTH="$OPTARG" ;;
        f) GTF_METHOD="$OPTARG" ;;
        F) GTF_URL="$OPTARG" ;;
        P) GTF_PATH="$OPTARG" ;;
        i) IDENTIFIER="$OPTARG" ;;
        h) usage; exit 0 ;;
        \?) log "ERROR" "Invalid option: -$OPTARG"; usage; exit 1 ;;
        :) log "ERROR" "Option -$OPTARG requires an argument."; usage; exit 1 ;;
    esac
done

# Validate required arguments
validate_required_args

# Create the working directory
mkdir -p "$LOCATION"
cd "$LOCATION"

# Create necessary directories
create_directories

# Handle data based on the data method
handle_data

# Rename files with _raw extension if trimming is disabled
rename_files

# Handle reference genome
handle_reference_genome

# Perform quality control and trimming
perform_quality_control_and_trimming

# Perform alignment, sorting, and indexing
perform_alignment

# Handle GTF file
handle_gtf

# Perform feature counting
perform_feature_counting

# End the timer
end_time=$(date +%s)
processing_time=$((end_time - start_time))
log "INFO" "⏰ End Time: $(date)"
log "INFO" "⏱️ Total Processing Time: $((processing_time / 60)) minutes and $((processing_time % 60)) seconds."
log "INFO" "✅ Pipeline execution completed!"
