# RNAlyze - Bioconda

### 🆕 New Version is coming soon! 🚨 [rnalyze v1.0.1]

**☑️ Support for Multiple Trimming Tools**

- Added support for Cutadapt alongside Trimmomatic, providing more flexibility and customization for trimming RNA-Seq data.
  
**☑️ Platform-Independent CPU Core Detection**
  
- The script now detects the number of CPU cores in a platform-independent way, supporting Linux, macOS, and Windows.
  
**☑️ Enhanced Argument Validation**

- Introduced a robust `validate_required_args` function to validate all input arguments, ensuring the script fails early with clear error messages if any required arguments are missing or invalid.
  
**☑️ Enhanced Reference Genome Handling**

- Added validation for reference genome file extensions (e.g., .fna, .fa, .fasta).
  
**☑️ Detailed Usage Documentation**

- The `-h` menu now provides a detailed, table-formatted guide with examples for different use cases.
  
- Clearer instructions for both beginners and advanced users.
  
### Overview
**RNAlyze** is an automated and scalable pipeline designed for RNA-Seq data processing, developed by **Mohamed Elsisi** and **Mohamed Elhwary**. It streamlines the entire RNA sequencing workflow, from raw data acquisition to feature quantification, with high efficiency and flexibility.

Now available on **Bioconda**, RNAlyze makes **Next-Generation Sequencing (NGS) expression analysis** easier and more automated than ever. Designed with bioinformaticians and researchers in mind, RNAlyze simplifies RNA-Seq data analysis with an automated, user-friendly pipeline, eliminating the need for complex command-line operations. Whether you're new to bioinformatics or an expert looking for an efficient workflow, RNAlyze ensures that anyone can easily process RNA-Seq data with minimal effort.

### 📢IMPORTANT NOTE: The pipeline only works in a conda environment with python version 3.10

### Key Features

**- Multi-Stage Pipeline:** Supports two processing modes: Full, and Alignment.

**- Flexible Alignment:** Choose from BWA, Bowtie2, and HISAT2 for read mapping and genome indexing.

**- Data Acquisition:** Handles both Download (via SRA) and Directory (pre-downloaded files) options.

**- Trimming & Quality Control:** The user can now choose from Trimmomatic and **Cutadapt** 🆕​ for quality trimming  and FastQC/MultiQC for quality assessment.

**- Reference Genome Handling:** Works with indexed/unindexed genomes from URLs or local paths.

**- Feature Quantification:** Utilizes FeatureCounts for generating count matrices.

**- Scalability:** Efficiently utilizes multi-threading for optimized performance.

**- Logging & Monitoring:** Detailed log outputs for error handling and reproducibility.

### Installation

**- To install:**

    conda install -c bioconda rnalyze 

or
    
    mamba install rnalyze

**To test if the package installed or not:**

    rnalyze -h
    
**- Usage:**

Run the pipeline with the desired options:

    rnalyze -p Full -t HISAT2 -d Download -s /path/to/SRR_Acc_List.txt -y SE \
        -r UnindexedURL -u http://example.com/ref.fna -T yes -a TruSeq3-SE \
        -l 5 -g Download -G http://example.com/annotation.gtf -i gene_id

### Command-Line Options

**General Parameters**

| Option              | Description                                                                                     | Inputs                                                                 |
|---------------------|-------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
| `-p, --pipeline`    | Specify pipeline type: Full, Alignment                                                          | `Full`, `Alignment`                                                   |
| `-t, --tool`        | Choose alignment tool: BWA, Bowtie2, HISAT2                                                    | `BWA`, `Bowtie2`, `HISAT2`                                            |
| `-d, --data`        | Data source: Download (via SRA) or Directory (local fastq.gz files)                             | `Download`, `Directory`                                               |
| `-s, --srr-path`    | Path to `SRR_Acc_List.txt` (required for Download mode)                                         | File path (e.g., `/path/to/SRR_Acc_List.txt`)                         |
| `-D, --data-dir`    | Directory containing raw sequencing files (required for Directory mode)                        | Directory path (e.g., `/path/to/data`)                                |
| `-y, --type`        | RNA-Seq data type: SE (Single-End) or PE (Paired-End)                                           | `SE`, `PE`                                                            |
| `-r, --ref`         | Reference genome format: UnindexedURL, IndexedURL, IndexedPath, UnindexedPath                   | `UnindexedURL`, `IndexedURL`, `IndexedPath`, `UnindexedPath`           |
| `-u, --ref-url`     | URL for downloading the reference genome (required if UnindexedURL or IndexedURL is selected)   | URL (e.g., `http://example.com/ref.fna`)                              |
| `-R, --ref-path`    | Path to the reference genome (required for local genome processing)                             | File path (e.g., `/path/to/ref_genome.fna`)                           |
| `-T, --trim`        | Enable trimming (yes or no) (default: no)                                                      | `yes`, `no`                                                           |
| `-x, --trim-tool`   | Selects trimming tool (if trimming is enabled).                                                | `Trimmomatic`, `Cutadapt`                                             |
| `-g, --gtf`         | GTF annotation source: Download or Path                                                        | `Download`, `Path`                                                    |
| `-G, --gtf-url`     | URL for downloading the GTF file (required if Download mode is selected)                        | URL (e.g., `http://example.com/annotation.gtf`)                       |
| `-P, --gtf-path`    | Local path to the GTF annotation file (required for Path mode)                                  | File path (e.g., `/path/to/annotation.gtf`)                           |
| `-i, --identifier`  | Gene identifier attribute in the GTF file (default: `gene_id`)                                 | (e.g., `gene_id`, `transcript_id`)                                    |
| `-h, --help`        | Show help message                                                                               | None                                                                  |

**Trimmomatic Parameters**

| Option              | Description                                                                                     | Inputs                                                                 |
|---------------------|-------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
| `-l, --leading`     | Trimmomatic LEADING threshold (default: 3)                                                     | Integer (e.g., `5`)                                                   |
| `-L, --trailing`    | Trimmomatic TRAILING threshold (default: 3)                                                    | Integer (e.g., `5`)                                                   |
| `-w, --sliding-window` | Trimmomatic SLIDINGWINDOW parameters (default: 4:25)                                        | String in format `window_size:quality_threshold` (e.g., `4:25`)       |
| `-m, --minlen`      | Trimmomatic MINLEN threshold (default: 36)                                                     | Integer (e.g., `36`)                                                  |
| `-a, --adapter-se`  | Adapter for Single-End data (required if --trim yes and --type SE)                              | Adapter name (e.g., `TruSeq3-SE`)                                     |
| `-A, --adapter-pe`  | Adapter for Paired-End data (required if --trim yes and --type PE)                              | Adapter name (e.g., `TruSeq3-PE`)                                     |

**🆕​ Cutadapt Parameters**

| Option              | Description                                                                                     | Inputs                                                                 |
|---------------------|-------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
| `-a, --adapter-se`  | Adapter sequence for Single-End data.                                                           | Adapter sequence                                                       |
| `-A, --adapter-pe`  | Adapter sequence for Paired-End data.                                                           | Adapter sequence                                                       |
| `-g, --front-se`    | 5' adapter sequence for Single-End data.                                                        | Adapter sequence                                                       |
| `-G, --front-pe`    | 5' adapter sequence for Paired-End data.                                                        | Adapter sequence                                                       |
| `-q, --quality-cutoff` | Trim bases below quality threshold.                                                          | Integer                                                                |
| `-M, --min-read-length` | Minimum read length after trimming.                                                         | Integer                                                                |

### How The Pipelines Works


**1. Full Pipeline**

- The Full pipeline automates the entire RNA-Seq workflow, from raw data acquisition to feature quantification. Here's how it works step-by-step:

*1. Data Acquisition:*

- Downloads raw sequencing data via SRA or processes pre-downloaded files from a directory.

*2. Quality Control:*

- Runs FastQC on raw reads to generate initial quality reports.
- Trims low-quality bases and adapters using Trimmomatic or Cutadapt based in your parameters.
- Re-runs FastQC on trimmed reads to ensure improved quality.

*3. Reference Genome Handling:*

- Downloads the reference genome from a URL if `UnindexedURL` or `IndexedURL` is specified.
- Uses a local reference genome if `IndexedPath` or `UnindexedPath` is provided.
- Indexes the genome if it is unindexed.

*4. Alignment:*

- Maps trimmed reads to the reference genome using the chosen alignment tool (BWA, Bowtie2, or HISAT2).
- Converts SAM files to BAM format and sorts/indexes them for downstream analysis.

*5. Feature Quantification:*

- Uses FeatureCounts to quantify gene expression levels based on the alignment results.
- Generates count matrices for downstream differential expression analysis.

*6. Output Generation:*

- Produces comprehensive output directories, including FastQC reports, alignment files, and feature count tables.

**2.Alignment Pipeline**

The Alignment pipeline focuses solely on mapping reads to the reference genome. It stops after generating sorted and indexed BAM files.

*1. Data Acquisition:*

- Downloads raw sequencing data via SRA or processes pre-downloaded files from a directory.

*2. Quality Control:*

- Runs FastQC on raw reads to generate initial quality reports.
- Trims low-quality bases and adapters using Trimmomatic or Cutadapt based on your parameters.

*3. Reference Genome Handling:*

- Downloads the reference genome from a URL if `UnindexedURL` or `IndexedURL` is specified.
- Uses a local reference genome if `IndexedPath` or `UnindexedPath` is provided.
- Indexes the genome if it is unindexed.

*4. Alignment:*

- Maps trimmed reads to the reference genome using the chosen alignment tool (BWA, Bowtie2, or HISAT2).
- Converts SAM files to BAM format and sorts/indexes them.

*5. Output Generation:*

- Produces alignment files (SAM, BAM, sorted BAM, and index files) for further downstream analysis.

### Example Workflows

*- Full Pipeline with Hisat2 & Trimmomatic (Single-End Data)*

      rnalyze -p Full -t HISAT2 -d Download -s /path/to/SRR_Acc_List.txt -y SE -r UnindexedURL -u http://example.com/ref.fna -T yes -x Trimmomatic -a TruSeq3-SE -l 5 -f Download -F http://example.com/annotation.gtf -i gene_id

*- Alignment Pipeline with BWA & Cutadapt (Paired-End Data)*

      rnalyze -p Alignment -t BWA -d Directory -D /path/to/data -y PE -r IndexedPath -R /path/to/indexed_genome -T yes -x Cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20 -M 30

*- Full Pipeline with Cutadapt (Paired-End Data)*

      rnalyze -p Full -t HISAT2 -d Directory -D /path/to/data -y PE -r UnindexedPath -R /path/to/ref_genome.fna -T yes -x Cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g AGATCGGAAGAGC -G AGATCGGAAGAGC -q 20 -M 30 -f Path -P /path/to/annotation.gtf -i gene_id

*- Alignment Pipeline with Cutadapt (Paired-End Data)*

       rnalyze -p Alignment -t Bowtie2 -d Download -s /path/to/SRR_Acc_List.txt -y SE -r IndexedURL -u http://example.com/indexed_genome.tar.gz -T yes -x Trimmomatic -a TruSeq3-SE -l 5

### Output Structure

- RNAlyze is designed to save your work in a highly organized and optimized directory structure. This ensures that all outputs are systematically stored, making it easy to locate and analyze results. Below is a detailed explanation of how the directories are created and used.

- RNAlyze automatically creates a main working directory named `RNAlyze-<timestamp>` in your current working directory. This directory contains all the subdirectories for storing pipeline outputs. The structure is as follows:

1. `RNAlyze-<timestamp>/`: The main working directory, named with a timestamp to ensure uniqueness.

    1.1. `Data/`: Stores raw and processed sequencing files (.fastq.gz files).

    1.2. `Fastqc_Reports/`: Contains FastQC reports for quality control.

      1.2.1. `Before_Trimming/`: FastQC reports for raw data (created if trim=yes).

     1.2.2. `After_Trimming/`: FastQC reports for trimmed data (created if trim=yes).

    1.3. `MultiQC_Reports/`: Contains MultiQC reports summarizing FastQC results.

     1.3.1. `Before_Trimming/`: MultiQC reports for raw data (created if trim=yes).

     1.3.2. `After_Trimming/`: MultiQC reports for trimmed data (created if trim=yes).

    1.4. `Mapping/`: Stores alignment outputs (.sam, .bam, sorted .bam, and index files).

    1.5. `Ref_Genome/`: Contains reference genome files (.fna, index files).

    1.6. `Annotation/`: Stores annotation files (.gtf files).

    1.7. `Featurecount/`: Contains gene expression count tables (.counts.txt, merged and preprocessed count matrices).

### Authors

**- Mohamed Elsisi:**

*M.elsisi@nu.edu.eg* 

**- Mohamed Elhwary:**

*M.elhwary@nu.edu.eg*

### License

**RNAlyze is released under the MIT License.**

- For issues and feature requests, please submit a report on GitHub or contact the authors directly.


