# RNAlyze
### Overview
**RNAlyze** is an automated and scalable pipeline designed for RNA-Seq data processing, developed by **Mohamed Elsisi** and **Mohamed Elhwary**. It streamlines the entire RNA sequencing workflow, from raw data acquisition to feature quantification, with high efficiency and flexibility.

Now available on **Bioconda**, RNAlyze makes **Next-Generation Sequencing (NGS) expression analysis** easier and more automated than ever. Designed with bioinformaticians and researchers in mind, RNAlyze simplifies RNA-Seq data analysis with an automated, user-friendly pipeline, eliminating the need for complex command-line operations. Whether you're new to bioinformatics or an expert looking for an efficient workflow, RNAlyze ensures that anyone can easily process RNA-Seq data with minimal effort.

### Key Features

**- Multi-Stage Pipeline:** Supports three processing modes: Full, Alignment, and Alignment&Featurecount.

**- Flexible Alignment:** Choose from BWA, Bowtie2, and HISAT2 for read mapping.

**- Data Acquisition:** Handles both Download (via SRA) and Directory (pre-downloaded files) options.

**- Trimming & Quality Control:** Uses Trimmomatic for quality trimming and FastQC/MultiQC for quality assessment.

**- Reference Genome Handling:** Works with indexed/unindexed genomes from URLs or local paths.

**- Feature Quantification:** Utilizes FeatureCounts for generating count matrices.

**- Scalability:** Efficiently utilizes multi-threading for optimized performance.

**- Logging & Monitoring:** Detailed log outputs for error handling and reproducibility.

### Installation

**- To install:**

    conda install -c bioconda rnalyze 

**- Usage:**

Run the pipeline with the desired options:

    rnalyze -p Full -t HISAT2 -d Download -s /path/to/SRR_Acc_List.txt -y SE \
        -r UnindexedURL -u http://example.com/ref.fna -T yes -a TruSeq3-SE \
        -l 5 -g Download -G http://example.com/annotation.gtf -i gene_id

### Command-Line Options

| Option        | Description |
|--------------|------------|
| `-p, --pipeline` | Specify pipeline type: Full, Alignment, Alignment&Featurecount |
| `-t, --tool` | Choose alignment tool: BWA, Bowtie2, HISAT2 |
| `-d, --data` | Data source: Download (via SRA) or Directory (local fastq.gz files) |
| `-s, --srr-path` | Path to `SRR_Acc_List.txt` (required for Download mode) |
| `-D, --data-dir` | Directory containing raw sequencing files (required for Directory mode) |
| `-y, --type` | RNA-Seq data type: SE (Single-End) or PE (Paired-End) |
| `-r, --ref` | Reference genome format: UnindexedURL, IndexedURL, IndexedPath, UnindexedPath |
| `-u, --ref-url` | URL for downloading the reference genome (required if UnindexedURL or IndexedURL is selected) |
| `-R, --ref-path` | Path to the reference genome (required for local genome processing) |
| `-T, --trim` | Enable trimming (yes or no) |
| `-g, --gtf` | GTF annotation source: Download or Path |
| `-G, --gtf-url` | URL for downloading the GTF file (required if Download mode is selected) |
| `-P, --gtf-path` | Local path to the GTF annotation file (required for Path mode) |
| `-i, --identifier` | Gene identifier attribute in the GTF file (default: `gene_id`) |
| `-h, --help` | Show help message |

### Example Workflows

*- Running Full Pipeline with HISAT2*

      rnalyze -p Full -t HISAT2 -d Directory -D /data -y PE \
        -r IndexedPath -R /genomes/hg38 \
        -g Path -P /annotations/gencode.gtf -i gene_id

*- Running Alignment Only with Bowtie2*

      rnalyze -p Alignment -t Bowtie2 -d Download -s /path/to/SRR_Acc_List.txt -y SE \
        -r UnindexedURL -u http://example.com/ref.fna

### Output Structure

**1- Fastqc_Reports/** → Raw and trimmed data quality control reports

**2- MultiQC_Reports/** → Combined reports summarizing quality control results

**3- Data/** → Processed raw and trimmed sequencing reads

**4- Mapping/** → Alignment output files (SAM, BAM, sorted BAM, index files)

**5- Annotation/** → Reference genome annotation files (GTF)

**6- Featurecount/** → Gene expression count tables (merged & preprocessed)

### Performance Optimization

**1- Multi-Threading:** The pipeline utilizes the maximum available CPU threads minus one for optimal performance.

**2- Automatic Cleanup:** Removes unnecessary files based on the selected pipeline mode to conserve storage.

**3- Error Handling:** Built-in validation and logging for missing or incorrect inputs.

### Authors

**- Mohamed Elsisi:**

*M.elsisi@nu.edu.eg* 

**- Mohamed Elhwary:**

*M.elhwary@nu.edu.eg*

### License

**RNAlyze is released under the MIT License.**

- For issues and feature requests, please submit a report on GitHub or contact the authors directly.


