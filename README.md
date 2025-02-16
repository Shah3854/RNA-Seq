# RNA-Seq Analysis Pipeline

This repository contains a Bash script for processing RNA-Seq data. The script performs quality control, filtering, alignment, and gene quantification using various bioinformatics tools.

## Features
- Quality control using **FastQC**
- Read filtering and trimming using **fastp**
- Alignment with **BWA-MEM**
- Gene quantification using **HTSeq-Count**
- Logs terminal output and processing status

## Prerequisites
Ensure the following dependencies are installed:
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [fastp](https://github.com/OpenGene/fastp)
- [BWA](http://bio-bwa.sourceforge.net/)
- [HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/)
- Bash shell environment (Linux/macOS)

## Installation
Clone the repository:
```bash
git clone https://github.com/yourusername/rna-seq-pipeline.git
cd rna-seq-pipeline
```

## Usage
Modify the `inputDir` and `file_ref` variables in the script as needed. Then, run the script:
```bash
chmod +x script.sh
./script.sh
```

### Input Data
Place paired-end FASTQ files (`*_R1.fastq.gz` and `*_R2.fastq.gz`) in the specified `inputDir`.

### Output Files
- `fastqc/`: Quality control reports
- `filtered_qc_report/`: Filtered FASTQ files and reports
- `SAM/`: Alignment output in SAM format
- `Counts/`: Gene expression counts
- `log_*.txt`: Processing logs
- `terminal_log_*.txt`: Terminal output logs

