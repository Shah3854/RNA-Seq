#!/bin/bash

inputDir=""
dir_qc="/fastqc"
dir_fastp="/filtered_qc_report"
file_ref="hg38_gencode"
log_file="log_$(date +'%Y%m%d%H%M%S').txt"
terminal_log="terminal_log_$(date +'%Y%m%d%H%M%S').txt"

Dir_map="/SAM"
[ -d $Dir_map ] || mkdir $Dir_map

Dir_counts="/Counts"
[ -d $Dir_counts ] || mkdir $Dir_counts

# Function to log messages
log() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') $1" >> "$log_file"
}

# Function to log terminal output
log_terminal() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') $1" >> "$terminal_log"
}

# Start logging
log "Script started."

# Redirect terminal output to log file
exec > >(tee -a "$terminal_log") 2>&1

for f in "$inputDir"/*.fastq.gz
do 
    if [ "$f" != "$inputDir/ls" ]
    then
        S=$(basename "$f")
        variable=$(echo "$S" | awk -F"_R" '{print $1;}')
        var1="$var1 $variable"
    fi
done

unqpre=$(echo "$var1" | xargs -n1 | sort -u | xargs)
echo "$unqpre"
log "Total number of Libraries are: $(echo "$unqpre" | wc -w)"

for f in $unqpre
do
    [ -d "$dir_qc" ] || mkdir -p "$dir_qc"

    filep1="$inputDir/${f}_R1.fastq.gz"
    filep2="$inputDir/${f}_R2.fastq.gz"

    fastqc -t 18 "$filep1" "$filep2" -o "$dir_qc"
    log "FastQC completed for $f."

    # The multi-QC will be at the end of analysis to avoid repetition
    [ -d "$dir_fastp" ] || mkdir "$dir_fastp"

    fastp -t 18 -i "$filep1" -I "$filep2" --disable_length_filtering --qualified_quality_phred 20 -o "${dir_fastp}/${f}_1_filtered.fastq.gz" -O "${dir_fastp}/${f}_2_filtered.fastq.gz" --html "${dir_fastp}/${f}.html"
    mv fastp.json "${dir_fastp}/${f}.json"
    log "Fastp completed for $f."

    
    bwa mem -t 18  $file_ref/GRCh38.p14.genome.fa ${dir_fastp}/${f}_1_filtered.fastq.gz ${dir_fastp}/${f}_2_filtered.fastq.gz  > "${Dir_map}/${f}.sam"
    log "BWA MEM completed for $f."
    
    htseq-count "${Dir_map}/${f}.sam" $file_ref/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf > "${Dir_counts}/${f}_counts.tsv"

done   
    
# End logging
log "Script completed."
log_terminal "Terminal output captured."
    

