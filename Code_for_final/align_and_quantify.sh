#!/bin/bash
#SBATCH -A MFROSE_LAB     
#SBATCH -p standard              
#SBATCH --nodes=1             
#SBATCH --ntasks=1           
#SBATCH --cpus-per-task=32     
#SBATCH --mem=92G

conda activate star

# Input directory and sample name prefix
input_dir=$1
sample_prefix=$2

# Define directories
GENOME_DIR="/pub/nshiroon/genomeDir"
RSEM_REF="/pub/nshiroon/rsem/ref"
RSEM_RESULTS_DIR="/dfs6/pub/nshiroon/RSEM_results"

# Create output directory if it doesn't exist
mkdir -p "$RSEM_RESULTS_DIR"

gzip -d $input_dir/CN*

# Collect all `_R1` and `_R2` reads into variables
read_files_r1=($(ls "$input_dir"/*_R1*.fastq))
read_files_r2=($(ls "$input_dir"/*_R2*.fastq))

# Generate sample name
sample_name="${sample_prefix}"

# Run STAR alignment for all paired reads
STAR --genomeDir $GENOME_DIR \
     --readFilesIn ${read_files_r1[@]} ${read_files_r2[@]} \
     --runThreadN 32 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped None \
     --outFilterMismatchNmax 3 \
     --outFilterMultimapNmax 1 \
     --outSAMattributes All \
     --quantMode TranscriptomeSAM GeneCounts \
     --outFileNamePrefix "$input_dir/${sample_name}_"
   
# Create RSEM directory and move BAM file
mkdir -p "$input_dir/rsem"
mv "$input_dir/${sample_name}_Aligned.toTranscriptome.out.bam" "$input_dir/rsem/"
cd "$input_dir/rsem"

# Run RSEM
rsem-calculate-expression --paired-end \
                          --alignments \
                          -p 32 \
                          --no-bam-output \
                          "${sample_name}_Aligned.toTranscriptome.out.bam" \
                          $RSEM_REF \
                          $sample_prefix

# Move RSEM results to results directory
cp ${sample_prefix}* $RSEM_RESULTS_DIR/
