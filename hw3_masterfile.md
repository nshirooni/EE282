# Homework 3

## Summarize Genome Assembly

To complete the following task of summarizing the genome assembly the following code which was used will be located in the **code** folder of this repository. The following script "**hw3_genome_summary.sh**" will complete the task. The data and outputs are located in the **genome** folder with the path **/genome/fasta**. The complete code is shown below:

```{bash}
srun -A class_ee282 --pty bash -i
conda activate ee282


mkdir -p ~/github/EE282/data/genome/fasta && \
wget https://ftp.flybase.net/releases/current/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz -P ~/github/EE282/data/genome/fasta && \
wget https://ftp.flybase.net/releases/current/dmel_r6.60/fasta/md5sum.txt -P ~/github/EE282/data/genome/fasta && \
cd ~/github/EE282/data/genome/fasta && \
grep "dmel-all-chromosome-r6.60.fasta.gz" md5sum.txt | md5sum -c > file_integrity.txt && \
cat file_integrity.txt && \
zcat dmel-all-chromosome-r6.60.fasta.gz | faSize -tab stdin \
| sort -rnk 2,2 \
| awk '/baseCount|seqCount|nBaseCount/' \
| tee genome_summary.txt 

```

First a node is claimed on the HPC3 to ensure no computation is done on the login node followed by the activation of the proper conda enviornment to ensure task completion

```{bash}
srun -A class_ee282 --pty bash -i
conda activate ee282
```

The next step is to create the proper directory for the fasta files and use wget to download both the fasta file for the chromosome and the .txt file md5sum for file integrity.

```{bash}
mkdir -p ~/github/EE282/data/genome/fasta && \
wget https://ftp.flybase.net/releases/current/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz -P ~/github/EE282/data/genome/fasta && \
wget https://ftp.flybase.net/releases/current/dmel_r6.60/fasta/md5sum.txt -P ~/github/EE282/data/genome/fasta 
```

A grep function is then applied only to check the chromosome fasta in the md5sum.txt file and the output shows that the file is good.

```{bash}
grep "dmel-all-chromosome-r6.60.fasta.gz" md5sum.txt | md5sum -c > file_integrity.txt && \
cat file_integrity.txt
# Successful integrity check
dmel-all-chromosome-r6.60.fasta.gz: OK

```

Lastly faSize is used with the -tab option and then edited with awk to only include the 3 features of interest to summarize the genome.

```{bash}
zcat dmel-all-chromosome-r6.60.fasta.gz | faSize -tab stdin \
| sort -rnk 2,2 \
| awk '/baseCount|seqCount|nBaseCount/' \
| tee genome_summary.txt 

# Successful output showing the 3 important features 
baseCount       143726002
nBaseCount      1152978
seqCount        1870

```

## Summarize an Annotation File

To complete the following task of summarizing an annotation file the following code which was used will be located in the **code** folder of this repository. The following script "**hw3_annotation_summary.sh**" will complete the task. The data and outputs are located in the **genome** folder with the path **/genome/gtf**. The complete code is shown below:

```{bash}
srun -A class_ee282 --pty bash -i
conda activate ee282

mkdir -p ~/github/EE282/data/genome/gtf && \
wget https://ftp.flybase.net/releases/current/dmel_r6.60/gtf/dmel-all-r6.60.gtf.gz -P ~/github/EE282/data/genome/gtf && \
wget https://ftp.flybase.net/releases/current/dmel_r6.60/gtf/md5sum.txt -P ~/github/EE282/data/genome/gtf && \
cd ~/github/EE282/data/genome/gtf && \
md5sum -c md5sum.txt > file_integrity.txt && \
cat file_integrity.txt && \
zcat dmel-all-r6.60.gtf.gz | head \
| tee >(bioawk -c gff '{print $3}' dmel-all-r6.60.gtf.gz | sort | uniq -c | sort -nr /dev/stdin \
| tee feature_summary.txt) \
| tee >(bioawk -c gff '$3=="gene" {print $1}' dmel-all-r6.60.gtf.gz | grep -E "^.*[2-4][LR]|^.*X|^.*Y|^4" | sort | uniq -c | sort -nr /dev/stdin \
| tee chromosome_summary.txt)

```

The first two steps will be similar to the first problem but in this case we are downloading the .gtf file and will also not need a grep function to filter through the md5sum.txt file since there is only the .gtf file.

```{bash}
srun -A class_ee282 --pty bash -i
conda activate ee282

mkdir -p ~/github/EE282/data/genome/gtf && \
wget https://ftp.flybase.net/releases/current/dmel_r6.60/gtf/dmel-all-r6.60.gtf.gz -P ~/github/EE282/data/genome/gtf && \
wget https://ftp.flybase.net/releases/current/dmel_r6.60/gtf/md5sum.txt -P ~/github/EE282/data/genome/gtf && \
cd ~/github/EE282/data/genome/gtf && \
md5sum -c md5sum.txt > file_integrity.txt && \
cat file_integrity.txt

#Correct output for integrity 
dmel-all-r6.60.gtf.gz: OK

```

To summarize the annotation file we first get a summary of the features. Bioawk is used to print the 3rd field which is correct field of interest for features, the data is then sorted to print most common to least common feature. Lastly to summarize the total number of genes per chromosome arm we index the gene and then print out the chromosome names and then apply a grep function to only include the chromosome arms of interest. The output is then printed.

```{bash}
zcat dmel-all-r6.60.gtf.gz | head \
| tee >(bioawk -c gff '{print $3}' dmel-all-r6.60.gtf.gz | sort | uniq -c | sort -nr /dev/stdin \
| tee feature_summary.txt) \
| tee >(bioawk -c gff '$3=="gene" {print $1}' dmel-all-r6.60.gtf.gz | grep -E "^.*[2-4][LR]|^.*X|^.*Y|^4" | sort | uniq -c | sort -nr /dev/stdin \
| tee chromosome_summary.txt)

#Output according to text files
 190038 exon
 163253 CDS
  46806 5UTR
  33741 3UTR
  30888 start_codon
  30828 stop_codon
  30802 mRNA
  17872 gene
   3059 ncRNA
    485 miRNA
    365 pseudogene
    312 tRNA
    270 snoRNA
    262 pre_miRNA
    115 rRNA
     32 snRNA
     
   4226 3R
   3649 2R
   3508 2L
   3481 3L
   2704 X
    114 4
    113 Y

```
