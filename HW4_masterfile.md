# Homework 4

## Part 1: Summarize partitions of a genome

### Filtering by sequence size and calculating total number of Nucleotides, Ns, and Sequences

The first step is to pull down the most current FlyBase genome and then download the fasta file for all chromosomes as shown below, this step is analogous to HW3:

```{bash}
srun -A class_ee282 --pty bash -i
conda activate ee282

mkdir -p ~/genome/fasta && \
wget https://ftp.flybase.net/releases/current/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz -P ~/genome/fasta && \
wget https://ftp.flybase.net/releases/current/dmel_r6.60/fasta/md5sum.txt -P ~/genome/fasta && \
cd ~/genome/fasta && \
grep "dmel-all-chromosome-r6.60.fasta.gz" md5sum.txt | md5sum -c > file_integrity.txt && \
cat file_integrity.txt
```

The next step is to filter the fasta file using bioawk for sequences larger and less than or equal to 100kb and then gzip the files:

```{bash}
bioawk -c fastx '(length($seq) > 100000) {print ">" $name ORS $seq}' dmel-all-chromosome-r6.60.fasta.gz  > dmel_greater_100.fasta 

bioawk -c fastx '(length($seq) <= 100000) {print ">" $name ORS $seq}' dmel-all-chromosome-r6.60.fasta.gz > dmel_less_100.fasta 

gzip dmel_*
```

Now we have the filtered fasta files which can be used downstream for the required tasks. The next step is to calculate a genome summary which can be achieved via faSize and using tab and awk for better formatting, the following below is the code to achieve this as well as the output:

```{bash}
zcat dmel_greater_100.fasta.gz | faSize -tab stdin \
| sort -rnk 2,2 \
| awk '/baseCount|seqCount|nBaseCount/' \
| tee genome_summary_greater_100.txt 

baseCount       137547960
nBaseCount      490385
seqCount        7


zcat dmel_less_100.fasta.gz | faSize -tab stdin \
| sort -rnk 2,2 \
| awk '/baseCount|seqCount|nBaseCount/' \
| tee genome_summary_less_100.txt   

baseCount       6178042
nBaseCount      662593
seqCount        1863


```

### Plotting sequence length distribution

The first step is to save these sequence lengths to a .txt file and to use bioawk to grab the sequence lengths:

```{bash}
bioawk -c fastx '{print $name, length($seq)}' dmel_greater_100.fasta.gz > sequence_lengths_greater.txt 

bioawk -c fastx '{print $name, length($seq)}' dmel_less_100.fasta.gz > sequence_lengths_less.txt
```

For plotting it is now required to load up R to make use of ggplot2 for a histogram distribution, first for sequences greater than 100kb:

```{r}
data <- read.table("~/genome/fasta/sequence_lengths_greater.txt", header = FALSE)
lengths <- data$V2

library(ggplot2)
library(scales)

ggplot(data.frame(lengths = lengths), aes(x = lengths)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  scale_x_log10(
    labels = scientific_format(), 
  ) +
  labs(
    title = "Sequence Length Distribution for greater than 100kb",
    x = "Log(Sequence Length)",
    y = "Frequency"
  ) 
```

![](https://hpc3.rcic.uci.edu/biojhub4/user/nshiroon/rstudio/graphics/94fa675b-9c7f-49da-9de5-89b0fe431520.png)

Now for smaller than or equal to 100kb

```{r}
data <- read.table("~/genome/fasta/sequence_lengths_less.txt", header = FALSE)
lengths <- data$V2

library(ggplot2)
library(scales)

ggplot(data.frame(lengths = lengths), aes(x = lengths)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  scale_x_log10(
    labels = scientific_format(), 
  ) +
  labs(
    title = "Sequence Length Distribution for less than or equal to 100kb",
    x = "Log(Sequence Length)",
    y = "Frequency"
  ) 
```

![](https://hpc3.rcic.uci.edu/biojhub4/user/nshiroon/rstudio/graphics/2beeac0c-e9e2-4e09-870f-c6bab80f9733.png)

### Plotting GC content distribution

The next step is to use bioawk to obtain the GC content and then make a similar histogram in R

```{bash}
bioawk -c fastx '{ print $name; print gc($seq) }' dmel_greater_100.fasta.gz > gc_content_greater_100.txt

bioawk -c fastx '{ print $name; print gc($seq) }' dmel_less_100.fasta.gz > gc_content_less_100.txt
```

Now plot in R first starting with sequences greater than 100kb:

```{r}
# Read the GC content data
gc_data <- read.table("gc_content_greater_100.txt")

library(ggplot2)

# Create histogram
ggplot(gc_data, aes(x = V2)) +
  geom_histogram(
    bins = 50, 
    fill = "blue", 
    color = "black"
  ) +
  labs(
    title = "GC Content Distribution for sequences greater than 100kb",
    x = "GC Content", 
    y = "Frequency"
  ) 
```

![](https://hpc3.rcic.uci.edu/biojhub4/user/nshiroon/rstudio/graphics/fb509ef6-6890-4ba6-9e22-7a16f80bcd40.png)

Now for sequences less than or equal to 100kb

```{bash}
gc_data <- read.table("gc_content_less_100.txt")

library(ggplot2)

# Create histogram
ggplot(gc_data, aes(x = V2)) +
  geom_histogram(
    bins = 50, 
    fill = "blue", 
    color = "black"
  ) +
  labs(
    title = "GC Content Distribution for sequences less than or equal to 100kb",
    x = "GC Content", 
    y = "Frequency"
  ) 

```

![](https://hpc3.rcic.uci.edu/biojhub4/user/nshiroon/rstudio/graphics/24406a05-ba29-4e1f-8ff8-50267ff3a87a.png)

### CDF plots

The last step for part 1 is to get the cumulative sequence sizes and plot them using the plotCDF2 function given in the class. First we start with sizes greater than 100kb:

```{bash}
bioawk -c fastx '{ print length($seq) }' dmel_greater_100.fasta.gz \
  | sort -nr \
  | awk 'BEGIN { print "Assembly\tLength" } { print "dmel_greater\t" $1 }' \
  > cumulative_seq_sizes_greater_100.length &
  plotCDF2 cumulative_seq_sizes_greater_100.length plot_cdf_greater_100.png
```

![](https://hpc3.rcic.uci.edu/biojhub4/user/nshiroon/rstudio/files/genome/fasta/plot_cdf_greater_100.png)

Then lastly for sizes equal to or less than 100kb

```{bash}
bioawk -c fastx '{ print length($seq) }' dmel_less_100.fasta.gz \
  | sort -nr \
  | awk 'BEGIN { print "Assembly\tLength" } { print "Dmel_Less\t" $1 }' \
  > cumulative_seq_sizes_less_100.length &
  plotCDF2 cumulative_seq_sizes_less_100.length plot_cdf_less_100.png
```

![](https://hpc3.rcic.uci.edu/biojhub4/user/nshiroon/rstudio/files/genome/fasta/plot_cdf_less_100.png)

## Part 2: Genome Assembly

### Assemble a genome using PacBio Hifi Reads

This step involves grabbing fasta files for Hifi reads and then using hifiasm to assemble our genome into a contig file then using awk to convert the resulting .gfa contig file to a .fa contig file for further downstream analysis. We ensure atleast 16 cores are used to speed up the process ran by hifiasm.

```{bash}
HIFI_READS=/data/class/ee282/public/ISO1_Hifi_AdaptorRem.40X.fasta.gz

mkdir HiFi_assembly

cd ~/HiFi_assembly

srun -A class_ee282 -c 16 --pty bash -i

conda activate ee282

hifiasm -o iso1.asm -t16 $HIFI_READS 2> asm.log

awk '/^S/{print ">"$2"\n"$3}' iso1.asm.bp.p_ctg.gfa | fold > iso1.asm.bp.p_ctg.fa
```

### Assembly Assessment

#### Calculating N50

Using the following code below the N50 can be calculated, this was tested on the r6 contig and scaffold N50 since these sizes are known to ensure the function worked as expected:

```{bash}
bioawk -c fastx '{print length($seq)}' iso1.asm.bp.p_ctg.fa | \
sort -rn | \
awk 'BEGIN {total=0} {total+=$1; lengths[NR]=$1} 
     END {
         half=total/2; sum=0;
         for (i=1; i<=NR; i++) {
             sum+=lengths[i];
             if (sum>=half) {
                 print "N50: " lengths[i];
                 break;
             }
         }
     }'

N50: 21715751

```

For the N50 we get about 21.7 Mb on our assembly and when comparing this to the contig of r6 that is 21.5 Mb

#### Comparing assemblies using plotCDF2

This step is analogous to the FIFO step that we had done in class, adapting that code for the dmel genome we downloaded earlier as well as including our assembly we get the following code:

```{bash}
dmel=~/genome/fasta/dmel-all-chromosome-r6.60.fasta.gz
iso1=~/HiFi_assembly/iso1.asm.bp.p_ctg.fa
createProject fifos ~/fifo
cd ~/fifo/fifos

mkfifo tmp/{dmelscaff,dmelctg,iso1ctg}_fifo


bioawk -c fastx ' { print length($seq) } ' $dmel \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nDmel_Scaff\t0" } { print "Dmel_Scaff\t" $1 } ' \
> tmp/dmelscaff_fifo \
| faSplitByN $dmel /dev/stdout 10 \
| bioawk -c fastx ' { print length($seq) } ' \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nDmel_Ctg\t0" } { print "Dmel_Ctg\t" $1 } ' \
> tmp/dmelctg_fifo &

bioawk -c fastx ' { print length($seq) } ' $iso1 \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nIso1_Ctg\t0" } { print "Iso1_Ctg\t" $1 } ' \
> tmp/iso1ctg_fifo &

plotCDF2 tmp/{dmelscaff,dmelctg,iso1ctg}_fifo ~/CDF.png

rm tmp/{dmelscaff,dmelctg,iso1ctg}_fifo
```

From that we get the following plot

![](https://hpc3.rcic.uci.edu/biojhub4/user/nshiroon/rstudio/files/CDF.png)

#### Calculating BUSCO scores

To accomplish this we used compleasm and set the mode to run as BUSCO, we also allow 16 threads since the computation can be heavy, this will be first done for our contig assembly and then compared to the contig from FlyBase as well:

```{bash}
conda deactivate

conda activate compleasm

compleasm run -a iso1.asm.bp.p_ctg.fa -o /pub/nshiroon/compleasm_iso1 -l diptera_odb10 -m busco -t 16

## lineage: diptera_odb10
S:99.63%, 3273
D:0.24%, 8
F:0.00%, 0
I:0.00%, 0
M:0.12%, 4
N:3285


dmel=~/genome/fasta/dmel-all-chromosome-r6.60.fasta.gz
faSplitByN $dmel ~/genome/fasta/dmel_contig.fasta 10

cd ~/genome/fasta

compleasm run -a dmel_contig.fasta -o /pub/nshiroon/compleasm_dmel -l diptera_odb10 -m busco -t 16

## lineage: diptera_odb10
S:99.70%, 3275
D:0.30%, 10
F:0.00%, 0
I:0.00%, 0
M:0.00%, 0
N:3285


```
