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
