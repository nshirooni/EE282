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

