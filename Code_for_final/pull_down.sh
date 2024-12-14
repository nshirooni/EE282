#First step is to transfer FASTQ's from lab computer onto my dfs6 storage in the hp3 using this code:

#FROM LAB COMPUTER TO HPC3
scp -r FASTQs_PooledRNAseq_Matt-Max-Alicia-06-2015/ nshiroon@hpc3.rcic.uci.edu:/pub/nshiroon

#next I create a directory for ref genome and pull it down from CRSP LAB storage
mkdir /pub/nshiroon/ref_genome
cp /share/crsp/lab/mfrose/share/RNAseq/DATA/ref/Mus_musculus.GRCm39.107.filtered_eGFP_tdTomato.gtf /pub/nshiroon/ref_genome

#Lastly pull down the .fa file from the lab computer
#FROM LAB COMPUTER TO HPC3
scp Mus_musculus.GRCm39.dna.primary_assembly_eGFP_tdTomato.fa nshiroon@hpc3.rcic.uci.edu:/pub/nshiroon/ref_genome/