STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir /pub/nshiroon/genomeDir \
--genomeFastaFiles /pub/nshiroon/ref_genome/Mus_musculus.GRCm39.dna.primary_assembly_eGFP_tdTomato.fa \
--sjdbGTFfile /pub/nshiroon/ref_genome/Mus_musculus.GRCm39.107.filtered_eGFP_tdTomato.gtf

rsem-prepare-reference --gtf /pub/nshiroon/ref_genome/Mus_musculus.GRCm39.107.filtered_eGFP_tdTomato.gtf \
                       --star \
                        -p 32 \
                        /pub/nshiroon/ref_genome/Mus_musculus.GRCm39.dna.primary_assembly_eGFP_tdTomato.fa \
                        /pub/nshiroon/rsem/ref 