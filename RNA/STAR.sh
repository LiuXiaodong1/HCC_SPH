

STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir /project/liuxd/liver/RNA-analysis/output/STAR/genomeDir \
--genomeFastaFiles /project/liuxd/liver/RNA-seq-data/ref/genome.fa \
--sjdbGTFfile /project/liuxd/liver/RNA-seq-data/ref/genome.gtf \
--sjdbOverhang 149


STAR --runThreadN 4 \
--genomeDir /project/liuxd/liver/RNA-analysis/output/STAR/genomeDir \
--readFilesIn /project/liuxd/liver/RNA-seq-data/CleanData/${1}_1.clean.fq.gz /project/liuxd/liver/RNA-seq-data/CleanData/${1}_2.clean.fq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 4 \
--outFileNamePrefix /project/liuxd/liver/RNA-analysis/output/STAR/$1 \
--quantMode TranscriptomeSAM GeneCounts 
