#rsem-prepare-reference --gtf /project/liuxd/liver/RNA-seq-data/ref/genome.gtf \
					     --STAR \
   				     	/project/liuxd/liver/RNA-seq-data/ref/genome.fa \
					/project/liuxd/liver/RNA-analysis/output/RSEM/ref/ref
rsem-calculate-expression -p 4 --paired-end \
					--bam \
					--estimate-rspd \
					--append-names \
					--output-genome-bam \
					/project/liuxd/liver/RNA-analysis/output/STAR/${1}Aligned.toTranscriptome.out.bam \
					/project/liuxd/liver/RNA-analysis/output/RSEM/ref/ref \
					/project/liuxd/liver/RNA-analysis/output/RSEM/exp/$1
