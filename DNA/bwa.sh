module load bwa/0.7.17
module load samtools/1.9
module load gatk/4.1.3.0


inputpath="/project/liuxd/liver/WES-analysis/data/"
filepath="/project/liuxd/liver/WES-analysis/input/"
outputpath="/project/liuxd/liver/WES-analysis/output/bwa/"

bwa mem -t 4 -R "@RG\\tID:$(zcat ${inputpath}${2}|head -n 1|cut -d":" -f1|sed "s/@//")\\tSM:${1}\\tLB:$(zcat ${inputpath}${2}|head -n 1|cut -d":" -f3)\\tPL:Illumina" ${filepath}human_g1k_v37_decoy.fasta $inputpath$2 $inputpath$3 | samtools view -Sb - > $outputpath${1}.bam && \
samtools sort -@ 4 -m 20G -O bam -o $outputpath${1}.sorted.bam $outputpath${1}.bam && \
gatk MarkDuplicates -I $outputpath${1}.sorted.bam -O $outputpath${1}.sorted.markdup.bam -M $outputpath${1}.sorted.markdup.bam.metrics && \
samtools index $outputpath${1}.sorted.markdup.bam && \
gatk BaseRecalibrator -I $outputpath${1}.sorted.markdup.bam \
   -R $filepath/human_g1k_v37_decoy.fasta \
   -L $filepath/v6.bed \
   --known-sites  $filepath/1000G_phase1.snps.high_confidence.b37.vcf \
   --known-sites $filepath/Mills_and_1000G_gold_standard.indels.b37.vcf \
   --known-sites $filepath/dbsnp_138.b37.vcf \
   -O $outputpath${1}.table && \
gatk ApplyBQSR \
   -R $filepath/human_g1k_v37_decoy.fasta \
   -I $outputpath/${1}.sorted.markdup.bam \
   -L $filepath/v6.bed \
   --bqsr-recal-file $outputpath/${1}.table \
   -O $outputpath/${1}.sorted.markdup.BQSR.bam && \
rm $outputpath/${1}.bam $outputpath/${1}.sorted.bam $outputpath/${1}.sorted.markdup.bam* $outputpath/${1}.table
