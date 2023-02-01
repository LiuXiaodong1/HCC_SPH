module load gatk/4.1.3.0
#use:${1}:Tumor name ${2}:normal.bam ${3}:tumor.bam 
inputpath="/project/liuxd/liver/WES-analysis/output/bwa/"
filepath="/project/liuxd/liver/WES-analysis/input/"
outputpath="/project/liuxd/liver/WES-analysis/output/mutect/"

gatk Mutect2 -R ${filepath}human_g1k_v37_decoy.fasta \
   -L ${filepath}v6.bed \
   -I ${inputpath}${2} \
   -I ${inputpath}${3} \
   -normal $(echo $2|sed "s/\\..*//") \
   -germline-resource ${filepath}af-only-gnomad.raw.sites.b37.vcf.gz \
   -pon ${outputpath}/DTpon.vcf.gz \
   --f1r2-tar-gz ${outputpath}${1}f1r2.tar.gz \
   -O ${outputpath}${1}unfiltered.vcf \
&& \
gatk LearnReadOrientationModel -I ${outputpath}${1}f1r2.tar.gz -O ${outputpath}${1}read-orientation-model.tar.gz \
&& \
gatk GetPileupSummaries \
    -I ${inputpath}${3} \
    -R ${filepath}human_g1k_v37_decoy.fasta \
    -V ${filepath}small_exac_common_3_b37.vcf.gz \
    -L ${filepath}v6.bed \
    -O ${outputpath}${1}getpileupsummaries.table \
&& \
gatk CalculateContamination \
    -I ${outputpath}${1}getpileupsummaries.table \
    -tumor-segmentation ${outputpath}${1}segments.table \
    -O ${outputpath}${1}calculatecontamination.table \
&& \
gatk FilterMutectCalls -V ${outputpath}${1}unfiltered.vcf \
   -R ${filepath}human_g1k_v37_decoy.fasta \
   --tumor-segmentation ${outputpath}${1}segments.table \
   --contamination-table ${outputpath}${1}calculatecontamination.table \
   --ob-priors ${outputpath}${1}read-orientation-model.tar.gz \
   -O ${outputpath}${1}filtered.vcf
