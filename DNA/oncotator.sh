. ~/.bashrc
conda activate oncotator
cd /project/liuxd/liver/WES-analysis/output/statistic
for i in $(ls);do
	oncotator -v --db-dir /project/liuxd/npc/hg19/input/oncotator_v1_ds_April052016 -i VCF ${i}.snp.recode.vcf ${i}.snp.maf hg19
	oncotator -v --db-dir /project/liuxd/npc/hg19/input/oncotator_v1_ds_April052016 -i VCF ${i}.indel.recode.vcf ${i}.indel.maf hg19
done
