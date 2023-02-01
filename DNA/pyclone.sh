cd /project/liuxd/liver/WES-analysis/code/pyclone/test
source ~/.bashrc 
conda activate pyclone
for i in DT04 DT06 DT07 DT09 DT10 DT12 DT13 DT14 DT17;do
python script/var_total_reads_from_vcf.py -v input/mpileup/${i}.vcf -m input/mutect/${i}/ -p ${i}


python script/seg_purity_file_from_sectors.py -s ${i}.sectors -p ${i} -d input/sequenza/

python script/assignCN_to_maf_multiSec.py -m ${i}_mutant_reads.tsv -s ${i}.seg  -p ${i} -d ${i}_read_depth.tsv 

sh script/pyclone_pipeline_1.sh ${i} ${i}.purity
done
