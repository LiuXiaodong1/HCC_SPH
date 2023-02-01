#This is the shell commands for each modules 
base_dir=$(pwd)
mpileup_dir=../liver-case/mpileup_results
mutect_dir=../liver-case/mutect_results
seqz_dir=../liver-case/sequenza_results
src_dir=../scripts
patient=pat_test
pyclone_out=output
citup_out=citup_out
mapscape_out=mapscape_out
citup_type="QIP"


#This script extract the union of all mutations from all sectors
python ${src_dir}/var_total_reads_from_vcf.py -v ${mpileup_dir}/WHT314.wes.vcf -p ${patient} -m ${mutect_dir}/

#extract segament info for each sectors 
python ${src_dir}/seg_purity_file_from_sectors.py -s ${patient}.sectors -p ${patient} -d ${seqz_dir}/

#This combines mutation and CNV info and formulate pyclone input files
python ${src_dir}/assignCN_to_maf_multiSec.py -m ${patient}_mutant_reads.tsv -s ${patient}.seg -p ${patient} -d ${patient}_read_depth.tsv

#Run pyclone
bash ${src_dir}/pyclone_pipeline_1.sh ${patient} ${patient}.purity 

#Run citup
mkdir ${citup_out}
mkdir ${mapscape_out}
cd ${citup_out}
Rscript ${base_dir}/${src_dir}/calculate_ITH_clonal_subclonal_1.R ${base_dir}/${pyclone_out}/${patient}
if [ $citup_type == "QIP" ];then
	bash ${base_dir}/${src_dir}/run_citup_qip.sh ${patient} 
	python ${base_dir}/${src_dir}/get_optimal_tree.py ${patient}_results.h5 ${patient} QIP
	#cd ${base_dir}/${mapscape_out}
	Rscript ${base_dir}/${src_dir}/plot_tree_1.R ${patient} qip ${base_dir}/${src_dir}
else
	bash ${base_dir}/${src_dir}/run_citup.sh ${patient} 
	python ${base_dir}/${src_dir}/get_optimal_tree.py ${patient}_results.h5 ${patient} iter
	Rscript ${base_dir}/${src_dir}/plot_tree_1.R ${patient} iter ${base_dir}/${src_dir}
fi

mv *.html ${base_dir}/${mapscape_out}
cd ${base_dir}
