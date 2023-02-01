source ~/.bashrc
conda activate pyclone
cd /project/liuxd/liver/WES-analysis/code/pyclone/test
#This is the shell commands for each modules 
base_dir=$(pwd)
src_dir=script
patient=pat_test
pyclone_out=output
mapscape_out=mapscape_out
citup_type="QIP"

patient="DT03"
citup_out=citup_out${patient}
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