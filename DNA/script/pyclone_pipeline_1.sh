#export PATH="/home/chuakp/anaconda3/bin:$PATH"
#
#export QT_QPA_PLATFORM='offscreen'
#
prior=major_copy_number
iters=20000
burnin=2000

# source activate $pycEnv
#source activate pyclone

# Sample name
sam=$1
sam_input_dir=${sam}_pyclone_input

# Purity file, make sure second column is purity
pur=$2

echo $(pwd)

mkdir -p output/${sam}

rm -f output/${sam}/purity.txt

#for sec in $(ls ${sam}*.tsv)
for sec in $(ls ${sam_input_dir}/*.tsv)
do
	# Remove prefix
	filename=$(basename ${sec})
	# Remove extension
	secName=${filename%.*}
	echo ${secName}
	purity=$(grep -P "${secName}\t" ${pur})
	purity=($purity)
	echo ${purity[1]} >> output/${sam}/purity.txt
done

PyClone run_analysis_pipeline --density pyclone_binomial --burnin ${burnin} --prior ${prior} --num_iters ${iters} --tumour_contents $(cat output/${sam}/purity.txt) --in_files $(ls ${sam_input_dir}/*.tsv) --working_dir output/${sam} --min_cluster_size 2 &> output/${sam}/setup_analysis.log 
