source ~/.bashrc
conda activate sequenza
module load samtools/1.9
inputpath="/project/liuxd/liver/WES-analysis/output/bwa/"
filepath="/project/liuxd/liver/WES-analysis/input/"
outputpath="/project/liuxd/liver/WES-analysis/output/sequenza/"

sequenza-utils bam2seqz -n $inputpath$2 -t $inputpath$3 --fasta $filepath/human_g1k_v37_decoy.fasta -gc $filepath/genome_gc50.wig.gz -o $outputpath/${1}out.seqz.gz

sequenza-utils seqz_binning --seqz $outputpath/${1}out.seqz.gz -w 50 -o $outputpath/${1}small.seqz.gz

Rscript /project/liuxd/liver/WES-analysis/code/sequenza/sequenza.R $1
