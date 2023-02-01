from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
matrices = matGen.SigProfilerMatrixGeneratorFunc("HCC", "GRCh37",".",plot=True,exome=False,chrom_based=False, tsb_stat=False,seqInfo=False, cushion=100,bed_file="v6.bed")
from SigProfilerExtractor import sigpro as sig
data = "output/SBS/HCC.SBS96.region"
sig.sigProfilerExtractor("table", "HCCoutputSBS",data,reference_genome="GRCh37",cpu=40,nmf_replicates=100,opportunity_genome = "GRCh37",minimum_signatures=1, maximum_signatures=5, resample = True,cosmic_version=3.2,get_all_signature_matrices= True)
data = "output/DBS/HCC.DBS78.region"
sig.sigProfilerExtractor("table", "HCCoutputDBS",data,reference_genome="GRCh37",cpu=40,nmf_replicates=100,opportunity_genome = "GRCh37",minimum_signatures=1, maximum_signatures=5, resample = True,cosmic_version=3.2,get_all_signature_matrices= True)
data = "output/ID/HCC.ID83.region"
sig.sigProfilerExtractor("table", "HCCoutputID",data,reference_genome="GRCh37",cpu=40,nmf_replicates=100,opportunity_genome = "GRCh37",minimum_signatures=1, maximum_signatures=5, resample = True,cosmic_version=3.2,get_all_signature_matrices= True)
