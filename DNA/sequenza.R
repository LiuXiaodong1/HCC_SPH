library(sequenza)
setwd("/project/liuxd/liver/WES-analysis/output/sequenza/")
args<-commandArgs(T)
seqname<-paste0(args[1],"small.seqz.gz")
data.file <-paste0("/project/liuxd/liver/WES-analysis/output/sequenza/",seqname)
extract <- sequenza.extract(data.file, verbose = FALSE)
CP <- sequenza.fit(extract)
sequenza.results(sequenza.extract = extract,
    cp.table = CP, sample.id = args[1],
    out.dir=args[1])
