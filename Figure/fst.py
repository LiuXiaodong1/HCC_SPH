import pandas as pd
import os
from itertools import combinations

os.chdir("E:/cancer genome/liver/分析/Fig/1/distance")

# alist代表variant allele的list，nlist是total allele的list。
# 一次只计算一个位点的fst，比如有两个个体，total allele分别是20,30，variant allele是10,12，那么fst_WC84([10,12],[20,30])；
# 群体遗传学应该是anc dri
# 群体的fst算平均值
def fst_WC84(alist, nlist):  # Bruce weir's book, 173, formula T1/T2,
    #    print alist,nlist
    if 0 in nlist or sum(nlist) == 2:
        return 0.0
    plist = []  # freqlist
    numpop = len(alist)  # number of population
    r = numpop * 1.0
    pbar = sum(alist) * 1.0 / sum(nlist)
    if pbar == 0 or pbar == 1:
        return 0
    MSP = 0.0
    i = 0
    while i < numpop:
        f = alist[i] * 1.0 / nlist[i]
        plist = plist + [f]  # creating plist
        MSP = MSP + nlist[i] * (f - pbar) * (f - pbar)
        i = i + 1
    MSP = MSP / (r - 1.0)

    nbar = sum(nlist) / r
    sumN2 = 0.0
    for x in nlist:
        sumN2 = sumN2 + x * x

    nc = (r * nbar - sumN2 / (r * nbar)) / (r - 1)  # calculate Nc
    SAsquare = MSP / nbar  # this is the formula on 171/172
    T1 = SAsquare - 1.0 / (nbar - 1) * (pbar * (1 - pbar) - (r - 1) / r * SAsquare)
    T2 = (nc - 1) / (nbar - 1) * pbar * (1 - pbar) + (1 + (r - 1) * (nbar - nc) / (nbar - 1)) * SAsquare / r
    #    print T1,T2
    Fst = T1 / T2
    return Fst


for patient in ["DT03","DT04","DT06","DT07","DT09","DT10","DT12","DT13","DT14","DT16","DT17","DT18","DT19"]:

    vcf = pd.read_table(patient+".vcf", skiprows=range(0, 107))
    variants = pd.read_table(patient+"_variants_pos.txt")
    variants_depth = pd.read_table(patient+"_variants_data.txt")
    tmp = pd.notna(variants_depth).all(1)
    filter1 = [x for x, i in enumerate(tmp) if i]
    filter2 = variants[
        variants["Reference_Allele"].isin(["A", "T", "C", "G"]) & variants["Tumor_Seq_Allele2"].isin(["A", "T", "C", "G"])]
    select = filter2.loc[[x for x in filter2.index if x not in filter1], :]

    tmp=vcf.loc[:, "#CHROM"].str.cat(vcf.loc[:, "POS"].astype(str), sep="_").isin(
        select.loc[:, "Chromosome"].str.cat(select.loc[:, "Start_position"].astype(str),sep="_"))

    vcf_select = vcf[tmp]
    vcf_select = vcf_select.sort_values(["#CHROM","POS"])
    select = select.sort_values(["Chromosome","Start_position"])

    sample_list=vcf_select.columns[vcf_select.columns.str.match("DT.*T")]


    combines=list(combinations(sample_list,2))

    sample1=[];sample2=[];fst=[]
    for combine in combines:
        tmp1=vcf_select[list(combine)].iloc[:,0].str.split(":")
        tmp2=vcf_select[list(combine)].iloc[:,1].str.split(":")
        ref1=[int(a[1].split(",")[0]) for a in tmp1]
        alt1=[int(a[1].split(",")[1]) for a in tmp1]
        all1=[ref1[a]+alt1[a] for a in range(len(ref1))]

        ref2=[int(a[1].split(",")[0]) for a in tmp2]
        alt2=[int(a[1].split(",")[1]) for a in tmp2]
        all2=[ref2[a]+alt2[a] for a in range(len(ref2))]

        fst_sum=0
        for i in range(len(all1)):
            all=[all1[i],all2[i]]
            alt=[ref1[i],ref2[i]]
            fst_single=fst_WC84(alt, all)
            fst_sum+=fst_single
        fst.append(fst_sum/len(all1))
        sample1.append(combine[0])
        sample2.append(combine[1])

    fst_res=pd.DataFrame({"sample1":sample1,"sample2":sample2,"fst":fst})
    fst_res.to_csv(patient+"_fst.tsv",sep="\t",index=False)
