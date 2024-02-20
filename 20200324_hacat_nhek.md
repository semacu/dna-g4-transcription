
Yuqi is looking for genes with a strong bg4 peak in the promoter (hacat compared to nhek) and with gene expression in hacat higher than nhek too. We are going to look at Robert's Nature Genetics 2016 dataset again.


## Datasets

We are going to start by using the datasets that Robert did already generated in the paper available at [GSE76688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76688)

- Diff BG4 NHEK vs HaCaT: GSE76688_diff.BG4-hek_vs_hacat.txt.gz
- Diff expression NHEK vs HaCaT: GSE76688_expr_diff.hek_vs_hacat.rnaseq.txt.gz

Here they used hg19 reference genome



### Diff BG4

Download and uncompress:

```bash
cd /scratchb/sblab/martin03/repository/20200217_yuqi/data
mkdir -p 20200324/{bg4,rnaseq,promoters,pqs,oqs} && cd 20200324/bg4
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76688/suppl/GSE76688_diff.BG4-hek_vs_hacat.txt.gz
pigz -d GSE76688_diff.BG4-hek_vs_hacat.txt.gz
wc -l GSE76688_diff.BG4-hek_vs_hacat.txt # 21815
```

Explore:

```r
library(data.table)

# define window width
options(width = 300)

# load data
bg4 <- fread("/scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/bg4/GSE76688_diff.BG4-hek_vs_hacat.txt")

nrow(bg4[FDR < 0.05 & logFC < 0]) # 19934 bg4 in hacat and not in nhek
nrow(bg4[FDR < 0.05 & logFC > 0]) # 272 bg4 not in hacat and in nhek
## see Fig S9b

nrow(bg4[FDR < 0.05 & logFC < -2]) # 3080 bg4 in hacat (strong) and not in nhek
bg4[FDR < 0.05 & logFC < -2][order(logFC)]
#                          locus     logFC     logCPM       PValue          FDR chrom     start       end
#   1:   chr22_23162020_23162232 -9.911586  6.1576823 5.807101e-07 1.275948e-06 chr22  23162020  23162232
#   2: chr14_106641650_106641818 -8.857277  4.6741706 2.511215e-06 4.864545e-06 chr14 106641650 106641818
#   3: chr14_106330025_106330074 -8.829415  3.9057259 5.015758e-07 1.115442e-06 chr14 106330025 106330074
#   4:   chr22_23247170_23247208 -8.673293  3.3687984 1.326776e-06 2.713509e-06 chr22  23247170  23247208
#   5: chr14_107169957_107170042 -8.663476  4.6684813 1.724140e-06 3.447015e-06 chr14 107169957 107170042
#  ---                                                                                                   
#3076:    chr1_67518726_67518779 -2.000712 -3.1213023 5.572570e-05 8.509033e-05  chr1  67518726  67518779
#3077:   chr15_67439246_67439713 -2.000450  0.2875853 5.394612e-11 3.118952e-10 chr15  67439246  67439713
#3078:   chr11_65149256_65150414 -2.000171  1.7339887 1.798280e-11 1.183339e-10 chr11  65149256  65150414
#3079:   chr14_90863157_90863346 -2.000135 -1.4677445 1.373165e-12 1.226627e-11 chr14  90863157  90863346
#3080:   chr10_43892649_43892985 -2.000011 -1.2529965 2.013686e-11 1.307337e-10 chr10  43892649  43892985

bg4 <- bg4[, c("chrom", "start", "end", "logCPM", "logFC", "FDR")]

# write
write.table(bg4, file = "/scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/bg4/GSE76688_diff.BG4-hek_vs_hacat.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
```



### Diff expression

Download and uncompress:

```bash
cd /scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/rnaseq
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76688/suppl/GSE76688_expr_diff.hek_vs_hacat.rnaseq.txt.gz
pigz -d GSE76688_expr_diff.hek_vs_hacat.rnaseq.txt.gz
wc -l GSE76688_expr_diff.hek_vs_hacat.rnaseq.txt # 13044
```

Explore:

```r
library(data.table)

# define window width
options(width = 300)

# load data
rnaseq <- fread("/scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/rnaseq/GSE76688_expr_diff.hek_vs_hacat.rnaseq.txt")

nrow(rnaseq[FDR < 0.05 & logFC < 0]) # 5383 genes expressed in hacat and not in nhek
nrow(rnaseq[FDR < 0.05 & logFC > 0]) # 5660 genes expressed not in hacat and in nhek
## see Fig S9a

nrow(rnaseq[FDR < 0.05 & logFC < -2]) # 614 genes expressed in hacat (strong) and not in nhek
nrow(rnaseq[FDR < 0.05 & logFC > 2]) # 1442 genes expressed not in hacat and in nhek (strong)

rnaseq[FDR < 0.05 & logFC < -2][order(logFC)]
#          logFC     logCPM        PValue           FDR      gene_id
#  1: -11.449619  4.3276019 2.218355e-321 6.182243e-320         ABP1
#  2: -11.431622  4.3098337 7.977115e-311 2.132080e-309        GATA2
#  3: -11.288067  4.1668099  0.000000e+00  0.000000e+00      SLCO1B3
#  4: -10.950797  3.8349010 7.440980e-265 1.612171e-263      SCGB1A1
#  5: -10.824188  3.7102674 4.974184e-239 9.389041e-238        CALB1
# ---                                                               
#610:  -2.008609  0.4958863  2.180237e-08  3.567985e-08 LOC100188947
#611:  -2.007254  3.0196547  2.728807e-41  9.930756e-41       FAM49A
#612:  -2.006563 -0.1539439  1.570233e-05  2.274101e-05       NLRP12
#613:  -2.006050  5.1702161 1.997505e-167 2.493153e-166       ZNF367
#614:  -2.003248  5.6294539 1.256259e-183 1.781020e-182      TRMT61A
```



### Promoters

Using GenomicFeatures:

```r
library(data.table)
library(GenomicFeatures)

# change width
options(width = 250)

# prepare coordinates table
txdb <- makeTxDbFromGFF("/scratcha/sblab/martin03/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf", format="gtf")

# gene promoters
gene_promoters <- data.table(data.frame(promoters(genes(txdb, columns="gene_id"), upstream=1000, downstream=0)))[, c("seqnames", "start", "end", "gene_id", "strand")][order(seqnames, start)]
gene_promoters[, start := ifelse(start < 0, 0, start)]
gene_promoters[, dummy := "."]
gene_promoters <- gene_promoters[, .(seqnames, start, end, gene_id, dummy, strand)]

# write
write.table(gene_promoters, file = "/scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/promoters/hg19_genepromoters.bed", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
```



### Promoters + Diff expression

```r
library(data.table)

# change width
options(width = 250)

# load promoters
gene_promoters <- fread("/scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/promoters/hg19_genepromoters.bed")
nrow(gene_promoters) # 25721

# load diff expression
rnaseq <- fread("/scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/rnaseq/GSE76688_expr_diff.hek_vs_hacat.rnaseq.txt")
nrow(rnaseq) # 13043

# merge
promoters_rnaseq <- merge(gene_promoters, rnaseq, by.x = "gene_id", by.y = "gene_id")
promoters_rnaseq <- promoters_rnaseq[,c("seqnames", "start", "end", "gene_id", "dummy", "strand", "logCPM", "logFC", "FDR")]

# write
write.table(promoters_rnaseq, file = "/scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/promoters/hg19_genepromoters_rnaseq.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
```



### PQS

Map G3L7 and G3L12 in hg19:

```bash
srun --mem 64G --pty /usr/bin/bash

cd /scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/pqs

ref=/scratcha/sblab/martin03/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

nohup fastaRegexFinder.py -f $ref -q | \
bedtools sort -i | \
sed 's/$/&\tG3L7/' > hg19.g3l7.bed &

nohup fastaRegexFinder.py -f $ref -r '([gG]{3,}\w{1,12}){3,}[gG]{3,}' -q | \
bedtools sort -i | \
sed 's/$/&\tG3L12/' > hg19.g3l12.bed &

wc -l hg19*
   # 705580 hg19.g3l12.bed
   # 361424 hg19.g3l7.bed

exit
```



### OQS

In HEK293T, Marsico2019 [GSE110582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110582)

- GSM3003539_Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed.gz
- GSM3003539_Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed.gz

```bash
cd /scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/oqs

wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3003nnn/GSM3003539/suppl/GSM3003539_Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3003nnn/GSM3003539/suppl/GSM3003539_Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed.gz

pigz -d GSM3003539_Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed.gz
pigz -d GSM3003539_Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed.gz

wc -l GSM3003539_Homo_all_w15_th-1_*
  # 216321 GSM3003539_Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed
  # 217951 GSM3003539_Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed

bedtools merge -i <(cat GSM3003539_Homo_all_w15_th-1_* | bedtools sort -i) | \
sed 's/$/&\tOQS_HEK293T/' > hg19_oqs_hek293t.bed

wc -l hg19_oqs_hek293t.bed # 428624
```



### bam and bigwig files

```bash
cd /scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324

mkdir bam && cd bam

# HaCaT bg4-chip rep 1 and 2
rsync -arvuP martin03@10.20.192.21:/archive/Groups/SBLab/fs05/berald01/repository/bam_clean/rhh_25cyc_BG4_12082015.bam* .
rsync -arvuP martin03@10.20.192.21:/archive/Groups/SBLab/fs05/berald01/repository/bam_clean/rhh175_ChIPwthacat_704_502_entst_26082015.bam* .

# HaCaT input rep 1 and 2
rsync -arvuP martin03@10.20.192.21:/archive/Groups/SBLab/fs05/berald01/repository/bam_clean/rhh155_25cyc_input_703_503_12082015.bam* .
rsync -arvuP martin03@10.20.192.21:/archive/Groups/SBLab/fs05/berald01/repository/bam_clean/rhh178_inputwthacat_703_517_entst_26082015.bam* .

# NHEK bg4-chip rep 1 and 2
rsync -arvuP martin03@10.20.192.21:/archive/Groups/SBLab/fs05/berald01/repository/bam_clean/HEKnp_Lonza_1472015_BG4.md.bam* .
rsync -arvuP martin03@10.20.192.21:/archive/Groups/SBLab/fs05/berald01/repository/bam_clean/HEKnp_Lonza_1572015_BG4.md.bam* .

# NHEK input rep 1
rsync -arvuP martin03@10.20.192.21:/archive/Groups/SBLab/fs05/berald01/repository/bam_clean/merged_14_and_15072015_input_heknplonza.md.bam* .

mkdir ../bw

for bam in *.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../bw/$bname.log --mem 32G --wrap "bamCoverage -b $bam -o ../bw/$bname.bw -of bigwig --binSize 1 -p 20 --normalizeUsing CPM"
done

tail ../bw/*.log
```

Copy data across in C02Q70MUFVH8:

```bash
cd /Users/martin03/Desktop
mkdir bw && cd bw
rsync -arvuP martin03@10.20.236.34:/scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/bw/*.bw .
```

Through the ftp:

```bash
cd /scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/bw
rsync -arvuP *.bw martin03@10.20.192.94:/Users/martin03/tmp/data/

## In 10.20.192.94:
cd /Users/martin03/tmp/data/

# list
curl -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl ftp://ftp2.cruk.cam.ac.uk/

# upload
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -T HEKnp_Lonza_1472015_BG4.md.bw ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -T HEKnp_Lonza_1572015_BG4.md.bw ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -T merged_14_and_15072015_input_heknplonza.md.bw ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -T rhh155_25cyc_input_703_503_12082015.bw ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -T rhh175_ChIPwthacat_704_502_entst_26082015.bw ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -T rhh178_inputwthacat_703_517_entst_26082015.bw ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -T rhh_25cyc_BG4_12082015.bw ftp://ftp2.cruk.cam.ac.uk

# list
curl -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl ftp://ftp2.cruk.cam.ac.uk/

# download to C02Q70MUFVH8
curl -C - -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -O ftp://ftp2.cruk.cam.ac.uk/HEKnp_Lonza_1472015_BG4.md.bw
curl -C - -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -O ftp://ftp2.cruk.cam.ac.uk/HEKnp_Lonza_1572015_BG4.md.bw
curl -C - -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -O ftp://ftp2.cruk.cam.ac.uk/merged_14_and_15072015_input_heknplonza.md.bw
curl -C - -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -O ftp://ftp2.cruk.cam.ac.uk/rhh155_25cyc_input_703_503_12082015.bw
curl -C - -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -O ftp://ftp2.cruk.cam.ac.uk/rhh175_ChIPwthacat_704_502_entst_26082015.bw
curl -C - -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -O ftp://ftp2.cruk.cam.ac.uk/rhh178_inputwthacat_703_517_entst_26082015.bw
curl -C - -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -O ftp://ftp2.cruk.cam.ac.uk/rhh_25cyc_BG4_12082015.bw

# delete files from ftp
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -Q '-DELE HEKnp_Lonza_1472015_BG4.md.bw' ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -Q '-DELE HEKnp_Lonza_1572015_BG4.md.bw' ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -Q '-DELE merged_14_and_15072015_input_heknplonza.md.bw' ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -Q '-DELE rhh155_25cyc_input_703_503_12082015.bw' ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -Q '-DELE rhh175_ChIPwthacat_704_502_entst_26082015.bw' ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -Q '-DELE rhh178_inputwthacat_703_517_entst_26082015.bw' ftp://ftp2.cruk.cam.ac.uk
curl -v -k --user sblabcollabftp01:SimarXMTR7 --ftp-pasv --ftp-ssl -Q '-DELE rhh_25cyc_BG4_12082015.bw' ftp://ftp2.cruk.cam.ac.uk

# delete files in 10.20.192.94
rm *.bw
```





## Overlaps

```bash
cd /scratchb/sblab/martin03/repository/20200217_yuqi/data/20200324/

wc -l promoters/hg19_genepromoters_rnaseq.bed bg4/GSE76688_diff.BG4-hek_vs_hacat.bed pqs/hg19.g3l7.bed oqs/hg19_oqs_hek293t.bed
  #  12143 promoters/hg19_genepromoters_rnaseq.bed
  #  21814 bg4/GSE76688_diff.BG4-hek_vs_hacat.bed
  # 361424 pqs/hg19.g3l7.bed
  # 428624 oqs/hg19_oqs_hek293t.bed


###################
# Promoters + BG4 #
###################
bedtools intersect \
-a <(bedtools sort -i promoters/hg19_genepromoters_rnaseq.bed) \
-b <(bedtools sort -i bg4/GSE76688_diff.BG4-hek_vs_hacat.bed) \
-wa -wb > promoters/hg19_genepromoters_rnaseq_bg4.bed

wc -l promoters/hg19_genepromoters_rnaseq_bg4.bed # 6230
head promoters/hg19_genepromoters_rnaseq_bg4.bed | column -t
# chr1  714069   715068   LOC100288069  .  -  1.70073972775118  0.712095787872963  0.0039704562640092    chr1  713994   714275   -0.853554282666274  -0.709646059455032  0.00373339431541032
# chr1  714069   715068   LOC100288069  .  -  1.70073972775118  0.712095787872963  0.0039704562640092    chr1  714278   714318   -2.60909002482863   -0.969023073327991  0.0193739743896211
# chr1  935553   936552   HES4          .  -  3.09112982242932  -2.80443711793276  2.3473977964665e-68   chr1  935498   935645   -1.78707680203653   -1.53492543033284   1.39842696190092e-05
# chr1  954503   955502   AGRN          .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948   955029   -2.5032361225441    -0.968191769355195  0.0297658617891449
# chr1  1051737  1052736  C1orf159      .  -  3.67125566964693  0.344444905264951  0.00722011094849972   chr1  1051475  1051856  -0.328385414486771  -1.44758307874753   1.63824755823749e-06
# chr1  1166629  1167628  B3GALT6       .  +  4.30678680909463  -0.40141837613466  4.13514585714299e-06  chr1  1167367  1167633  -0.886141152070085  -1.46019700416577   1.99914289732987e-05
# chr1  1167448  1168447  SDF4          .  -  6.76736345215231  0.5740699422656    2.62729325589239e-39  chr1  1167367  1167633  -0.886141152070085  -1.46019700416577   1.99914289732987e-05
# chr1  1242994  1243993  PUSL1         .  +  4.60417093250567  0.373993489189774  4.00692632192369e-06  chr1  1243860  1243910  -3.04932914008297   -1.85395385456334   0.000804129625838727
# chr1  1243270  1244269  ACAP3         .  -  5.28654175499842  1.09017004823457   5.19272693025329e-37  chr1  1243860  1243910  -3.04932914008297   -1.85395385456334   0.000804129625838727
# chr1  1260068  1261067  CPSF3L        .  -  5.91161970474368  0.403550074415224  3.41849560995884e-15  chr1  1260021  1260224  -1.68316560531449   -1.38000960200687   2.56310576542915e-06

grep -P "MYC\t" promoters/hg19_genepromoters_rnaseq_bg4.bed | column -t
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128747406  128747778  -0.338059360217115  -1.80318345742238  1.40085883244716e-12
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128748029  128748502  -0.517575932025211  -2.20561688113167  1.99150639797418e-16

grep -P "KRAS\t" promoters/hg19_genepromoters_rnaseq_bg4.bed | column -t
# chr12  25403866  25404865  KRAS  .  -  6.20379552442811  0.334966162950006  1.40841443791145e-10  chr12  25403914  25404009  -2.1288132579731  -1.54795516419889  0.000164267503568793

# Columns:
# 1: chr promoter
# 2: start promoter
# 3: end promoter
# 4: gene name
# 5: dummy
# 6: strand promoter/gene
# 7: logCPM diff expression (NHEK/HaCaT)
# 8: logFC diff expression (NHEK/HaCaT)
# 9: FDR diff expression (NHEK/HaCaT)
# 10: chr bg4
# 11: start bg4
# 12: end bg4
# 13: logCPM diff bg4 (NHEK/HaCaT)
# 14: logFC diff bg4 (NHEK/HaCaT)
# 15: FDR diff bg4 (NHEK/HaCaT)



##########################
# Promoters + BG4 + G3L7 #
##########################
bedtools intersect \
-a promoters/hg19_genepromoters_rnaseq_bg4.bed \
-b pqs/hg19.g3l7.bed \
-wa -wb > promoters/hg19_genepromoters_rnaseq_bg4_pqsg3l7.bed

wc -l promoters/hg19_genepromoters_rnaseq_bg4_pqsg3l7.bed # 5655
head promoters/hg19_genepromoters_rnaseq_bg4_pqsg3l7.bed | column -t
# chr1  935553   936552   HES4      .  -  3.09112982242932  -2.80443711793276  2.3473977964665e-68   chr1  935498   935645   -1.78707680203653   -1.53492543033284   1.39842696190092e-05  chr1  935742   935790   chr1_935742_935790_for    48   +  GGGAGGGGCCGCTGGGGGCGAACGGGGCCCGGGACCCCCGGGGCTGGG                                                                                                        G3L7
# chr1  954503   955502   AGRN      .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948   955029   -2.5032361225441    -0.968191769355195  0.0297658617891449    chr1  954852   954886   chr1_954852_954886_for    34   +  ggggcctgggggggcggggccgggaggggcgggg                                                                                                                      G3L7
# chr1  954503   955502   AGRN      .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948   955029   -2.5032361225441    -0.968191769355195  0.0297658617891449    chr1  955018   955040   chr1_955018_955040_for    22   +  GGGATCGGGGCCGGGTCTGGGG                                                                                                                                  G3L7
# chr1  954503   955502   AGRN      .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948   955029   -2.5032361225441    -0.968191769355195  0.0297658617891449    chr1  955155   955188   chr1_955155_955188_for    33   +  ggggggccgcggcgggggaggggcgccTGCGGG                                                                                                                       G3L7
# chr1  954503   955502   AGRN      .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948   955029   -2.5032361225441    -0.968191769355195  0.0297658617891449    chr1  955215   955250   chr1_955215_955250_rev    35   -  CCCCGCCCCCTCCCGGGCCGCCCCTTTCCCGCCCC                                                                                                                     G3L7
# chr1  954503   955502   AGRN      .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948   955029   -2.5032361225441    -0.968191769355195  0.0297658617891449    chr1  955312   955353   chr1_955312_955353_for    41   +  gggcgggccgggggagggaggagggagggggcgggaggggg                                                                                                               G3L7
# chr1  954503   955502   AGRN      .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948   955029   -2.5032361225441    -0.968191769355195  0.0297658617891449    chr1  955365   955406   chr1_955365_955406_for    41   +  gggggaggaggagggcgcggggagggggtaggggggcgggg                                                                                                               G3L7
# chr1  954503   955502   AGRN      .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948   955029   -2.5032361225441    -0.968191769355195  0.0297658617891449    chr1  955429   955461   chr1_955429_955461_rev    32   -  ccctggccccgcccccccccgcccgcccgccc                                                                                                                        G3L7
# chr1  1051737  1052736  C1orf159  .  -  3.67125566964693  0.344444905264951  0.00722011094849972   chr1  1051475  1051856  -0.328385414486771  -1.44758307874753   1.63824755823749e-06  chr1  1051763  1051799  chr1_1051763_1051799_for  36   +  GGGAGACGGGGCGGGGCGTCGTCGGGGTCTCCCGGG                                                                                                                    G3L7
# chr1  1051737  1052736  C1orf159  .  -  3.67125566964693  0.344444905264951  0.00722011094849972   chr1  1051475  1051856  -0.328385414486771  -1.44758307874753   1.63824755823749e-06  chr1  1052466  1052616  chr1_1052466_1052616_rev  150  -  cccttccaacccctcagaccccccaaaacccccagatcccccaaccccccatatacccccaacccctcagaccccccaaccccccagaccccccaacctcccaggcccccagaccccccagCTCCCTTGCTCCCGGGAGCTCCCAGGCCC  G3L7

grep -P "MYC\t" promoters/hg19_genepromoters_rnaseq_bg4_pqsg3l7.bed | column -t
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128747406  128747778  -0.338059360217115  -1.80318345742238  1.40085883244716e-12  chr8  128747402  128747429  chr8_128747402_128747429_for  27  +  GGGAGGCGTGGGGGTGGGACGGTGGGG                             G3L7
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128747406  128747778  -0.338059360217115  -1.80318345742238  1.40085883244716e-12  chr8  128748168  128748222  chr8_128748168_128748222_rev  54  -  ccccaccttccccaccctccccaccctccccaTAAGCGCCCCTCCCGGGTTCCC  G3L7
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128748029  128748502  -0.517575932025211  -2.20561688113167  1.99150639797418e-16  chr8  128747402  128747429  chr8_128747402_128747429_for  27  +  GGGAGGCGTGGGGGTGGGACGGTGGGG                             G3L7
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128748029  128748502  -0.517575932025211  -2.20561688113167  1.99150639797418e-16  chr8  128748168  128748222  chr8_128748168_128748222_rev  54  -  ccccaccttccccaccctccccaccctccccaTAAGCGCCCCTCCCGGGTTCCC  G3L7

grep -P "KRAS\t" promoters/hg19_genepromoters_rnaseq_bg4_pqsg3l7.bed | column -t
# chr12  25403866  25404865  KRAS  .  -  6.20379552442811  0.334966162950006  1.40841443791145e-10  chr12  25403914  25404009  -2.1288132579731  -1.54795516419889  0.000164267503568793  chr12  25403982  25404010  chr12_25403982_25404010_for  28  +  gggcggtgtgggaagagggaagaggggg                          G3L7
# chr12  25403866  25404865  KRAS  .  -  6.20379552442811  0.334966162950006  1.40841443791145e-10  chr12  25403914  25404009  -2.1288132579731  -1.54795516419889  0.000164267503568793  chr12  25404027  25404079  chr12_25404027_25404079_for  52  +  ggggagaaggagggggccgggccgggccggcgggggaggagcgggggccggg  G3L7

# Columns:
# 1: chr promoter
# 2: start promoter
# 3: end promoter
# 4: gene name
# 5: dummy
# 6: strand promoter/gene
# 7: logCPM diff expression (NHEK/HaCaT)
# 8: logFC diff expression (NHEK/HaCaT)
# 9: FDR diff expression (NHEK/HaCaT)
# 10: chr bg4
# 11: start bg4
# 12: end bg4
# 13: logCPM diff bg4 (NHEK/HaCaT)
# 14: logFC diff bg4 (NHEK/HaCaT)
# 15: FDR diff bg4 (NHEK/HaCaT)
# 16: chr pqs
# 17: start pqs
# 18: end pqs
# 19: id pqs
# 20: length pqs
# 21: strand pqs
# 22: sequence pqs
# 23: type pqs



################################
# Promoters + BG4 + G3L7 + OQS #
################################
bedtools intersect \
-a promoters/hg19_genepromoters_rnaseq_bg4_pqsg3l7.bed \
-b oqs/hg19_oqs_hek293t.bed \
-wa -wb > promoters/hg19_genepromoters_rnaseq_bg4_pqsg3l7_oqs.bed

wc -l promoters/hg19_genepromoters_rnaseq_bg4_pqsg3l7_oqs.bed # 7544
head promoters/hg19_genepromoters_rnaseq_bg4_pqsg3l7_oqs.bed | column -t
# chr1  954503  955502  AGRN  .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948  955029  -2.5032361225441  -0.968191769355195  0.0297658617891449  chr1  954852  954886  chr1_954852_954886_for  34  +  ggggcctgggggggcggggccgggaggggcgggg         G3L7  chr1  954838  954958  OQS_HEK293T
# chr1  954503  955502  AGRN  .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948  955029  -2.5032361225441  -0.968191769355195  0.0297658617891449  chr1  954852  954886  chr1_954852_954886_for  34  +  ggggcctgggggggcggggccgggaggggcgggg         G3L7  chr1  955318  955471  OQS_HEK293T
# chr1  954503  955502  AGRN  .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948  955029  -2.5032361225441  -0.968191769355195  0.0297658617891449  chr1  955018  955040  chr1_955018_955040_for  22  +  GGGATCGGGGCCGGGTCTGGGG                     G3L7  chr1  954838  954958  OQS_HEK293T
# chr1  954503  955502  AGRN  .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948  955029  -2.5032361225441  -0.968191769355195  0.0297658617891449  chr1  955018  955040  chr1_955018_955040_for  22  +  GGGATCGGGGCCGGGTCTGGGG                     G3L7  chr1  955318  955471  OQS_HEK293T
# chr1  954503  955502  AGRN  .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948  955029  -2.5032361225441  -0.968191769355195  0.0297658617891449  chr1  955155  955188  chr1_955155_955188_for  33  +  ggggggccgcggcgggggaggggcgccTGCGGG          G3L7  chr1  954838  954958  OQS_HEK293T
# chr1  954503  955502  AGRN  .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948  955029  -2.5032361225441  -0.968191769355195  0.0297658617891449  chr1  955155  955188  chr1_955155_955188_for  33  +  ggggggccgcggcgggggaggggcgccTGCGGG          G3L7  chr1  955318  955471  OQS_HEK293T
# chr1  954503  955502  AGRN  .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948  955029  -2.5032361225441  -0.968191769355195  0.0297658617891449  chr1  955215  955250  chr1_955215_955250_rev  35  -  CCCCGCCCCCTCCCGGGCCGCCCCTTTCCCGCCCC        G3L7  chr1  954838  954958  OQS_HEK293T
# chr1  954503  955502  AGRN  .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948  955029  -2.5032361225441  -0.968191769355195  0.0297658617891449  chr1  955215  955250  chr1_955215_955250_rev  35  -  CCCCGCCCCCTCCCGGGCCGCCCCTTTCCCGCCCC        G3L7  chr1  955318  955471  OQS_HEK293T
# chr1  954503  955502  AGRN  .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948  955029  -2.5032361225441  -0.968191769355195  0.0297658617891449  chr1  955312  955353  chr1_955312_955353_for  41  +  gggcgggccgggggagggaggagggagggggcgggaggggg  G3L7  chr1  954838  954958  OQS_HEK293T
# chr1  954503  955502  AGRN  .  +  9.49301838920165  0.946806876691389  1.61488189675275e-66  chr1  954948  955029  -2.5032361225441  -0.968191769355195  0.0297658617891449  chr1  955312  955353  chr1_955312_955353_for  41  +  gggcgggccgggggagggaggagggagggggcgggaggggg  G3L7  chr1  955318  955471  OQS_HEK293T

grep -P "MYC\t" promoters/hg19_genepromoters_rnaseq_bg4_pqsg3l7_oqs.bed | column -t
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128747406  128747778  -0.338059360217115  -1.80318345742238  1.40085883244716e-12  chr8  128747402  128747429  chr8_128747402_128747429_for  27  +  GGGAGGCGTGGGGGTGGGACGGTGGGG                             G3L7  chr8  128747399  128747504  OQS_HEK293T
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128747406  128747778  -0.338059360217115  -1.80318345742238  1.40085883244716e-12  chr8  128747402  128747429  chr8_128747402_128747429_for  27  +  GGGAGGCGTGGGGGTGGGACGGTGGGG                             G3L7  chr8  128748111  128748216  OQS_HEK293T
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128747406  128747778  -0.338059360217115  -1.80318345742238  1.40085883244716e-12  chr8  128748168  128748222  chr8_128748168_128748222_rev  54  -  ccccaccttccccaccctccccaccctccccaTAAGCGCCCCTCCCGGGTTCCC  G3L7  chr8  128747399  128747504  OQS_HEK293T
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128747406  128747778  -0.338059360217115  -1.80318345742238  1.40085883244716e-12  chr8  128748168  128748222  chr8_128748168_128748222_rev  54  -  ccccaccttccccaccctccccaccctccccaTAAGCGCCCCTCCCGGGTTCCC  G3L7  chr8  128748111  128748216  OQS_HEK293T
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128748029  128748502  -0.517575932025211  -2.20561688113167  1.99150639797418e-16  chr8  128747402  128747429  chr8_128747402_128747429_for  27  +  GGGAGGCGTGGGGGTGGGACGGTGGGG                             G3L7  chr8  128747399  128747504  OQS_HEK293T
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128748029  128748502  -0.517575932025211  -2.20561688113167  1.99150639797418e-16  chr8  128747402  128747429  chr8_128747402_128747429_for  27  +  GGGAGGCGTGGGGGTGGGACGGTGGGG                             G3L7  chr8  128748111  128748216  OQS_HEK293T
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128748029  128748502  -0.517575932025211  -2.20561688113167  1.99150639797418e-16  chr8  128748168  128748222  chr8_128748168_128748222_rev  54  -  ccccaccttccccaccctccccaccctccccaTAAGCGCCCCTCCCGGGTTCCC  G3L7  chr8  128747399  128747504  OQS_HEK293T
# chr8  128747315  128748314  MYC  .  +  9.03573538887895  -0.325504191957396  3.56894745921627e-38  chr8  128748029  128748502  -0.517575932025211  -2.20561688113167  1.99150639797418e-16  chr8  128748168  128748222  chr8_128748168_128748222_rev  54  -  ccccaccttccccaccctccccaccctccccaTAAGCGCCCCTCCCGGGTTCCC  G3L7  chr8  128748111  128748216  OQS_HEK293T

grep -P "KRAS\t" promoters/hg19_genepromoters_rnaseq_bg4_pqsg3l7_oqs.bed | column -t
# chr12  25403866  25404865  KRAS  .  -  6.20379552442811  0.334966162950006  1.40841443791145e-10  chr12  25403914  25404009  -2.1288132579731  -1.54795516419889  0.000164267503568793  chr12  25403982  25404010  chr12_25403982_25404010_for  28  +  gggcggtgtgggaagagggaagaggggg                          G3L7  chr12  25404054  25404189  OQS_HEK293T
# chr12  25403866  25404865  KRAS  .  -  6.20379552442811  0.334966162950006  1.40841443791145e-10  chr12  25403914  25404009  -2.1288132579731  -1.54795516419889  0.000164267503568793  chr12  25404027  25404079  chr12_25404027_25404079_for  52  +  ggggagaaggagggggccgggccgggccggcgggggaggagcgggggccggg  G3L7  chr12  25404054  25404189  OQS_HEK293T

# Columns:
# 1: chr promoter
# 2: start promoter
# 3: end promoter
# 4: gene name
# 5: dummy
# 6: strand promoter/gene
# 7: logCPM diff expression (NHEK/HaCaT)
# 8: logFC diff expression (NHEK/HaCaT)
# 9: FDR diff expression (NHEK/HaCaT)
# 10: chr bg4
# 11: start bg4
# 12: end bg4
# 13: logCPM diff bg4 (NHEK/HaCaT)
# 14: logFC diff bg4 (NHEK/HaCaT)
# 15: FDR diff bg4 (NHEK/HaCaT)
# 16: chr pqs
# 17: start pqs
# 18: end pqs
# 19: id pqs
# 20: length pqs
# 21: strand pqs
# 22: sequence pqs
# 23: type pqs
# 24: chr oqs
# 25: start oqs
# 26: end oqs
# 27: type oqs
```





## Clean-up

```bash
cd /scratchb/sblab/martin03/repository/20200217_yuqi/data
rm -rf 20200324
```
