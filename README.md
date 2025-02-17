# RNAseq_Medicago

 Differential Expression analysis of two Medicago mutants

## 0. Data Preparation

### 0.1 Reference

Medicago truncatula v4.0 genome and associated annotation (from https://phytozome-next.jgi.doe.gov/info/Mtruncatula_Mt4_0v1) were used as reference. The reference genome was indexed using `novoindex`.

```bash
cd /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/
curl --cookie jgi_session=/api/sessions/ad1b3f6d89afedaa910f96b575be09c5 --output download.20241120.095555.zip -d "{\"ids\":{\"Phytozome-285\":[\"53112a7f49607a1be0055936\",\"53112a7e49607a1be0055934\"]}}" -H "Content-Type: application/json" https://files-download.jgi.doe.gov/filedownload/
/rds/user/cx264/hpc-work/program/novocraft/novoindex Mtruncatula_285_Mt4.0.index Mtruncatula_285_Mt4.0.fa
# get gtf file from gff3
/rds/user/cx264/rds-csc_programmes-FTKWLWDeHys/programs/gffread/gffread-0.12.7.Linux_x86_64/gffread Mtruncatula_285_Mt4.0v1.gene.gff3 -T -o Mtruncatula_285_Mt4.0v1.gene.gtf
```

### 0.2 Directories

```bash
for i in {1..3}; do mkdir Jester_mock_12hrs_rep${i} Jester_spotinoc_12hrs_rep${i} nin-1_mock_12hrs_rep${i} nin-1_spotinoc_12hrs_rep${i} nin-1_mock_24hrs_rep${i} nin-1_spotinoc_24hrs_rep${i} nsp1-1_mock_12hrs_rep${i} nsp1-1_spotinoc_12hrs_rep${i} nsp1-1_mock_24hrs_rep${i} nsp1-1_spotinoc_24hrs_rep${i} nsp2-2_mock_12hrs_rep${i} nsp2-2_spotinoc_12hrs_rep${i} nsp2-2_mock_24hrs_rep${i} nsp2-2_spotinoc_24hrs_rep${i}; done
```

## 1. Reads Quality Control and Trimming

Raw RNA-seq reads were subjected to `Trimmomatic` (v0.39) to filter out reads of low quality. Specifically, bases from the start and end of a read with quality bellow 20 were trimmed, and reads with length shorter than 40 bp were discarded, by setting `LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40`.

```bash
for file in /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/0.data/*; do
  if [ -d "$file" ]; then
    sample_name=$(basename "$file")
    fq="/rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/rawdata/${sample_name}.fastq.gz"
    cd ${sample_name}

      # Run Trimmomatic
    echo "Processing ${sample_name}"
    java -jar /usr/local/software/trimmomatic/0.39/trimmomatic-0.39.jar SE -threads 4 -summary "${sample_name}.summary" "$fq" ${sample_name}.fq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40
    echo "${sample_name} finished"
    cd /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/0.data/
  fi
done

echo "Trimmomatic process completed."
```

## 2.Reads Mapping

`Novoalign` V4.03.08 with default parameters was used to align filtered reads to reference genome. The alignment in bam format was sorted using `novosort`. `qualimap rnaseq` function in `qualimap` (v2.3) was used to check quality of the alignment. `htseq-count` (2.0.5) was used to calculate expression counts.

```bash
cd /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/1.mapping/

export PATH=/rds/user/cx264/rds-csc_programmes-FTKWLWDeHys/programs/java/jdk-17.0.10/bin:/rds/project/rds-FTKWLWDeHys/programs/mambaforge/bin:$PATH

# Run novoalign iterate through each sample
for file in /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/0.data/*; do
  if [ -d "$file" ]; then
    sample_name=$(basename "$file")
    fq="/rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/0.data/${sample_name}/${sample_name}.fq.gz"
    cd ${sample_name}

    # Run novoalign
    echo "Processing ${sample_name}"


    # Change to the working directory
    cd /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/1.mapping/"${sample_name}"
    echo "Processing ${sample_name}"
    echo "${fq}"

    # Map with novoalign
    /rds/user/cx264/hpc-work/program/novocraft/novoalign -d /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0.index -f "${fq}" -o BAM "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:illumina" -k > 1.raw.bam

    # Create a tmp directory for novosort
    mkdir -p tmp

    # Sort and process with novosort
    /rds/user/cx264/hpc-work/program/novocraft/novosort -m 32G --threads 16 --tmpdir ./tmp --output 2.sorted.bam --index --bai --markDuplicates --stats 2.duplicate.summary 1.raw.bam

    # Run mapinsights bamqc
    /rds/user/cx264/hpc-work/program/mapinsights/mapinsights bamqc -r /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0.fa -i 2.sorted.bam -o ./
    # Run qualimap
    /rds/user/cx264/rds-csc_programmes-FTKWLWDeHys/programs/qualimap_v2.3/qualimap rnaseq -bam 2.sorted.bam -gtf /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.gtf -outdir ./qualimap/ -outformat pdf


    # Expression Counts
    # HTseq-count
    /rds/user/cx264/hpc-work/program/miniforge3/envs/python3.8/bin/htseq-count --type transcript --counts_output 3.sorted.bam.count.tsv --nprocesses 16 --max-reads-in-buffer 1000000 2.sorted.bam /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.gtf

    echo "${sample_name} finished"
    cd /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/1.mapping/

  fi
done

echo "Alignment process completed."
```

## 3.Clustering

```bash
cd /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/2.clustering

/rds/project/rds-FTKWLWDeHys/programs/mambaforge/bin/multiBamSummary BED-file --BED /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.bed --bamfiles ../1.mapping/Jester_mock_12hrs_rep1/2.sorted.bam ../1.mapping/nsp1-1_mock_24hrs_rep1/2.sorted.bam ../1.mapping/Jester_mock_12hrs_rep2/2.sorted.bam ../1.mapping/nsp1-1_mock_24hrs_rep2/2.sorted.bam ../1.mapping/Jester_mock_12hrs_rep3/2.sorted.bam ../1.mapping/nsp1-1_mock_24hrs_rep3/2.sorted.bam ../1.mapping/Jester_spotinoc_12hrs_rep1/2.sorted.bam ../1.mapping/nsp1-1_spotinoc_12hrs_rep1/2.sorted.bam ../1.mapping/Jester_spotinoc_12hrs_rep2/2.sorted.bam ../1.mapping/nsp1-1_spotinoc_12hrs_rep2/2.sorted.bam ../1.mapping/Jester_spotinoc_12hrs_rep3/2.sorted.bam ../1.mapping/nsp1-1_spotinoc_12hrs_rep3/2.sorted.bam ../1.mapping/nin-1_mock_12hrs_rep1/2.sorted.bam ../1.mapping/nsp1-1_spotinoc_24hrs_rep1/2.sorted.bam ../1.mapping/nin-1_mock_12hrs_rep2/2.sorted.bam ../1.mapping/nsp1-1_spotinoc_24hrs_rep2/2.sorted.bam ../1.mapping/nin-1_mock_12hrs_rep3/2.sorted.bam ../1.mapping/nsp1-1_spotinoc_24hrs_rep3/2.sorted.bam ../1.mapping/nin-1_mock_24hrs_rep1/2.sorted.bam ../1.mapping/nsp2-2_mock_12hrs_rep1/2.sorted.bam ../1.mapping/nin-1_mock_24hrs_rep2/2.sorted.bam ../1.mapping/nsp2-2_mock_12hrs_rep2/2.sorted.bam ../1.mapping/nin-1_mock_24hrs_rep3/2.sorted.bam ../1.mapping/nsp2-2_mock_12hrs_rep3/2.sorted.bam ../1.mapping/nin-1_spotinoc_12hrs_rep1/2.sorted.bam ../1.mapping/nsp2-2_mock_24hrs_rep1/2.sorted.bam ../1.mapping/nin-1_spotinoc_12hrs_rep2/2.sorted.bam ../1.mapping/nsp2-2_mock_24hrs_rep2/2.sorted.bam ../1.mapping/nin-1_spotinoc_12hrs_rep3/2.sorted.bam ../1.mapping/nsp2-2_mock_24hrs_rep3/2.sorted.bam ../1.mapping/nin-1_spotinoc_24hrs_rep1/2.sorted.bam ../1.mapping/nsp2-2_spotinoc_12hrs_rep1/2.sorted.bam ../1.mapping/nin-1_spotinoc_24hrs_rep2/2.sorted.bam ../1.mapping/nsp2-2_spotinoc_12hrs_rep2/2.sorted.bam ../1.mapping/nin-1_spotinoc_24hrs_rep3/2.sorted.bam ../1.mapping/nsp2-2_spotinoc_12hrs_rep3/2.sorted.bam ../1.mapping/nsp1-1_mock_12hrs_rep1/2.sorted.bam ../1.mapping/nsp2-2_spotinoc_24hrs_rep1/2.sorted.bam ../1.mapping/nsp1-1_mock_12hrs_rep2/2.sorted.bam ../1.mapping/nsp2-2_spotinoc_24hrs_rep2/2.sorted.bam ../1.mapping/nsp1-1_mock_12hrs_rep3/2.sorted.bam ../1.mapping/nsp2-2_spotinoc_24hrs_rep3/2.sorted.bam --outFileName 1.correlation.matrix --metagene -p max --outRawCounts 2.rawCounts.tab --verbose --labels Jester_mock_12hrs_rep1 nsp1-1_mock_24hrs_rep1 Jester_mock_12hrs_rep2 nsp1-1_mock_24hrs_rep2 Jester_mock_12hrs_rep3 nsp1-1_mock_24hrs_rep3 Jester_spotinoc_12hrs_rep1 nsp1-1_spotinoc_12hrs_rep1 Jester_spotinoc_12hrs_rep2 nsp1-1_spotinoc_12hrs_rep2 Jester_spotinoc_12hrs_rep3 nsp1-1_spotinoc_12hrs_rep3 nin-1_mock_12hrs_rep1 nsp1-1_spotinoc_24hrs_rep1 nin-1_mock_12hrs_rep2 nsp1-1_spotinoc_24hrs_rep2 nin-1_mock_12hrs_rep3 nsp1-1_spotinoc_24hrs_rep3 nin-1_mock_24hrs_rep1 nsp2-2_mock_12hrs_rep1 nin-1_mock_24hrs_rep2 nsp2-2_mock_12hrs_rep2 nin-1_mock_24hrs_rep3 nsp2-2_mock_12hrs_rep3 nin-1_spotinoc_12hrs_rep1 nsp2-2_mock_24hrs_rep1 nin-1_spotinoc_12hrs_rep2 nsp2-2_mock_24hrs_rep2 nin-1_spotinoc_12hrs_rep3 nsp2-2_mock_24hrs_rep3 nin-1_spotinoc_24hrs_rep1 nsp2-2_spotinoc_12hrs_rep1 nin-1_spotinoc_24hrs_rep2 nsp2-2_spotinoc_12hrs_rep2 nin-1_spotinoc_24hrs_rep3 nsp2-2_spotinoc_12hrs_rep3 nsp1-1_mock_12hrs_rep1 nsp2-2_spotinoc_24hrs_rep1 nsp1-1_mock_12hrs_rep2 nsp2-2_spotinoc_24hrs_rep2 nsp1-1_mock_12hrs_rep3 nsp2-2_spotinoc_24hrs_rep3 --numberOfProcessors 16

/rds/project/rds-FTKWLWDeHys/programs/mambaforge/bin/plotCorrelation --corData 1.correlation.matrix --corMethod spearman --whatToPlot heatmap --plotFile 3.spearman.heatmap1.pdf --plotTitle spearman --labels Jester_mock_12hrs_rep1 nsp1-1_mock_24hrs_rep1 Jester_mock_12hrs_rep2 nsp1-1_mock_24hrs_rep2 Jester_mock_12hrs_rep3 nsp1-1_mock_24hrs_rep3 Jester_spotinoc_12hrs_rep1 nsp1-1_spotinoc_12hrs_rep1 Jester_spotinoc_12hrs_rep2 nsp1-1_spotinoc_12hrs_rep2 Jester_spotinoc_12hrs_rep3 nsp1-1_spotinoc_12hrs_rep3 nin-1_mock_12hrs_rep1 nsp1-1_spotinoc_24hrs_rep1 nin-1_mock_12hrs_rep2 nsp1-1_spotinoc_24hrs_rep2 nin-1_mock_12hrs_rep3 nsp1-1_spotinoc_24hrs_rep3 nin-1_mock_24hrs_rep1 nsp2-2_mock_12hrs_rep1 nin-1_mock_24hrs_rep2 nsp2-2_mock_12hrs_rep2 nin-1_mock_24hrs_rep3 nsp2-2_mock_12hrs_rep3 nin-1_spotinoc_12hrs_rep1 nsp2-2_mock_24hrs_rep1 nin-1_spotinoc_12hrs_rep2 nsp2-2_mock_24hrs_rep2 nin-1_spotinoc_12hrs_rep3 nsp2-2_mock_24hrs_rep3 nin-1_spotinoc_24hrs_rep1 nsp2-2_spotinoc_12hrs_rep1 nin-1_spotinoc_24hrs_rep2 nsp2-2_spotinoc_12hrs_rep2 nin-1_spotinoc_24hrs_rep3 nsp2-2_spotinoc_12hrs_rep3 nsp1-1_mock_12hrs_rep1 nsp2-2_spotinoc_24hrs_rep1 nsp1-1_mock_12hrs_rep2 nsp2-2_spotinoc_24hrs_rep2 nsp1-1_mock_12hrs_rep3 nsp2-2_spotinoc_24hrs_rep3 --colorMap hot_r --plotNumbers --skipZeros --removeOutliers --zMin 0.6
```

## 4.Differential Expression Analysis

For DE analysis, my code is mainly from edgeR user's guide (https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf). This user's guide is a must-read in order to know the rationale behind the code.
We start from concatenating all counts:

```bash
paste -d"\t" ../1.mapping/Jester_mock_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/Jester_mock_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/Jester_mock_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/Jester_spotinoc_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/Jester_spotinoc_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/Jester_spotinoc_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nin-1_mock_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nin-1_mock_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nin-1_mock_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nin-1_mock_24hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nin-1_mock_24hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nin-1_mock_24hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_24hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_24hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_24hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_mock_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_mock_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_mock_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_mock_24hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_mock_24hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_mock_24hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_24hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_24hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_24hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_mock_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_mock_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_mock_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_mock_24hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_mock_24hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_mock_24hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_spotinoc_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_spotinoc_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_spotinoc_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_spotinoc_24hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_spotinoc_24hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_spotinoc_24hrs_rep3/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12"\t"$14"\t"$16"\t"$18"\t"$20"\t"$22"\t"$24"\t"$26"\t"$28"\t"$30"\t"$32"\t"$34"\t"$36"\t"$38"\t"$40"\t"$42"\t"$44"\t"$46"\t"$48"\t"$50"\t"$52"\t"$54"\t"$56"\t"$58"\t"$60"\t"$62"\t"$64"\t"$66"\t"$68"\t"$70"\t"$72"\t"$74"\t"$76"\t"$78"\t"$80"\t"$82"\t"$84}' > 01.raw.counts.tsv
```

### Jester_spotinoc_12hrs VS Jester_mock_12hrs

```R
cd /data/pathology/cxia/projects/Giles/Jinpeng/3.DEanalysis/01.Jester_spotinoc_12hrs_VS_Jester_mock_12hrs
paste -d"\t" ../1.mapping/Jester_mock_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/Jester_mock_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/Jester_mock_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/Jester_spotinoc_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/Jester_spotinoc_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/Jester_spotinoc_12hrs_rep3/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12}' > 01.raw.counts.tsv

mamba activate R4.2.3
R
library(EnhancedVolcano)
library(gplots)
library(RColorBrewer)
library("ggplot2")
library(edgeR)
library(stringr)
library(statmod)
library(readr)

my_count_matrix <- read_tsv("01.raw.counts.tsv")
# we have two groups (Jester_mock_12hrs, Jester_spotinoc_12hrs), and each group has three samples
# so we design experiment:
group <- c(1,1,1,2,2,2)
sample <- factor(c("1","2","3","1","2","3"))
y <- DGEList(counts=my_count_matrix[,c(2:7)], group=group,genes=my_count_matrix[,1])

# then we filter genes that are not expressed in either experimental condition
# you may check and change the options for your purpose
# since these data are sparse, we apply less strict criteria for filtering
keep <- filterByExpr(y, min.count = 0.01, min.total.count = 0.01, large.n = 2, min.prop = 0.01)
table(keep)
#FALSE  TRUE
#47486  3413
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)
# TMM normalization is applied to account for the compositional biases
y<-normLibSizes(y)
y$samples
#plot clustering of treatment and control
pdf("1.1.LRT.MDS.pdf")
plotMDS(y, col=rep(1:2, each=3))
dev.off()
```

Following is an MDS plot shows the relative similarities of the six samples. This gives you another way to evaluate your experiment quality.

((MDS plot))

Next, we need to define our design matrix based on the experimental design.

```R
design <- model.matrix(~sample+group)
rownames(design) <- colnames(y)
design
#                           (Intercept) sample2 sample3 group
#Jester_mock_12hrs_rep1               1       0       0     1
#Jester_mock_12hrs_rep2               1       1       0     1
#Jester_mock_12hrs_rep3               1       0       1     1
#Jester_spotinoc_12hrs_rep1           1       0       0     2
#Jester_spotinoc_12hrs_rep2           1       1       0     2
#Jester_spotinoc_12hrs_rep3           1       0       1     2
#attr(,"assign")
#[1] 0 1 1 2
#attr(,"contrasts")
#attr(,"contrasts")$sample
#[1] "contr.treatment"
```

Then we do a `Dispersion estimation` step. Dispersion estimation is performed to account for variability in gene expression data that is not explained by the mean expression levels alone. This variability can be due to biological differences, technical noise, or other factors.

```R
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
#[1] 0.151497
pdf("1.2.LRT.dispersion.BCV.pdf")
plotBCV(y)
dev.off()
```

![](C:\Users\cx264\AppData\Roaming\marktext\images\2025-01-21-16-18-30-image.png)

```R
fit <- glmQLFit(y, design, robust=TRUE)
pdf("1.3.LRT.QLdispersion.pdf")
plotQLDisp(fit)
dev.off()
```

((QTLdispersion figure))

Then we do test for significant differential expression in each gene using the QLF-test.

```R
qlf <- glmQLFTest(fit)
topTags(qlf)
cpm(y)[rownames(topTags(qlf)),]
summary(decideTests(qlf,lfc=0.5))
#       group
#Down       4
#NotSig  3408
#Up         1

# following is to save results with change information, FDR
allTags2 <- topTags(qlf, n = nrow(qlf$genes), adjust.method = "BH", sort.by = "none", p.value = 1)
DEG_glmQLF <- as.data.frame(allTags2)
# differential expression was defined as FDR < 0.05 and log2FC>0.5 or
log2FC<-0.5
k1 <- (DEG_glmQLF$FDR < 0.05) & (DEG_glmQLF$logFC < -0.5)
k2 <- (DEG_glmQLF$FDR < 0.05) & (DEG_glmQLF$logFC > 0.5)
DEG_glmQLF$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "Not_significant"))
table(DEG_glmQLF$change)
DEG_glmQLF$CPM <- cpm(y)[rownames(qlf),]
write.csv(DEG_glmQLF, "1.4.QLF.DE.results.Jester_12hrs_spotinoc_vs_mock.FINAL.csv", quote = F)


## Visualization
deseq_results <- read.csv("1.4.QLF.DE.results.Jester_12hrs_spotinoc_vs_mock.FINAL.csv",sep=",", header = T)
diffGenes <- deseq_results[deseq_results$FDR < 0.05,][,c(2,9,10,11,12,13,14)]
clustRows <- hclust(as.dist(1-cor(t(diffGenes[,c(2,3,4,5,6,7)]),method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(diffGenes[,c(2,3,4,5,6,7)],method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9)
module.color <- module.color[as.vector(module.assign)]
#myheatcolors1 <- blueorange(75)
myheatcolors1 <- colorRampPalette(c("blue", "orange"))(25)
pdf("1.5.QLF.DiffGene.heatmap.pdf")
heatmap.2(data.matrix(diffGenes[,c(2,3,4,5,6,7)]),
    Rowv=as.dendrogram(clustRows),
    Colv=as.dendrogram(clustColumns),
    RowSideColors=module.color,
    col=myheatcolors1, scale='row', labRow=NA,
    density.info="none", trace="none",
    cexRow=1, cexCol=1, margins=c(8,20), keysize=1, srtCol=45)
dev.off()
png("1.5.QLF.DiffGene.heatmap.png",width = 636, height = 980)
heatmap.2(data.matrix(diffGenes[,c(2,3,4,5,6,7)]),
    Rowv=as.dendrogram(clustRows),
    Colv=as.dendrogram(clustColumns),
    RowSideColors=module.color,
    col=myheatcolors1, scale='row', labRow=NA,
    density.info="none", trace="none",
    cexRow=1, cexCol=1, margins=c(8,20), keysize=0.5, srtCol=45)
dev.off()

pdf("1.6.QLF.DiffGene.volcano.pdf", 7, 7)
EnhancedVolcano(deseq_results, lab = rownames(deseq_results), x = 'logFC', y = 'FDR', pCutoff = 5e-2,
    FCcutoff=1.0, ylim = c(0, 3), xlim = c(-8, 6), pointSize = 1.0, labSize = 0,
    colAlpha = 1, legendIconSize=2.0, legendLabSize = 10) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
png("1.6.QLF.DiffGene.volcano.png", 800, 800)
EnhancedVolcano(deseq_results, lab = rownames(deseq_results), x = 'logFC', y = 'FDR', pCutoff = 5e-2,
    FCcutoff=1.0, ylim = c(0, 3), xlim = c(-8, 6), pointSize = 1.0, labSize = 0,
    colAlpha = 1, legendIconSize=2.0, legendLabSize = 10) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
```

![](C:\Users\cx264\AppData\Roaming\marktext\images\2025-01-21-16-17-21-image.png)

### nin-1_spotinoc_12hrs VS nin-1_mock_12hrs

```bash
paste -d"\t" ../1.mapping/nin-1_mock_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nin-1_mock_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nin-1_mock_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_12hrs_rep3/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12}' > 02.nin-1_spotinoc_12hrs_VS_mock/01.raw.counts.tsv
```

### nin-1_spotinoc_24hrs_VS_mock

```bash
paste -d"\t" ../1.mapping/nin-1_mock_24hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nin-1_mock_24hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nin-1_mock_24hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_24hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_24hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nin-1_spotinoc_24hrs_rep3/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12}' > 03.nin-1_spotinoc_24hrs_VS_mock/01.raw.counts.tsv
```

nsp1-1_spotinoc_12hrs_VS_mock

```bash
paste -d"\t" ../1.mapping/nsp1-1_mock_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_mock_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_mock_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_12hrs_rep3/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12}' > 04.nsp1-1_spotinoc_12hrs_VS_mock/01.raw.counts.tsv
```

nsp1-1_spotinoc_24hrs_VS_mock

```bash
paste -d"\t" ../1.mapping/nsp1-1_mock_24hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_mock_24hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_mock_24hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_24hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_24hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp1-1_spotinoc_24hrs_rep3/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12}' > 05.nsp1-1_spotinoc_24hrs_VS_mock/01.raw.counts.tsv
```

```bash
paste -d"\t" ../1.mapping/nsp2-2_mock_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_mock_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_mock_12hrs_rep3/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_spotinoc_12hrs_rep1/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_spotinoc_12hrs_rep2/3.sorted.bam.count.tsv ../1.mapping/nsp2-2_spotinoc_12hrs_rep3/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12}' > 06.nsp2-2_spotinoc_12hrs_VS_mock/01.raw.counts.tsv
```

# X204SC24110826-Z01-F001

This is another batch, 36 samples.

## 0.Data Preparation

```bash
cd /rds/project/rds-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/02.X204SC24110826-Z01-F001
wget https://objectstorage.uk-london-1.oraclecloud.com/p/v2XIprSU1KrkEYuelgeTVUQkjA7HVOx09Y-r7FtRp-nBlElrsUS7k8rGeRCPkBIy/n/cnyr09qj8zbo/b/england-data/o/out/CP2021051100119/X204SC24110826-Z01-F001/X204SC24110826-Z01-F001.tar
wget https://objectstorage.uk-london-1.oraclecloud.com/p/v2XIprSU1KrkEYuelgeTVUQkjA7HVOx09Y-r7FtRp-nBlElrsUS7k8rGeRCPkBIy/n/cnyr09qj8zbo/b/england-data/o/out/CP2021051100119/X204SC24110826-Z01-F001/MD5.txt
tar -xvf X204SC24110826-Z01-F001.tar

#reference genome annotation gtf to bed12 for ReQC analyses
/rds-d7/user/cx264/rds-csc_programmes-FTKWLWDeHys/programs/gffread/gffread-0.12.7.Linux_x86_64/gffread --bed -o /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.bed12 /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.gff3
```

## 1.Trimming
Be aware that two or more libraries were sequenced for some samples.

```bash
#!/bin/bash

# we are using trimmomatic.jar from /usr/local/software/trimmomatic/0.39/trimmomatic-0.39.jar
input_dir="/rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/02.X204SC24110826-Z01-F001/X204SC24110826-Z01-F001/01.RawData"
output_dir="/rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/02.X204SC24110826-Z01-F001/X204SC24110826-Z01-F001/02.CleanData"

# Run Trimmomatic iterate through each sample
for fq1 in "$input_dir"/*/*_1.fq.gz; do
    fq2="${fq1/_1.fq.gz/_2.fq.gz}" # Get the corresponding paired-end file
    cd "${output_dir}"

    # Extract sample name
    sample_name=$(basename "$fq1" | cut -d'_' -f1-3)

    # Set directory
    if [ -d "$sample_name" ]; then
        # If the directory exists, change to that directory
        cd "$sample_name"
        echo "Changed directory to $sample_name"
            # Set output files
        output_fwd_paired="${sample_name}_02_1P.fq.gz"
        output_fwd_unpaired="${sample_name}_02_1U.fq.gz"
        output_rev_paired="${sample_name}_02_2P.fq.gz"
        output_rev_unpaired="${sample_name}_02_2U.fq.gz"
        echo "Processing ${sample_name}_02"
        java -jar /usr/local/software/trimmomatic/0.39/trimmomatic-0.39.jar PE -threads 8 -summary "${sample_name}_02.summary" "$fq1" "$fq2" \
        "$output_fwd_paired" "$output_fwd_unpaired" "$output_rev_paired" "$output_rev_unpaired" \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:60
        echo "${sample_name}_02 finished"
        cd "${output_dir}"
    else
        # If the directory does not exist, create it
        mkdir -p "${sample_name}"
        echo "Directory $sample_name created"
        # Change to the newly created directory
        cd "${sample_name}"
        echo "Changed directory to $sample_name"
        output_fwd_paired="${sample_name}_1P.fq.gz"
        output_fwd_unpaired="${sample_name}_1U.fq.gz"
        output_rev_paired="${sample_name}_2P.fq.gz"
        output_rev_unpaired="${sample_name}_2U.fq.gz"
        # Run Trimmomatic
        echo "Processing ${sample_name}"
        java -jar /usr/local/software/trimmomatic/0.39/trimmomatic-0.39.jar PE -threads 8 -summary "${sample_name}.summary" "$fq1" "$fq2" \
        "$output_fwd_paired" "$output_fwd_unpaired" "$output_rev_paired" "$output_rev_unpaired" \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:60
        cd "${output_dir}"
    fi
    echo "${sample_name} finished"
done

echo "Trimmomatic process completed."
```

## 2.Mapping

```bash
#!/bin/bash
#SBATCH --partition=icelake
#!SBATCH --partition=cclake,icelake,sapphire,icelake-himem,cclake-himem
#SBATCH -A EVESVANDENAKKER-SL2-CPU
#SBATCH --job-name=novoalign
#SBATCH --output=00.o
#SBATCH --error=00.e
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cx264@cam.ac.uk

cd /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/02.X204SC24110826-Z01-F001/X204SC24110826-Z01-F001/03.Mapping/

export JAVA_HOME=/rds/user/cx264/rds-csc_programmes-FTKWLWDeHys/programs/java/jdk-11.0.22+7/
export PATH=$JAVA_HOME/bin:/rds/project/rds-FTKWLWDeHys/programs/mambaforge/bin:/rds/project/rds-FTKWLWDeHys/programs/mambaforge/bin:$PATH

# Step 0: Set working and input directories
working_dir="/rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/02.X204SC24110826-Z01-F001/X204SC24110826-Z01-F001/03.Mapping/"
input_dir="/rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/02.X204SC24110826-Z01-F001/X204SC24110826-Z01-F001/02.CleanData/"

# Step 1: Pre-processing reads, since some samples have multiple libraries sequenced
# Iterate over samples
for sample_dir in "$input_dir"/*/; do
    # Extract the sample name from the directory path
    sample_name=$(basename "$sample_dir")
    cd ${input_dir}/${sample_name}

    # Find all the forward and reverse read files
    fq1_files=("${input_dir}"/"${sample_name}"/*_1P.fq.gz)
    fq2_files=("${input_dir}"/"${sample_name}"/*_2P.fq.gz)

    # Check if there is more than one file to concatenate
    if [ "${#fq1_files[@]}" -gt 1 ] && [ "${#fq2_files[@]}" -gt 1 ]; then
        # Concatenate the forward reads
        zcat "${fq1_files[@]}" | gzip > "${sample_dir}/${sample_name}_cat_1P.fq.gz"
        # Concatenate the reverse reads
        zcat "${fq2_files[@]}" | gzip > "${sample_dir}/${sample_name}_cat_2P.fq.gz"

        fq1=("${input_dir}"/"${sample_name}"/*_cat_1P.fq.gz)
        fq2=("${input_dir}"/"${sample_name}"/*_cat_2P.fq.gz)
    else
        # If only one pair, use the files directly
        fq1="${fq1_files[0]}"
        fq2="${fq2_files[0]}"
    fi

    # Step 2: Mapping reads to reference genome with novoalign
    echo "Processing Sample ${sample_name}"

    # Change to the working directory
    cd "${working_dir}"

    # Create a directory for the sample and move into it
    mkdir -p "${sample_name}" && cd "${sample_name}"
    echo "Processing ${sample_name}"
    echo "${fq1}" "${fq2}"

    /rds/user/cx264/hpc-work/program/novocraft/novoalign -d /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0.index -f "${fq1}" "${fq2}" -o BAM "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:illumina" -i PE 250,50 -k > 1.raw.bam

    # Create a tmp directory for novosort
    mkdir -p tmp

    # Sort and process with novosort
    /rds/user/cx264/hpc-work/program/novocraft/novosort -m 32G --threads 16 --tmpdir ./tmp --output 2.sorted.bam --index --bai --markDuplicates --stats 2.duplicate.summary 1.raw.bam

    # Step 3: Quality-control of alignments
    # Run mapinsights bamqc
    /rds/user/cx264/hpc-work/program/mapinsights/mapinsights bamqc -r /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0.fa -i 2.sorted.bam -o ./
    # Run qualimap
    /rds/user/cx264/rds-csc_programmes-FTKWLWDeHys/programs/qualimap_v2.3/qualimap rnaseq -bam 2.sorted.bam -gtf /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.gtf -outdir ./qualimap/ -outformat pdf

    # Step 4: Expression counts
    # HTseq-count
    /rds/user/cx264/rds-csc_programmes-FTKWLWDeHys/programs/mambaforge/envs/python3.8/bin/htseq-count --type transcript --counts_output 3.sorted.bam.count.tsv --nprocesses 16 --max-reads-in-buffer 1000000 2.sorted.bam /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.gtf

    echo "${sample_name} completed"
    cd "${working_dir}"
done

Concatenate all samples and only keep genes that are expressed in at least two samples.
```python
import pandas as pd
import os
import sys

def combine_expression_data(base_directory, output_file):
    """
    Combines expression data from multiple tab-delimited files located in sample directories into a single TSV file.

    Args:
        base_directory (str): The base directory containing sample subdirectories.
        output_file (str): The path to save the combined output TSV file.
    """
    dataframes = []

    # Traverse through each sample directory
    for sample_dir in os.listdir(base_directory):
        sample_path = os.path.join(base_directory, sample_dir)

        # Check if it's a directory
        if os.path.isdir(sample_path):
            # Look for .tsv files in the current sample directory
            for filename in os.listdir(sample_path):
                if filename.endswith('.tsv'):
                    file_path = os.path.join(sample_path, filename)

                    # Read the TSV file without headers
                    df = pd.read_csv(file_path, sep='\t', header=None)

                    # Extract SampleID from filename (assuming it's in the format SampleID.sorted.bam.count.tsv)
                    sample_id = file_path.split('/')[-1]  # Adjust this if your naming convention differs

                    # Rename columns: first column as GeneID and second column as SampleID
                    df.columns = ['GeneID', sample_id]

                    df.set_index('GeneID', inplace=True)
                    dataframes.append(df)

    # Concatenate all DataFrames along columns (axis=1)
    combined_df = pd.concat(dataframes, axis=1)
    combined_df.reset_index(inplace=True)

    # Filter out genes with no expression
    combined_df = combined_df[(combined_df.iloc[:, 1:] != 0).any(axis=1)]
    # Filter out genes expressed in only one sample
    non_zero_counts = (combined_df.iloc[:, 1:] > 0).sum(axis=1)
    combined_df = combined_df[non_zero_counts > 1]

    # Save the combined DataFrame to a new TSV file
    combined_df.to_csv(output_file, sep='\t', index=False)

    print(f"Combined expression data saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <base_directory> <output_file>")
        sys.exit(1)

    base_directory = sys.argv[1]
    output_file = sys.argv[2]

    combine_expression_data(base_directory, output_file)
```

```bash
cd /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/02.X204SC24110826-Z01-F001/X204SC24110826-Z01-F001/03.Mapping
python3 concate.expression.counts.py /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/02.X204SC24110826-Z01-F001/X204SC24110826-Z01-F001/03.Mapping raw.counts2.tsv
```


## 3.QC

In this pipeline, I do quality control throughout many steps. For example, after `trimmomatic` step, you can check how many reads (forward and/or reverse) survived the trimming criteria. Also in mapping step, I used `qualimap` and `mapinsights` to check quality/summary of alignment results, statistics including 'total reads aligned','Strand ratio','Mean mapping quality', etc.

Here I use two more approaches to evaluate quality of the experiment, `RSeQC` and Clustering. Details on `RSeQC` is [RSeQC: An RNA-seq Quality Control Package &#8212; RSeQC documentation](https://rseqc.sourceforge.net/#). Some specific statistics I found valuable are sequencing saturation, mapped reads distribution, reads quality, nucleotide composition bias, PCR bias and GC bias. 

### 3.1ReQC

```bash
#!/bin/bash

cd /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/02.X204SC24110826-Z01-F001/X204SC24110826-Z01-F001/03.Mapping
working_dir="/rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/02.X204SC24110826-Z01-F001/X204SC24110826-Z01-F001/03.Mapping"


for sample_dir in "$working_dir"/*/; do
    # Set-up directory path
    sample_name=$(basename "$sample_dir")
    cd ${working_dir}/${sample_name}
    mkdir RSeQC

    # Step1: Calculate the distributions of clipped nucleotides across reads
    clipping_profile.py -i 2.sorted.bam -s "PE" -o RSeQC/01.${sample_name}

    # Step2: Calculate the distributions of deletions across reads
    # set -l, length of original reads
    python3 /rds/project/rds-FTKWLWDeHys/programs/mambaforge/bin/deletion_profile.py -i 2.sorted.bam -l 150 -o RSeQC/02.${sample_name}

    # Step3: Calculate the RNA-seq reads coverage over gene body.
    # note gene model in 12-column bed format
    geneBody_coverage.py -r /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.bed12 -i 2.sorted.bam -o RSeQC/03.${sample_name}
    # Step4: reads were stranded for strand-specific RNA-seq data
    infer_experiment.py -r /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.bed12 -i 2.sorted.bam > RSeQC/04.${sample_name}.experiment.txt
    # Step5: Calculate inner distance between read pairs.
    inner_distance.py -i 2.sorted.bam -o RSeQC/05.${sample_name} -r /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.bed12

    # Step6: Calculate the distributions of inserted nucleotides across reads.
    python3 /rds/project/rds-FTKWLWDeHys/programs/mambaforge/bin/insertion_profile.py -s "PE" -i 2.sorted.bam -o RSeQC/06.${sample_name}  
    # Step7: splice junctions
    junction_annotation.py -i 2.sorted.bam -o RSeQC/07.${sample_name} -r /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.bed12 
    # Step8: check if current sequencing depth is deep enough to perform alternative splicing analyses
    junction_saturation.py -i 2.sorted.bam -o RSeQC/08.${sample_name} -r /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.bed12
    # Step9: Calculate the distribution of mismatches across reads.
    python3 /rds/project/rds-FTKWLWDeHys/programs/mambaforge/bin/mismatch_profile.py -i 2.sorted.bam -l 150 -o RSeQC/09.${sample_name}
    # Step10: calculate how mapped reads were distributed over genome feature (like CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions)
    read_distribution.py -i 2.sorted.bam -r /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.bed12 > RSeQC/10.${sample_name}.read_distribution.txt
    # Step11: reads duplication rate
    read_duplication.py -i 2.sorted.bam -o RSeQC/11.${sample_name}
    # Step12: GC content distribution of reads
    read_GC.py -i 2.sorted.bam -o RSeQC/12.${sample_name}
    # Step13: check the nucleotide composition bias
    read_NVC.py -i 2.sorted.bam -o RSeQC/13.${sample_name}
    # Step14: read quality
    read_quality.py -i 2.sorted.bam -o RSeQC/14.${sample_name}
    # Step15: check if the current sequencing depth was saturated or not (or if the RPKM values were stable or not) in terms of genes’ expression estimation
    # note default non-strand-specific RNA-seq    
    RPKM_saturation.py -r /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.bed12 -i 2.sorted.bam -o RSeQC/15.${sample_name} 

    echo "${sample_name} completed"
    cd "${working_dir}"
done
```

### 3.2Clustering

## 4.Differential Expression analysis
### 4.1 