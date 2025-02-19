# RNAseq_Medicago

Differential Expression analysis of Medicago mutants
Project iD: X204SC24110826-Z01-F001. We have 36 samples, 10 conditions, some of them have only one or two replicates.

## Table of Contents
- [0. Data Preparation](https://github.com/chongjing/RNAseq_Medicago#0-data-preparation)
- [1. Trimming](https://github.com/chongjing/RNAseq_Medicago#1-trimming)
- [2. Mapping](https://github.com/chongjing/RNAseq_Medicago#2-mapping)
- [3. QC](https://github.com/chongjing/RNAseq_Medicago#3-qc)
- [4. Differential Expression Analysis](https://github.com/chongjing/RNAseq_Medicago#4-differential-expression-analysis)
- [5. Gene Expression Clustering](https://github.com/chongjing/RNAseq_Medicago/tree/main#5-gene-expression-clustering)
## 0. Data Preparation

### 0.1 Reference

Medicago truncatula v4.0 genome and associated annotation (from https://phytozome-next.jgi.doe.gov/info/Mtruncatula_Mt4_0v1) were used as reference. The reference genome was indexed using `bwa-mem2, samtools, novoindex` for subsequent analyses.

```bash
cd /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/
curl --cookie jgi_session=/api/sessions/ad1b3f6d89afedaa910f96b575be09c5 --output download.20241120.095555.zip -d "{\"ids\":{\"Phytozome-285\":[\"53112a7f49607a1be0055936\",\"53112a7e49607a1be0055934\"]}}" -H "Content-Type: application/json" https://files-download.jgi.doe.gov/filedownload/
/rds/user/cx264/hpc-work/program/novocraft/novoindex Mtruncatula_285_Mt4.0.index Mtruncatula_285_Mt4.0.fa
bwa-mem2 index Mtruncatula_285_Mt4.0.fa
samtools faidx Mtruncatula_285_Mt4.0.fa

# get gtf file from gff3
/rds/user/cx264/rds-csc_programmes-FTKWLWDeHys/programs/gffread/gffread-0.12.7.Linux_x86_64/gffread Mtruncatula_285_Mt4.0v1.gene.gff3 -T -o Mtruncatula_285_Mt4.0v1.gene.gtf
# reference genome annotation gtf to bed12 for ReQC analyses
/rds-d7/user/cx264/rds-csc_programmes-FTKWLWDeHys/programs/gffread/gffread-0.12.7.Linux_x86_64/gffread --bed -o /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.bed12 /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.gff3
```

### 0.Raw Data Download

```bash
cd /rds/project/rds-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/02.X204SC24110826-Z01-F001
wget https://objectstorage.uk-london-1.oraclecloud.com/p/v2XIprSU1KrkEYuelgeTVUQkjA7HVOx09Y-r7FtRp-nBlElrsUS7k8rGeRCPkBIy/n/cnyr09qj8zbo/b/england-data/o/out/CP2021051100119/X204SC24110826-Z01-F001/X204SC24110826-Z01-F001.tar
wget https://objectstorage.uk-london-1.oraclecloud.com/p/v2XIprSU1KrkEYuelgeTVUQkjA7HVOx09Y-r7FtRp-nBlElrsUS7k8rGeRCPkBIy/n/cnyr09qj8zbo/b/england-data/o/out/CP2021051100119/X204SC24110826-Z01-F001/MD5.txt
tar -xvf X204SC24110826-Z01-F001.tar
```
After data downloading (and unzip), I recommend to verify data integrity using `md5sum -c MD5.txt`. Make sure all files are 'OK'.

## 1. Trimming
After we get raw sequencing reads, I usually directly go to reads trimming using `trimmomatic`, since low quality data will be filtered out, and adapters have been removed when raw data released. But you can run `fastqc` to check adapters, compare the qualities between raw reads and `trimmomatic` processed reads. 
  Raw reads were filtered and trimmed to get high quality reads using Trimmomatic (v0.39). Briefly, the bases with quality less than 20 at the start or end of a read were cut off, and reads with length shorter than 60bp were dropped.
Since there are many samples, I use a `for` loop to iterate.
Be aware that two libraries were sequenced for some samples.

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
        # Need to set specific quality filtering criteria in the following
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
        # Need to set specific quality filtering criteria in the following
        java -jar /usr/local/software/trimmomatic/0.39/trimmomatic-0.39.jar PE -threads 8 -summary "${sample_name}.summary" "$fq1" "$fq2" \
        "$output_fwd_paired" "$output_fwd_unpaired" "$output_rev_paired" "$output_rev_unpaired" \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:60
        cd "${output_dir}"
    fi
    echo "${sample_name} finished"
done

echo "Trimmomatic process completed."
```
After this step we get a summary file for each sample. This file reports `number of input reads pairs, number (and percentage) of both forward and reverse reads survived`. I recommend to save these information, and provide these in a supplementary table when prepare manuscript.

Note, if there are too many dropped reads during this step (e.g. >10%), it may suggest your data quality might be low, you need doublecheck.

## 2. Mapping
After we get high quality reads, we proceed to align these reads to reference genome. After testing several mapper programs, I found `novoalign` perform better (e.g. more mapped reads, higher mapping quality) on many datasets (fungi, plant, nematodes). But, I never tested this on polyploid. Again I recommend to test several mapper programs on your own datasets.

 `novoalign` is a commercial program, but it's free for academia, you can request a license from `sales@novocraft.com`.
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
    # FPKM normalized expression
    /data/pathology/program/stringtie-3.0.0.Linux_x86_64/stringtie -p 32 -G /rds-d7/user/cx264/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/01.SpatialTranscriptomics/00.ref/M.truncatula_genome_v4/Mtruncatula_285_Mt4.0v1.gene.gtf -e -B -A 4.sorted.FPKM.tsv 2.sorted.bam

    echo "${sample_name} completed"
    cd "${working_dir}"
done
```
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
                    sample_id = sample_path.split('/')[-1]  # Adjust this if your naming convention differs

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

In this script, I use `qualimap` and `mapinsingts` to check the quality of alignment for each sample. 
In the last part of this script, I use `htseq-count` to quantify the expression (to count number of reads mapped to each gene/transcript). We call these "raw counts", and these raw counts will be used for subsequent differential expression analysis.

## 3. QC

In this pipeline, I do quality control throughout many steps. For example, after `trimmomatic` step, you can check how many reads (forward and/or reverse) survived the trimming criteria. Also in mapping step, I used `qualimap` and `mapinsights` to check quality/summary of alignment results, statistics including 'total reads aligned','Strand ratio','Mean mapping quality', etc. The results are stored in `./qualimap/genome_results.txt` and `./Overall_mapping_summary.log`. You can retrieve alignment statistics (e.g.: number (percentage) of mapped reads; mapping quality) and put in a supplementary table for your manuscript.

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
<table>
  <tr>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/03.Mapping/WT_EV_2/RSeQC/05.WT_EV_2.inner_distance_plot.jpeg" alt="Image 1" width="400"/></td>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/03.Mapping/WT_EV_2/RSeQC/03.WT_EV_2.geneBodyCoverage.curves.jpeg" alt="Image 2" width="400"/></td>
  </tr>
  <tr>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/03.Mapping/WT_EV_2/RSeQC/10.Reads.Distribution.over.GeneBody.jpg" alt="Image 3" width="400"/></td>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/03.Mapping/WT_EV_2/RSeQC/15.WT_EV_2.saturation.jpeg" alt="Image 4" width="400"/></td>
  </tr>
</table>

### 3.2TOmicsVis_Clustering
Following is a way to visualize processed RNAseq results (mainly based on expression counts). Some of the featurs are helpful for quality control, e.g. correlation heatmap of replicates/treatments, and several clustering plots (PCA, tsne).

```r
library(TOmicsVis)
library(readr)

setwd("/data/pathology/cxia/projects/Giles/Jinpeng/X204SC24110826-Z01-F001/05.Clustering")

counts <- read_tsv("../04.DE_analysis/raw.counts.order.tsv")
#please note the replicates need to be listed together for clustering, as the plot only shows in the original order as in the input
group <- read_tsv("../04.DE_analysis/groups.tsv")
group <- as.data.frame(group)

# Correlation Heatmap for samples/groups based on Pearson/Spearman algorithm
pdf("061.corr_heatmap.pearson.pdf")
corr_heatmap(
  data = counts,
  corr_method = "pearson",
  cell_shape = "square",
  fill_type = "full",
#  multi_shape = FALSE,
  lable_size = 2,
  axis_angle = 45,
  axis_size = 9,
  lable_digits = 1,
  color_low = "blue",
  color_mid = "white",
  color_high = "red",
  outline_color = "white",
  ggTheme = "theme_bw"
)
corr_heatmap(
  data = counts,
  corr_method = "spearman",
  cell_shape = "square",
  fill_type = "full",
#  multi_shape = FALSE,
  lable_size = 2,
  axis_angle = 45,
  axis_size = 9,
  lable_digits = 1,
  color_low = "blue",
  color_mid = "white",
  color_high = "red",
  outline_color = "white",
  ggTheme = "theme_bw"
)
corr_heatmap(
  data = counts,
  corr_method = "kendall",
  cell_shape = "square",
  fill_type = "full",
#  multi_shape = FALSE,
  lable_size = 2,
  axis_angle = 45,
  axis_size = 9,
  lable_digits = 1,
  color_low = "blue",
  color_mid = "white",
  color_high = "red",
  outline_color = "white",
  ggTheme = "theme_bw"
)
dev.off()

# pca_analysis
pdf("062.pca.pdf")
pca_plot(
  sample_gene = counts,
  group_sample = group,
  multi_shape = FALSE,
  xPC = 1,
  yPC = 2,
  point_size = 4,
  text_size = 2,
  fill_alpha = 0.10,
  border_alpha = 0.00,
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_bw"
)
dev.off()

# tsne plot
pdf("063.tsne.pdf")
tsne_plot(
  sample_gene = counts,
  group_sample = group,
  seed = 1,
  multi_shape = FALSE,
  point_size = 3,
  point_alpha = 0.8,
  text_size = 5,
  text_alpha = 0.60,
  fill_alpha = 0.20,
  border_alpha = 0.00,
  sci_fill_color = "Sci_AAAS",
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_bw"
)
dev.off()

# UMAP
pdf("064.umap.pdf")
umap_plot(
  sample_gene = counts,
  group_sample = group,
  seed = 1,
  multi_shape = FALSE,
  point_size = 3,
  point_alpha = 1,
  text_size = 5,
  text_alpha = 0.60,
  fill_alpha = 0.20,
  border_alpha = 0.00,
  sci_fill_color = "Sci_AAAS",
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_bw"
)
dev.off()

# dendro_plot
pdf("065.samples_dendrop.pdf")
dendro_plot(
  data = counts,
  dist_method = "euclidean",
  hc_method = "ward.D2",
  tree_type = "rectangle",
  k_num = 5,
  palette = "npg",
  color_labels_by_k = TRUE,
  horiz = FALSE,
  label_size = 1,
  line_width = 1,
  rect = TRUE,
  rect_fill = TRUE,
  xlab = "Samples",
  ylab = "Height",
  ggTheme = "theme_light"
)
dev.off

# Heatmap group for visualizing grouped gene expression data, too busy for whole-gene sets, better for DE genes only!
pdf("065.heatmap_group.pdf")
heatmap_group(
  sample_gene = counts,
  group_sample = group,
  scale_data = "row",
  clust_method = "complete",
  border_show = TRUE,
  border_color = "#ffffff",
  value_show = TRUE,
  value_decimal = 2,
  value_size = 5,
  axis_size = 8,
  cell_height = 10,
  low_color = "#00880055",
  mid_color = "#ffffff",
  high_color = "#ff000055",
  na_color = "#ff8800",
  x_angle = 45
)
heatmap_cluster(
  data = counts,
  dist_method = "euclidean",
  hc_method = "average",
  k_num = 6,
  show_rownames = FALSE,
  palette = "RdBu",
  cluster_pal = "Set1",
  border_color = "#ffffff",
  angle_col = 45,
  label_size = 10,
  base_size = 12,
  line_color = "#0000cd",
  line_alpha = 0.2,
  summary_color = "#0000cd",
  summary_alpha = 0.8
)
dev.off()
```

<table>
  <tr>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/05.Clustering/061.corr_heatmap.pearson.jpeg" alt="Image 1" width="400"/></td>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/05.Clustering/062.pca.jpeg" alt="Image 2" width="400"/></td>
  </tr>
  <tr>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/05.Clustering/062.tsne.jpeg" alt="Image 3" width="400"/></td>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/05.Clustering/065.samples_dendrop.jpeg" alt="Image 4" width="400"/></td>
  </tr>
</table>

## 4. Differential Expression Analysis
For DE analysis, my code is mainly from `edgeR` user's guide (https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf). This user's guide is a must-read in order to know the rationale behind the code.
Differential Expression analysis using `edgeR` under `mamba activate R4.3.2`.

### 4.1 WTEV_vs_WTCC
```R
library(EnhancedVolcano)
library(gplots)
library(RColorBrewer)
library("ggplot2")
library(edgeR)
library(stringr)
library(statmod)
library(readr)
library(ClusterGVis) #clutering
library(TOmicsVis)

my_count_matrix <- read_tsv("../raw.counts.tsv")
my_count_matrix1 <- my_count_matrix[,c(1,17,33,7,34,35)]
group <- c(1,1,2,2,2)
sample <- factor(c("1","2","1","2","3"))
y <- DGEList(counts=my_count_matrix1[,c(2:6)], group=group,genes=my_count_matrix[,1])
#no filtering to keep all genes for all conditions
#keep <- filterByExpr(y, min.count = 0.01, min.total.count = 0.01, large.n = 2, min.prop = 0.01)
keep <- filterByExpr(y, group=group)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)
# [1] 25120     5
y <- normLibSizes(y)
y$samples

pdf("1.1.LRT.MDS.pdf")
plotMDS(y, col=c("1","1","2","2","2"))
dev.off()

# continue from previous section
design <- model.matrix(~sample+group)
rownames(design) <- colnames(y)
design
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
# [1] 0.322674
pdf("1.2.LRT.dispersion.BCV.pdf")
plotBCV(y)
dev.off()
fit <- glmQLFit(y, design, robust=TRUE)
pdf("1.3.LRT.QLdispersion.pdf")
plotQLDisp(fit)
dev.off()
qlf <- glmQLFTest(fit)
topTags(qlf)
cpm(y)[rownames(topTags(qlf)),]
summary(decideTests(qlf,lfc=0.5))
allTags2 <- topTags(qlf, n = nrow(qlf$genes), adjust.method = "BH", sort.by = "none", p.value = 1)
DEG_glmQLF <- as.data.frame(allTags2)
log2FC <- 0.5
k1 <- (DEG_glmQLF$FDR < 0.05) & (DEG_glmQLF$logFC < -0.5)
k2 <- (DEG_glmQLF$FDR < 0.05) & (DEG_glmQLF$logFC > 0.5)
DEG_glmQLF$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "Not_significant"))
table(DEG_glmQLF$change)

DEG_glmQLF$CPM <- cpm(y)[rownames(qlf),]
write.csv(DEG_glmQLF, "1.4.QLF.DE.results.csv", quote = F)

## Visualization
deseq_results <- read.csv("1.4.QLF.DE.results.csv",sep=",", header = T)
diffGenes <- deseq_results[deseq_results$FDR < 0.05,][,c(2,9,10,11,12,13)]
clustRows <- hclust(as.dist(1-cor(t(diffGenes[,c(2,3,4,5,6)]),method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(diffGenes[,c(2,3,4,5,6)],method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9)
module.color <- module.color[as.vector(module.assign)]
#myheatcolors1 <- blueorange(75)
myheatcolors1 <- colorRampPalette(c("blue", "orange"))(25)
pdf("1.5.QLF.DiffGene.heatmap.pdf")
heatmap.2(data.matrix(diffGenes[,c(2,3,4,5,6)]),
    Rowv=as.dendrogram(clustRows),
    Colv=as.dendrogram(clustColumns),
    RowSideColors=module.color,
    col=myheatcolors1, scale='row', labRow=NA,
    density.info="none", trace="none",
    cexRow=1, cexCol=1, margins=c(8,20), keysize=1, srtCol=45)
dev.off()
png("1.5.QLF.DiffGene.heatmap.png",width = 636, height = 980)
heatmap.2(data.matrix(diffGenes[,c(2,3,4,5,6)]),
    Rowv=as.dendrogram(clustRows),
    Colv=as.dendrogram(clustColumns),
    RowSideColors=module.color,
    col=myheatcolors1, scale='row', labRow=NA,
    density.info="none", trace="none",
    cexRow=1, cexCol=1, margins=c(8,20), keysize=0.5, srtCol=45)
dev.off()

pdf("1.6.QLF.DiffGene.volcano.pdf", 7, 7)
EnhancedVolcano(deseq_results, lab = rownames(deseq_results), x = 'logFC', y = 'FDR', pCutoff = 5e-2,
    FCcutoff=1.0, ylim = c(0, 7.5), xlim = c(-10, 15), pointSize = 1.0, labSize = 0,
    colAlpha = 1, legendIconSize=2.0, legendLabSize = 10) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
png("1.6.QLF.DiffGene.volcano.png", 800, 800)
EnhancedVolcano(deseq_results, lab = rownames(deseq_results), x = 'logFC', y = 'FDR', pCutoff = 5e-2,
    FCcutoff=1.0, ylim = c(0, 7.5), xlim = c(-10, 15), pointSize = 1.0, labSize = 0,
    colAlpha = 1, legendIconSize=2.0, legendLabSize = 10) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()

# gene ranking dotplot, function from TOmicsVis
# data should be (1st-col: Genes, 2nd-col: log2FoldChange, 3rd-col: Pvalue, 4th-col: FDR)
pdf("1.7.GeneRanking.pdf")
gene_rank_plot(
  data = deseq_results[,c(2,3,6,7)],
  log2fc = 1,
  palette = "Spectral",
  top_n = 10,
  genes_to_label = NULL,
  label_size = 5,
  base_size = 12,
  title = "Gene ranking dotplot",
  xlab = "Ranking of differentially expressed genes",
  ylab = "Log2FoldChange"
)
dev.off()

q()
```
### 4.2 WTEV_vs_WTRE
```R
library(EnhancedVolcano)
library(gplots)
library(RColorBrewer)
library("ggplot2")
library(edgeR)
library(stringr)
library(statmod)
library(readr)
library(ClusterGVis) #clutering
library(TOmicsVis)

my_count_matrix <- read_tsv("../raw.counts.tsv")
my_count_matrix1 <- my_count_matrix[,c(1,4,12,14,7,34,35)]
group <- c(1,1,1,2,2,2)
sample <- factor(c("1","2","3","1","2","3"))
y <- DGEList(counts=my_count_matrix1[,c(2:7)], group=group,genes=my_count_matrix[,1])
#no filtering to keep all genes for all conditions
#keep <- filterByExpr(y, min.count = 0.01, min.total.count = 0.01, large.n = 2, min.prop = 0.01)
keep <- filterByExpr(y, group=group)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)
# [1] 23484     6
y <- normLibSizes(y)
y$samples

pdf("1.1.LRT.MDS.pdf")
plotMDS(y, col=c("1","1","1","2","2","2"))
dev.off()
```
Following is an MDS plot shows the relative similarities of the six samples. This gives you another way to evaluate your experiment quality.
![MDS plot](https://github.com/chongjing/RNAseq_Medicago/blob/main/04.DE_analysis/02.WTEV_vs_WTRE/1.1.LRT.MDS.jpeg)

```R
design <- model.matrix(~sample+group)
rownames(design) <- colnames(y)
design
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
# [1] 0.04480926
pdf("1.2.LRT.dispersion.BCV.pdf")
plotBCV(y)
dev.off()
fit <- glmQLFit(y, design, robust=TRUE)
pdf("1.3.LRT.QLdispersion.pdf")
plotQLDisp(fit)
dev.off()
qlf <- glmQLFTest(fit)
topTags(qlf)
cpm(y)[rownames(topTags(qlf)),]
summary(decideTests(qlf,lfc=0.5))
allTags2 <- topTags(qlf, n = nrow(qlf$genes), adjust.method = "BH", sort.by = "none", p.value = 1)
DEG_glmQLF <- as.data.frame(allTags2)
log2FC <- 0.5
k1 <- (DEG_glmQLF$FDR < 0.05) & (DEG_glmQLF$logFC < -0.5)
k2 <- (DEG_glmQLF$FDR < 0.05) & (DEG_glmQLF$logFC > 0.5)
DEG_glmQLF$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "Not_significant"))
table(DEG_glmQLF$change)
#           DOWN Not_significant              UP
#           2532           18898            2054

DEG_glmQLF$CPM <- cpm(y)[rownames(qlf),]
write.csv(DEG_glmQLF, "1.4.QLF.DE.results.csv", quote = F)

## Visualization
deseq_results <- read.csv("1.4.QLF.DE.results.csv",sep=",", header = T)
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
    FCcutoff=1.0, ylim = c(0, 7.5), xlim = c(-10, 15), pointSize = 1.0, labSize = 0,
    colAlpha = 1, legendIconSize=2.0, legendLabSize = 10) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
png("1.6.QLF.DiffGene.volcano.png", 800, 800)
EnhancedVolcano(deseq_results, lab = rownames(deseq_results), x = 'logFC', y = 'FDR', pCutoff = 5e-2,
    FCcutoff=1.0, ylim = c(0, 7.5), xlim = c(-10, 15), pointSize = 1.0, labSize = 0,
    colAlpha = 1, legendIconSize=2.0, legendLabSize = 10) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
```

<table>
  <tr>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/04.DE_analysis/02.WTEV_vs_WTRE/1.5.QLF.DiffGene.heatmap.jpeg" alt="Image 1" width="400"/></td>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/04.DE_analysis/02.WTEV_vs_WTRE/1.6.QLF.DiffGene.volcano.png" alt="Image 2" width="400"/></td>
  </tr>
</table>

```R

# gene ranking dotplot, function from TOmicsVis
# data should be (1st-col: Genes, 2nd-col: log2FoldChange, 3rd-col: Pvalue, 4th-col: FDR)
pdf("1.7.GeneRanking.pdf",9,7)
gene_rank_plot(
  data = deseq_results[,c(2,3,6,7)],
  log2fc = 1,
  palette = "Spectral",
  top_n = 10,
  genes_to_label = NULL,
  label_size = 5,
  base_size = 12,
  title = "Gene ranking dotplot",
  xlab = "Ranking of differentially expressed genes",
  ylab = "Log2FoldChange"
)
dev.off()

q()
```
![A gene rank plot](https://github.com/chongjing/RNAseq_Medicago/blob/main/04.DE_analysis/02.WTEV_vs_WTRE/1.7.GeneRanking.jpeg)

## 5. Gene Expression Clustering
`ClusterGVis` documents states normalized matrix or data frame containing gene expressions are accepted, so we get FPKM first, and then use Differentially Expressed genes for clustering trends.
### 5.1 get FPKM

```R
library(GenomicFeatures)
library(dplyr)
library(GenomicRanges)
library(tibble)

raw_counts <- read.csv("raw.counts.tsv", sep="\t", row.names = 1)
gtf_file <- "Mtruncatula_285_Mt4.0v1.gene.gtf"

# Create a TxDb object from the GTF file
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
# Extract gene lengths (in kilobases)
exons_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths <- sapply(exons_by_gene, function(x) sum(width(x)))  # Total exon length per gene in base pairs
gene_lengths <- data.frame(GeneID = names(gene_lengths), gene_length_kb = gene_lengths / 1000)
# Step 3: Calculate total mapped reads per sample
# Total mapped reads is the sum of raw counts for each sample
total_mapped_reads <- colSums(raw_counts) / 1e6  # Convert to millions

# Step 4: Calculate FPKM
fpkm_data <- raw_counts %>%
  rownames_to_column(var = "GeneID") %>%  # Convert row names to a column
  left_join(gene_lengths, by = "GeneID") %>%  # Join with gene lengths
  column_to_rownames(var = "GeneID")  # Convert GeneID column back to row names
# Calculate FPKM for each sample
options(scipen = 999)
fpkm <- sapply(colnames(raw_counts), function(sample) {
  result <- (raw_counts[, sample] / fpkm_data$gene_length_kb) / total_mapped_reads[sample]
  format(result, digits = 5, nsmall = 5)
})
options(scipen = 0)

# Convert to a data frame
fpkm <- as.data.frame(fpkm)
rownames(fpkm) <- rownames(raw_counts)

# Save the FPKM data to a file
write.csv(fpkm, "fpkm_normalized_counts.csv", row.names = TRUE, quote=F)
```
### 5.2 Clustering for expression matrix
```R
library(ClusterGVis)
library(Biobase)
library(Mfuzz)

FPKM_counts <- read.csv("07.FPKM.genes4clustering.csv", sep="\t", row.names = 1)
FPKM_counts <- FPKM_counts[,c(24,26,36,4,19,21,29,1,18,27,30,15,17,20,25)]
#define a suitable cluster numbers
pdf("081.getClusters.pdf")
getClusters(FPKM_counts)
cm <- clusterData(FPKM_counts, cluster.method = "mfuzz", cluster.num = 8)
ct <- clusterData(FPKM_counts, cluster.method = "TCseq", cluster.num = 8)
ck <- clusterData(FPKM_counts, cluster.method = "kmeans", cluster.num = 8)
visCluster(cm, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(ct, plot.type = "line", ms.col = c("green", "orange", "red"))
visCluster(ck, plot.type = "line")
dev.off()


#line plot
pdf("082.line_heatmap.pdf")
visCluster(object = cm,
           plot.type = "both",
           ms.col = c("green","orange","red"),
           column_names_rot = 45)
dev.off()

# save cluster information and membership (if exit)
write.csv(cm$wide.res, "09.cm.8clusters.csv",row.names = TRUE, quote = F)
write.csv(ct$wide.res, "09.ct.8clusters.csv",row.names = TRUE, quote = F)
write.csv(ck$wide.res, "09.ck.8clusters.csv",row.names = TRUE, quote = F)
```
<table>
  <tr>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/05.Clustering/081.getClusters.jpeg" alt="Image 1" width="300"/></td>
    <td><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/05.Clustering/082.ck.cluster.jpeg" alt="Image 2" width="400"/></td>
  </tr>
  <tr>
    <td colspan="2"><img src="https://github.com/chongjing/RNAseq_Medicago/blob/main/05.Clustering/082.line_heatmap.jpeg" alt="Image 3" width="750"/></td>
  </tr>
</table>