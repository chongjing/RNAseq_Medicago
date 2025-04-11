# Alignment using STAR

### 1. Reference

For this specific project, we use concatenated genomes of host (TAIR10) and nematode as reference.

```bash
cd /data/pathology/cxia/projects/0.ref/09.Arabidopsis_TAIR10_Hschatti1.2

# concatenate genomes
cat /data/pathology/cxia/projects/0.ref/2.Arabidopsis_TAIR10/TAIR10_chr_all.fasta /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa > TAIR10_Hschachtii.genome.fasta
# concatenate annotations
cat /data/pathology/cxia/projects/0.ref/2.Arabidopsis_TAIR10/TAIR10_GFF3_genes.gtf /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.annotations.gtf > TAIR10_Hschachtii.gtf

# index
/data/pathology/program/STAR/bin/Linux_x86_64/STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /data/pathology/cxia/projects/0.ref/09.Arabidopsis_TAIR10_Hschatti1.2/TAIR10_Hschachtii_index \
  --genomeFastaFiles /data/pathology/cxia/projects/0.ref/09.Arabidopsis_TAIR10_Hschatti1.2/TAIR10_Hschachtii.genome.fasta \
  --sjdbGTFfile /data/pathology/cxia/projects/0.ref/09.Arabidopsis_TAIR10_Hschatti1.2/TAIR10_Hschachtii.gtf \
  --sjdbOverhang 99

```

### 2. Trimming
```bash
cd /data/pathology/cxia/projects/Sebastian/05.X204SC25031192-Z01-F001/X204SC25031192-Z01-F001/02.CleanData

# we are using trimmomatic.jar from /usr/local/software/trimmomatic/0.39/trimmomatic-0.39.jar
input_dir="/data/pathology/cxia/projects/Sebastian/05.X204SC25031192-Z01-F001/X204SC25031192-Z01-F001/01.RawData"
output_dir="/data/pathology/cxia/projects/Sebastian/05.X204SC25031192-Z01-F001/X204SC25031192-Z01-F001/02.CleanData"

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
        java -jar /data/pathology/program/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -summary "${sample_name}_02.summary" "$fq1" "$fq2" \
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
        java -jar /data/pathology/program/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -summary "${sample_name}.summary" "$fq1" "$fq2" \
        "$output_fwd_paired" "$output_fwd_unpaired" "$output_rev_paired" "$output_rev_unpaired" \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:60
        cd "${output_dir}"
    fi
    echo "${sample_name} finished"
done

echo "Trimmomatic process completed."
```

### 3. Aligment

```bash
#!/bin/bash
# we are using STAR for alignment
input_dir="/data/pathology/cxia/projects/Sebastian/05.RNAseq_X204SC25023739.r3/X204SC25023739-Z01-F001/02.CleanData"
output_dir="/data/pathology/cxia/projects/Sebastian/05.RNAseq_X204SC25023739.r3/X204SC25023739-Z01-F001/03.Mapping"

# Run STAR iterate through each sample
for fq1 in "$input_dir"/*/*_1P.fq.gz; do
    fq2="${fq1/_1P.fq.gz/_2P.fq.gz}" # Get the corresponding paired-end file
    cd "${output_dir}"

    # Extract sample name
    sample_name=$(basename "$fq1" | cut -d'_' -f1-3)

    # Set directory
    mkdir -p "${sample_name}"
    echo "Directory $sample_name created"
    # Change to the newly created directory
    cd "${sample_name}"
    echo "Changed directory to $sample_name"

    # Run STAR
    echo "Processing ${sample_name}"
    /data/pathology/program/STAR/bin/Linux_x86_64/STAR \
      --runThreadN 16 \
      --genomeDir /data/pathology/cxia/projects/0.ref/09.Arabidopsis_TAIR10_Hschatti1.2/TAIR10_Hschachtii_index \
      --readFilesIn ${fq1} ${fq2} \
      --readFilesCommand zcat \
      --outFilterMultimapNmax 1 \
      --outSAMmultNmax 1 \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix ${sample_name}

    # Run mapinsights bamqc
    /data/pathology/program/mapinsights/mapinsights/mapinsights bamqc -r /data/pathology/cxia/projects/0.ref/09.Arabidopsis_TAIR10_Hschatti1.2/TAIR10_Hschachtii.genome.fasta -i ${sample_name}Aligned.sortedByCoord.out.bam -o ./
    # Run QualiMap
    /data/pathology/program/QualiMap/qualimap_v2.3/qualimap rnaseq -bam ${sample_name}Aligned.sortedByCoord.out.bam -gtf /data/pathology/cxia/projects/0.ref/09.Arabidopsis_TAIR10_Hschatti1.2/TAIR10_Hschachtii.gtf --java-mem-size=32G -outdir ./qualimap/ -outformat pdf

    # Get expression raw counts
    /data/pathology/program/Miniforge3/bin/htseq-count --type transcript --counts_output 3.sorted.bam.count.tsv --nprocesses 16 --max-reads-in-buffer 1000000 ${sample_name}Aligned.sortedByCoord.out.bam /data/pathology/cxia/projects/0.ref/09.Arabidopsis_TAIR10_Hschatti1.2/TAIR10_Hschachtii.gtf
    # Get expression FPKM TPM
    /data/pathology/program/stringtie-3.0.0.Linux_x86_64/stringtie -p 32 -G /data/pathology/cxia/projects/0.ref/09.Arabidopsis_TAIR10_Hschatti1.2/TAIR10_Hschachtii.gtf -e -B -A 4.sorted.FPKM.tsv ${sample_name}Aligned.sortedByCoord.out.bam

   cd "${output_dir}"

   echo "${sample_name} finished"
done

echo "STAR alignment process completed."


```

 