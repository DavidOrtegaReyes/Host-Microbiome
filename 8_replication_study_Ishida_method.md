```bash
##Ishida et al. post 16S rRNA taxonomic classification method

##Files
#Sample-1
F-CS-20231201-1400-I-1_20231207-S29_L001_R1_001.fastq.gz
F-CS-20231201-1400-I-1_20231207-S29_L001_R2_001.fastq.gz
#Sample-2
F-RK-20231202-2350-I-1_20231207-S31_L001_R1_001.fastq.gz
F-RK-20231202-2350-I-1_20231207-S31_L001_R2_001.fastq.gz
#Sample-3
F-SH-20231124-0030-I-1_20231207-S33_L001_R1_001.fastq.gz
F-SH-20231124-0030-I-1_20231207-S33_L001_R2_001.fastq.gz
#Sample-4
F-TK-20231124-1400-I-1_20231207-S35_L001_R1_001.fastq.gz
F-TK-20231124-1400-I-1_20231207-S35_L001_R2_001.fastq.gz
#Sample-5
F-YO-20231201-1500-I-1_20231207-S37_L001_R1_001.fastq.gz
F-YO-20231201-1500-I-1_20231207-S37_L001_R2_001.fastq.gz

## Use the Bowtie-2 program to remove reads mapped to the PhiX 174 sequence and the Genome Reference Consortium human build 38 (GRCh38).
# Build index for reference files
bowtie2-build PhiX174.fasta PhiX174
bowtie2-build hg38.fa hg38
# Remove reads mapped to the PhiX 174 sequence and the Genome Reference Consortium human build 38 (GRCh38)
bowtie2 -x PhiX174 -U F-YO-20231201-1500-I-1_20231207-S37_L001_R2_001.fastq.gz --un F-YO-T-2.fastq 
bowtie2 -x hg38 -U F-YO-T-2.fastq --un F-YO-T-2_removeHuman.fastq

## Trim the 3' region of each read with a PHRED quality score less than 17. Remove trimmed reads less than 150 bp in length with an average quality score less than 25 
#TRAILING:17: This option removes bases from the end of each read if their quality falls below a specified threshold. In this case, any bases with a quality score below 17 will be trimmed. MINLEN:150: This option discards reads that are shorter than the specified length. In this case, any reads shorter than 150 bases will be discarded. AVGQUAL:25: This option discards reads if their average quality score falls below the specified threshold. In this case, any reads with an average quality score below 25 will be discarded.
java -jar trimmomatic-0.39.jar SE -phred33 F-YO-T-2_removeHuman.fastq F-YO-T-2_removeHuman_Trim.fastq HEADCROP:0 TRAILING:17 MINLEN:150 AVGQUAL:25

## Remove potential chimeric sequences by reference-based chimaera checking in USEARCH (ver. 5.2.32) and the Genomes OnLine Database (GOLD)
#Chimeric sequences are artificial DNA or RNA sequences that are formed by the merging of two or more different sequences during the sequencing process. These chimeras can arise due to various factors, such as PCR artifacts, sequencing errors, or contamination.
vsearch --fastq_filter F-YO-T-2_removeHuman_Trim.fastq --fastaout F-YO-T-2_removeHuman_Trim.fasta
vsearch --uchime_ref F-YO-T-2_removeHuman_Trim.fasta --db gold.fa --nonchimeras F-YO-T-2_removeHuman_Trim_nonchimeric.fasta


###QIIME
#QIIME is a widely used bioinformatics tool for analyzing microbial communities through high-throughput sequencing data.

## Analyze the non-chimeric sequences using the QIIME software 
# Assign the sequences to operational taxonomic units (OTUs) at the genus level using open-reference OTU picking with a 97% pairwise identity threshold and the Greengenes reference database
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
conda activate qiime2-amplicon-2023.9

#'FeatureData[Sequence]' This type of data typically represents biological sequences, such as DNA or RNA sequences, that are used to characterize microbial communities. 'SampleData[SequencesWithQuality]' represents sequencing data that includes both the sequence reads and their associated quality scores.
#Greengenes reference
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path gg_12_10.fasta \
  --output-path greengenes_reference_sequences.qza

#Target sequences from each sample
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path F-CS-T-1_removeHuman_Trim_nonchimeric.fasta \
  --output-path F-CS-T-1.qza

vsearch --fasta2fastq F-CS-T-1_removeHuman_Trim_nonchimeric.fasta --fastq_asciiout 33 --fastqout F-CS-T-1_removeHuman_Trim_nonchimeric.fastq 
gzip F-RK-T-2_removeHuman_Trim_nonchimeric.fastq

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path SampleData.qza \
  --input-format SingleEndFastqManifestPhred33

#Perform denoising on 16S rRNA gene sequencing data: remove sequencing errors and artifacts. The output includes representative sequences, a feature table, and per-sample statistics.
qiime deblur denoise-16S \
  --i-demultiplexed-seqs SampleData.qza \
  --p-trim-length 120 \
  --p-min-reads 1 \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --o-stats deblur-stats.qza

#Perform open-reference clustering of features in a microbial community dataset. It takes the feature table, representative sequences, and reference sequences as input, and clusters the features based on a specified percent identity threshold. The output includes a clustered feature table, clustered representative sequences, and new reference sequences.
qiime vsearch cluster-features-open-reference \
  --i-table table-deblur.qza \
  --i-sequences rep-seqs-deblur.qza \
  --i-reference-sequences greengenes_reference_sequences.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table table-or.qza \
  --o-clustered-sequences rep-seqs-or.qza \
  --o-new-reference-sequences new-ref-seqs-or.qza

##Taxonomic classification
#Taxonomic reference from Greengenes
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path gg_12_10_taxonomy.txt \
  --output-path ref-taxonomy.qza

#fit-classifier-naive-bayes: is training a Naive Bayes classifier in QIIME 2 using the reference sequences and their associated taxonomic assignments. The trained classifier can then be used to classify unknown sequences and assign them taxonomic labels based on their similarity to the reference sequences.
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads new-ref-seqs-or.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifierFINAL.qza

#Classify
qiime feature-classifier classify-sklearn \
  --i-classifier classifierFINAL.qza \
  --i-reads rep-seqs-or.qza \
  --o-classification taxonomyClusALL.qza

qiime metadata tabulate \
  --m-input-file taxonomyClusALL.qza \
  --o-visualization taxonomyClusALL.qzv

#Bar-Plot
qiime taxa barplot \
  --i-table table-or.qza \
  --i-taxonomy taxonomyClusALL.qza \
  --o-visualization taxa-bar-plotsClusALL.qzv