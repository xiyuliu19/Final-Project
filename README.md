# Final-Project
TRGN 510
Tumor Microsatellite Instability classification of melanoma
***
## Description  
The deliverable of this project is a HTML Notebook which shows the pipelines of transferring FASTQ files to MAF files, and the implement of an existing pipeline—MSIpred for melanoma microsatellite instability classification from tumor mutation annotation data using a support vector machine. What’s more, there are also plots of variation analysis created by R in this project.  
## Datasets  
The dataset for this project is FASTQ files made up of melanoma cases. In the process, reference files are also needed for aligning and annotation.
## Proposed Analysis  
The core process of this project is getting tumor microsatellite instability classification based on FASTQ files.  
To realize it, the first step is to transfer FASTQ files to BAM files using BWA-MEM for aligning.  
Next, using Strelka to call variant, and this step can transfer BAM files to VCF files.  
Then, annotating the VCF files by SNPEFF, and getting MAF files.  
Now, ***the first part of variation analysis can be done in R studio.***  
        1.Library "MAFtools"  
        2.Analyze the variant classification, variant type, SNV class, variants per case, and top 10 mutated genes.  
        3.Create plots, such as waterfall map and gene word cloud.  
***The second part of this project in dealing with tumor microsatellite instability classification is performed by MSIpred.***   
  
## Proposed timeline & major milestones  
Timelines:   
Due to 11/08/2019: Identify the pipeline for tumor microsatellite instability classification and the cancer type for analysis.   
Due to 11/13/2019: Complete the alignment and variant calling part by successfully transferring FASTQ files to VCF files. (Milestone 1)   
Due to 11/20/2019: Complete the annotation part by successfully transferring VCF files to MAF files. And finish the variant analysis work by R. (Milestone 2)   
Due to 11/27/2019: Complete the microsatellite instability classification. (Milestone 3)   
Due to 12/05/2019: Complete the HTML Notebook. Post the final project.   
***
## Progress Update   
### 11/12/2019   
#### Dataset   
Original_Data=/home/xiyuliu/unitTest_data   
gunzip   
FASTQ_Directory=/home/xiyuliu/source   
For normal:   
R1=normal_L002_R1_001.fastq   
R2=normal_L002_R2_001.fastq   
For tumor:   
R1=tumor_L001_R1_001.fastq   
R2=tumor_L001_R2_001.fastq   
#### Reference Data   
**For alignment**:   
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz   
gunzip hg38.fa.gz   
Ref_ali=/home/xiyuliu/hg38.fa   
**For annotation**:   
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz   
gunzip GRCh38.primary_assembly.genome.fa.gz   
Ref_ann=/home/xiyuliu/GRCh38.primary_assembly.genome.fa   
#### Software installation   
From fastq to bam:  **BWA-MEM**   
git clone https://github.com/lh3/bwa.git   
cd bwa  
make  
./bwa index hg38.fa  
./bwa mem hg38.fa normal_L002_R1_001.fastq normal_L002_R2_001.fastq | gzip -3 > normal-aln-pe.sam.gz   
./bwa mem hg38.fa tumor_L001_R1_001.fastq tumor_L001_R2_001.fastq | gzip -3 > tumor-aln-pe.sam.gz   
