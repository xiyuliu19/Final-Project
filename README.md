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
* Milestone 1   
#### Dataset   
Original_Data=/home/xiyuliu/unitTest_data   
```
gunzip normal_L002_R1_001.fastq.gz
gunzip normal_L002_R2_001.fastq.gz
gunzip tumor_L001_R1_001.fastq.gz
gunzip tumor_L001_R2_001.fastq.gz
```
FASTQ_Directory=/home/xiyuliu/source   
   
For normal:   
R1=normal_L002_R1_001.fastq   
R2=normal_L002_R2_001.fastq   
   
For tumor:   
R1=tumor_L001_R1_001.fastq   
R2=tumor_L001_R2_001.fastq  
   
#### Reference Data   
**For alignment**:   
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz   
gunzip hg38.fa.gz  
```
Ref_ali=/home/xiyuliu/hg38.fa   
   
**For annotation**:   
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz   
gunzip GRCh38.primary_assembly.genome.fa.gz   
```
Ref_ann=/home/xiyuliu/GRCh38.primary_assembly.genome.fa   
   
#### Software Implement   
From fastq to sam:  **`BWA-MEM`**   
```
git clone https://github.com/lh3/bwa.git   
cd bwa  
make  
./bwa index hg38.fa  
./bwa mem hg38.fa normal_L002_R1_001.fastq normal_L002_R2_001.fastq | gzip -3 > normal-aln-pe.sam.gz   
./bwa mem hg38.fa tumor_L001_R1_001.fastq tumor_L001_R2_001.fastq | gzip -3 > tumor-aln-pe.sam.gz   
gunzip normal-aln-pe.sam.gz   
gunzip tumor-aln-pe.sam.gz   
```
From sam to bam: **`samtools`**   
```
samtools view -b -S normal-aln-pe.sam > normal-aln-pe.bam   
samtools view -b -S tumor-aln-pe.sam > tumor-aln-pe.bam   
```
***
In the next 7 days:   
1. Running Strelka to make the bam files to VCF files.    
2. Running vcf2maf to make the VCF files to MAF files.    
3. Finishing the variant analysis part using maftools within R.    
***

From bam to vcf: **`Strelka`**   
Installation   
```
wget https://github.com/Illumina/strelka/releases/download/v2.8.2/strelka-2.8.2.centos5_x86_64.tar.bz2   
tar xvfj strelka-2.8.2.centos5_x86_64.tar.bz2   
```
Preparation: **`samtools`**   
```
samtools sort normal-aln-pe.bam > normal-sorted.bam   
samtools sort tumor-aln-pe.bam > tumor-sorted.bam   
samtools index normal-sorted.bam   
samtools index tumor-sorted.bam   
```
Run   
```
/home/xiyuliu/bwa/strelka-2.8.2.centos5_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam normal-sorted.bam --tumorBam tumor-sorted.bam --referenceFasta hg38.fa --runDir demo_somatic --exome   
./home/xiyuliu/bwa/strelka-2.8.2.centos5_x86_64/bin/demo_somatic/runWorkflow.py -m local -j 8   
gunzip somatic.indels.vcf.gz   
gunzip somatic.snvs.vcf.gz   
```
***
Annotation 1：**`SnpEff`**  
Environment:   
**`java`**   
Installation   
```
wget https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip   
unzip snpEff_latest_core.zip   
```
Download reference database:   
```
java -jar snpEff.jar download GRCh37.75   
```
Run:   
```
java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh37.75 somatic.indels.vcf > somatic.indels_snpeff.vcf   
java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh37.75 somatic.snvs.vcf > somatic.snvs_snpeff.vcf   
```
Output files:    
somatic.indels_snpeff.vcf   
somatic.snvs_snpeff.vcf   
***
Annotation 2: **`oncotator`**   
Environment:   
**`Python`**  
Installation   
```
wget https://github.com/broadinstitute/oncotator/archive/v1.9.2.0.tar.gz   
pip install v1.9.2.0.tar.gz --user   
```
Download reference database:   
```
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/oncotator/oncotator_v1_ds_April052016.tar.gz   
tar zxvf oncotator_v1_ds_April052016.tar.gz
```
Run:   
```
Oncotator -v --db-dir oncotator_v1_ds_April052016 --input_format=VCF --output_format=TCGAMAF somatic.indels.vcf somatic_indels.maf hg19   
Oncotator -v --db-dir oncotator_v1_ds_April052016 --input_format=VCF --output_format=TCGAMAF somatic.snvs.vcf somatic_snvs.maf hg19   
```
Output files:   
somatic_indels.maf   
somatic_snvs.maf   
***
Annotation 3: **`Annovar`**   
Environment:   
**`Perl`**    
Installation   
Download requires xiyuliu@usc.edu email registration   
Get the link to download ANNOVAR via email   
```
tar xvfz annovar.latest.tar.gz   
cd annovar   
```
Download reference database:   
```
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/   
```
Input file format conversion   
```
./convert2annovar.pl -format vcf4old somatic.indels.vcf > somatic.indels.avinput   
./convert2annovar.pl -format vcf4old somatic.snvs.vcf > somatic.snvs.avinput   
```
Run
```
./annotate_variation.pl -geneanno -dbtype refGene -out somatic_indels -build hg38 somatic.indels.avinput humandb/   
./annotate_variation.pl -geneanno -dbtype refGene -out somatic_snvs -build hg38 somatic.snvs.avinput humandb/ 
```
 ***
