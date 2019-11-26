# Final-Project
TRGN 510
Tumor Microsatellite Instability classification of melanoma
***
## Description  
The deliverable of this project is a HTML Notebook which shows the pipelines of transferring FASTQ files to MAF files. It also contains the application process of an existing pipeline—MSIpred for melanoma microsatellite instability classification from tumor mutation annotation data using a support vector machine and the MSI prediction result.   
## Datasets  
The dataset of this project is consist of FASTQ files made up of melanoma cases. In the process, reference files are also need for aligning and annotation.
## Proposed Analysis  
The core process of this project is getting tumor microsatellite instability classification based on FASTQ files.  
To realize it, the first step is to transfer FASTQ files to BAM files using `BWA-MEM` for aligning.  
Next, using `Strelka` to call variant, and this step can transfer BAM files to VCF files.  
Then, annotating the VCF files by `VEP`, and using `VCF2MAF` to get MAF files.  
Lastly, using `MSIpred` to get the MSI prediction result.   
   
## Proposed timeline & major milestones  
Timelines:   
Due to 11/08/2019: Identify the pipeline for tumor microsatellite instability classification and the cancer type for analysis.   
Due to 11/13/2019: Complete the alignment and variant calling part by successfully transferring FASTQ files to VCF files. (Milestone 1)   
Due to 11/20/2019: Complete the annotation part by successfully transferring VCF files to MAF files. (Milestone 2)    
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
2. Running VEP and vcf2maf to make the VCF files to MAF files.       
***
### 11/20/2019   
* Milestone 2    
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
Extracting and generating a VCF file with only "PASS":   
```
awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' somatic.snvs.vcf > somatic_snvs_PASS.vcf   
```
***
Annotation VCF: **`VEP`**   
Environment:   
**`Perl`**    
Installation   
```
git clone https://github.com/Ensembl/ensembl-vep.git   
cd ensembl-vep   
cpanm Try::Tiny  
cpanm DBD::mysql  
curl -L http://cpanmin.us | perl - --notest -l $PERL_PATH LWP::Simple LWP::Protocol::https Archive::Extract Archive::Tar Archive::Zip CGI DBI Time::HiRes --force   
perl INSTALL.pl   
y  #Do you want to install any cache files (y/n)?   
y  #Cache directory /home/xiyuliu/.vep does not exists - do you want to create it (y/n)?   
320  #The following species/files are available; which do you want (can specify multiple separated by spaces or 0 for all)?   
y  #Do you want to install any FASTA files (y/n)?   
85  #FASTA files for the following species are available; which do you want (can specify multiple separated by spaces, "0" to install for species specified for cache download)?   
y  #Do you want to install any plugins (y/n)?   
y  #Plugins directory /home/xiyuliu/.vep/Plugins does not exists - do you want to create it (y/n)?   
0  #The following plugins are available; which do you want (can specify multiple separated by spaces or 0 for all)?   
```
From VCF to MAF: **`vcf2maf`**   
Installation   
Downloading the source code of the latest stable release from https://github.com/mskcc/vcf2maf/releases.   
Uploading it into the sever.
```
scp /Users/Liu/vcf2maf-1.6.17.zip xiyuliu@trgn.usc.edu:/home/xiyuliu/vcf2maf/
#Connecting the sever
cd vcf2maf
gunzip vcf2maf-1.6.17.zip
cd vcf2maf-1.6.17
```
Preparation   
```
wget ftp://ftp.broadinstitute.org/pub/ExAC_release/current/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz   
```
Run   
```
perl vcf2maf.pl --input-vcf somatic_snvs_PASS.vcf --output-maf somatic_snvs_PASS.maf --ref-fasta hg38.fa --filter-vcf ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --vep-path ../ensembl-vep --vep-data ~/.vep --ncbi-build GRCh38   
```
Checking the maf file: finding the pathogenic variant (NM_004333.6(BRAF):c.1799T>A (p.Val600Glu))   
```
grep "c.1799T>A" somatic_snvs_PASS.vep.vcf
```
***
####Special section (not necessary)   
My struggling process of trying other annotation methods   
***
Try 1：**`SnpEff`**  
Environment:   
**`java`**   
Installation   
```
wget https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip   
unzip snpEff_latest_core.zip   
```
Download reference database:   
```
java -jar snpEff.jar download GRCh38.86   
```
Run:   
```
java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh38.86 somatic.indels.vcf > somatic_indels_snpeff.vcf   
java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh38.86 somatic.snvs.vcf > somatic_snvs_snpeff.vcf   
```
However, the `vcf2maf` can not ideal with input files annotated by snpEff.    
So, this annotaion method is dropped.   
***
Try 2: **`oncotator`**   
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
However, as it can be seen, the oncotator can only use the reference database of hg19 because it did not update after 2015. So I can not use it because my reference database used for alignment is hg38.   
***
Try 3: **`Annovar`**   
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
However, the `vcf2maf` can not ideal with input files annotated by annovar.    
So, this annotaion method is also dropped.    
***
In the next 7 days:
Running `MSIpred` to get the MSI prediction result.   
***
### 11/26/2019   
* Milestone 3  
Tumor Microsatellite Instability classification of melanoma: **`MSIpred`**   
Install Environment：
```
pip install -U scikit-learn==0.19.1 --user
pip install -U intervaltree==2.1.0 --user   
pip install -U pandas==0.20.3 --user   
```
Install MSIpred package:   
```
python setup.py install   
```
Download Database:   
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz   
```
Preparation of the input MAF file:   
Because in my maf file (somatic_snvs_PASS.maf), there are individual items in the chromosome column that contain some special values other than "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY, I need to filter them first, otherwise an error will be reported.   
Creating a python script (filter.py) to filter out the special values.   
```
vi filter.py   
```
```
#!/usr/bin/python   
   
from __future__ import print_function

maf = "somatic_snvs_PASS.maf"

chrs = {"chr1":1,"chr2":1,"chr3":1,"chr4":1,"chr5":1,"chr6":1,"chr7":1,"chr8":1,"chr9":1,"chr10":1,
		"chr11":1,"chr12":1,"chr13":1,"chr14":1,"chr15":1,"chr16":1,"chr17":1,"chr18":1,"chr19":1,
		"chr20":1,"chr21":1,"chr22":1,"chrX":1,"chrY":1}

with open(maf) as mfin,open("somatic_snvs_PASS_filter_by_chr.maf","w") as outf:
	h1 = mfin.readline().strip()
	h2 = mfin.readline().strip()
	print(h1,file=outf)
	print(h2,file=outf)
	for line in mfin:
		content = line.strip().split("\t")
		if content[4] in chrs:
			print(line.strip(),file=outf)
```
```
chmod 755 filter.py   
./filter.py   
```
Run:  
```
#Run python    
>>> import MSIpred as mp   
>>> snvs_maf = mp.Raw_Maf(maf_path='somatic_snvs_PASS_filter_by_chr.maf')   
>>> snvs_maf.create_tagged_maf(ref_repeats_file='simpleRepeat.txt',tagged_maf_file = 'tagged_snvs.maf')   
>>> tagged_snvs_maf = mp.Tagged_Maf(tagged_maf_path='tagged_snvs.maf')   
>>> snvs_features = tagged_snvs_maf.make_feature_table(exome_size=44)   
>>> predicted_MSI = mp.msi_prediction(feature_table=snvs_features,svm_model=None)   
>>> predicted_MSI.to_csv('MSIpred_prediction.csv')
>>> quit()
```
***
Now, I have the MSI prediction result for my sample. The result is **MSS**.  
In the next 7 days:
1. Creating some plots using R.   
2. Creating the HTML Notebook for posting.   



