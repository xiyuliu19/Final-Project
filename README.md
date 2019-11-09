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
1.	Library "MAFtools"  
2.Analyze the variant classification, variant type, SNV class, variants per case, and top 10 mutated genes.  
3.	Create plots, such as waterfall map and gene word cloud.  
***The second part of this project in dealing with tumor microsatellite instability classification is performed by MSIpred.***   
  
## Proposed timeline & major milestones  
Timelines:   
Due to 11/08/2018: Identify the pipeline for tumor microsatellite instability classification and the cancer type for analysis.   
Due to 11/13/2018: Complete the alignment and variant calling part by successfully transferring FASTQ files to VCF files. (Milestone 1)   
Due to 11/20/2018: Complete the annotation part by successfully transferring VCF files to MAF files. And finish the variant analysis work by R. (Milestone 2)   
Due to 11/27/2018: Complete the microsatellite instability classification. (Milestone 3)   
Due to 12/05/2018: Complete the HTML Notebook. Post the final project.   
