# Final-Project
TRGN 510
Implement a tumor mutation burden pipeline
***
## Description  
The deliverable of this project is a Shiny Server. The user can query the tumor mutation burden value of each case based on their Sample ID. What's more, the user can further query gene that is highly expressed in the TMB high group.  
## Datasets  
There are two datasets used in this project. They are all downloaded from the TCGA database Repository. The format of the first dataset is MAF file. The data category is simple nucleotide variation. The data type is masked somatic mutation. The experimental strategy is WXS. Based on the evalution of the pipeline of transfering FASTQ file to VCF file, I choose MuTect2 Variant Aggregation and Masking workflow. My research area is brain cancer and the dataset is from the project of TCGA-LGG which consists 513 cases. The format of the second dataset is txt file, which contains the gene expression information based on the same project TCGA-LGG. The data category is transcriptome profiling, and data type is gene expression quantification. The experimental strategy is RNA-seq, and the workflow type is HTSeq - FPKM.  
The first dataset can be downloaded in this URI: https://portal.gdc.cancer.gov/files/1e0694ca-fcde-41d3-9ae3-47cfaf527f25  
The second dataset can be downloaded in this URI: https://portal.gdc.cancer.gov/repository?facetTab=files&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22brain%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-LGG%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22HTSeq%20-%20FPKM%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22transcriptome%20profiling%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22maf%22%2C%22txt%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Gene%20Expression%20Quantification%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D&searchTableTab=files
## Proposed Analysis  
The major part of this project is the implement of a TMB pipeline to calulate the TMB value for each case. As there are no golden standard for the measurement of TMB until now, there are only a few commercial pipelines available on the market for users to apply, such as TruSight Oncology 500 of illumina and FoundationOne CDx of Foundation Medicine, which we have no access to use it by ourselves. Therefore, I explored an open source pipeline based on the paper "XX" to calulate TMB of each case.  
After I get the further TMB file contains the TMB value for each case, I want to analyze the data in two levels.  
1. variantions (based on the MAF file): In R studio, library "MAFtools", analyze the variant classification, variant type, SNV class, variants per case, and top 10 mutated genes. Create plots, such as waterfall map and gene word cloud.  
2. Relationship between TMB and gene expression （based on the further TMB file and the txt file）: This part is concentrate on differential analysis. In R studio, library "limma", analyze the gene expression in different TMB groups and output the differences in all genes. Create heat map.  
## Proposed timeline & major mialstones  
Timelines:    
Due to 11/06/2018: Identify the original dataset and the major things I want to analysis based on the TMB value.  
Due to 11/13/2018: Identify the pipeline used to calculate the TMB value. Complete the dataset organization and complete the first part of the analysis based on the MAF file. (Milestone 1)  
Due to 11/20/2018: Complete TMB calculation. (Milestone 2)  
Due to 11/27/2018: Complete the analysis of gene expression part. Integrate all of the results and plots. (Milestone 3)  
Due to 12/05/2018: Complete the Shiny Server. Post the final project.  
## User Interface  
Search bar and Plots  
