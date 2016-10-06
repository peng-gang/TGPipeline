# TGPipeline
Pipeline for Target Sequencing Data Analysis

##  Usage
First download the source file from "src" directory. Then change the content in res.txt according to your computer. res.txt indluces the software and its location used in the pipeline. 
After changing "res.txt", typing the following command to run the pipeline: 

*`Rscript TGPipeline.R -fq /fastq/file/dir/ -o /output/dir/ -build hg38 -trim 20 -ref /reference/sequence.fa -dbsnp /dbsnp/info.vcf -ampbed amplicon.bed -seqbed TargetSequence.bed`*

Check manual for details.

If you have any questions, please email me: michael dot gang dot peng AT gmail dot com
