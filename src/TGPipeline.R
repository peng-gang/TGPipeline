# Target Sequencing Pipeline

library("parallel")
library("ggplot2")

source("arguments.R")
source("trim.R")
source("alignment.R")
source("sortedBam.R")
source("gatk.R")
source("QualityControl.R")
source("freebayes.R")
source("annovar.R")

# arguments
args <- commandArgs(TRUE)
args.list <- parserArgs(args)


#resource file
if(is.null(args.list[["res"]])){
  res <- "res.txt"
}else{
  res <- args.list[["res"]]
}

# fastq file direcotry
if(is.null(args.list[["fq"]])){
  stop("Please set the directory of fastq files.")
} else {
  fastq.dir <- args.list[["fq"]]
}

# output directory
if(is.null(args.list[["o"]])){
  stop("Please set the output directory.")
} else {
  out.dir <- args.list[["o"]]
}

if(is.null(args.list[["ref"]])){
  stop("Please set the reference file.")
} else {
  ref <- args.list[["ref"]]
}

dbsnp <- args.list[["dbsnp"]]


# trim length, no trimming if not set or set to 0
if(is.null(args.list[["trim"]])){
  trim.length <- 0
} else {
  trim.length <- as.integer(args.list[["trim"]])
}

# bed file for each amplicon (overlap between amplicons)
amp.bed <- args.list[["ampbed"]]
# bed file for target sequencing (no overlap)
seq.bed <- args.list[["seqbed"]]

if(is.null(args.list[["build"]])){
  build <- "hg38"
} else {
  build <- args.list[["build"]]
}



# read resource file
res.info.tmp <- read.table(res, stringsAsFactors = FALSE, sep = "\t")
res.info <- res.info.tmp[,2]
names(res.info) <- res.info.tmp[,1]


# trimming
if(trim.length > 0){
  in.dir <- paste0(out.dir, "trimmed/")
  cmd <- paste0("mkdir ", in.dir)
  system(cmd)
  
  message(" ")
  message("*********************")
  trimfastqParallel(fastq.dir, in.dir, trim.length, res.info)
} else {
  in.dir <- fastq.dir
}


# alignment
message(" ")
message("*********************")
message("Reads Alignment")
fastq.file <- list.files(path = in.dir, pattern = "*.fastq.gz$")
alignment.dir <- paste0(out.dir, "alignment/")
cmd <- paste0("mkdir ", alignment.dir)
system(cmd)
alignmentParallel(fastq.file, ref, in.dir, alignment.dir, res.info)


# sort reads and created bam file
message(" ")
message("*********************")
message("Reads Sorting and Create Bam File (indexed)")
sortedbamParallel(alignment.dir, alignment.dir, res.info)


# QC
message(" ")
message("*********************")
message("Quality Control")
bamfile.qc <- list.files(path = alignment.dir, pattern = "*.bam$")
sample.name.qc <- substr(bamfile.qc, 1, nchar(bamfile.qc)-4)

qc.dir <- paste0(out.dir, "qc/")
cmd <- paste0("mkdir ", qc.dir)
system(cmd)
basecoverage(alignment.dir, seq.bed, qc.dir, sample.name.qc, res.info)
mapsummary(alignment.dir, qc.dir, sample.name.qc, res.info)
uniformity(alignment.dir, amp.bed, qc.dir, sample.name.qc, res.info)

# GATK
message(" ")
message("*********************")
message("Variant call with GATK")
gatk.dir <- paste0(out.dir, "gatk/")
cmd <- paste0("mkdir ", gatk.dir)
system(cmd)
gatkParallel(alignment.dir, ref, gatk.dir, res.info, dbsnp = dbsnp, bedfile = seq.bed)

# FreeBayes
message(" ")
message("*********************")
message("Variant call with freebayes")
freebayes.dir <- paste0(out.dir, "freebayes/")
cmd <- paste0("mkdir ", freebayes.dir)
system(cmd)
freebayesParallel(alignment.dir,ref, freebayes.dir, res.info,  bedfile = seq.bed)

# ANNOVAR
message(" ")
message("*********************")
annovar.dir <- paste0(out.dir, "annovar/")
cmd <- paste0("mkdir ", annovar.dir)
system(cmd)

annovar(gatk.dir, freebayes.dir, annovar.dir, res.info, "annovar", "hg38",
        gatk.pattern=".vcf$", freebayes.pattern = ".vcf$")




