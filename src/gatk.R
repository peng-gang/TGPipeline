gatk <- function(bamfile, ref,  bam.dir, out.dir, res.info, dbsnp = NULL, bedfile = NULL){
  # bamfile: bam file name
  # ref: reference file
  # bam.dir: bam file directory
  # out.dir: output directory
  # res.info: resource information
  # dbsnp: dbsnp file
  # bedfile: bed file for target sequence region
  
  #variatn calling
  out.vcf <- paste(out.dir, substr(bamfile, 1, nchar(bamfile) - 3), "g.vcf", sep = "")
  cmd <- paste0("java -Xmx4g -jar ", res.info["gatk"], " -T HaplotypeCaller ")
  cmd <- paste(cmd, "--genotyping_mode DISCOVERY ", sep = "")
  cmd <- paste(cmd, "-R ", ref, sep = "")
  cmd <- paste(cmd, " -I ", bam.dir, bamfile, sep = "")
  cmd <- paste(cmd, " -o ", out.vcf, sep = "")
  if(!is.null(dbsnp)){
    cmd <- paste(cmd, " --dbsnp ", dbsnp, sep = "")
  }
  if(!is.null(bedfile)){
    cmd <- paste(cmd, " -L ", bedfile);
  }
  cmd <- paste(cmd, " -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --maxReadsInRegionPerSample 15000")
  
  message(cmd)
  system(cmd)
}


gatkParallel <- function(bam.dir, ref, out.dir, res.info, dbsnp = NULL, bedfile = NULL){
  # bam.dir: bam file directory
  # ref: reference file
  # out.dir: output directory
  # res.info: resorce information
  # dbsnp: dbsnp file
  # bedfile: bed file for target seqeunce region
  X <- list.files(path = bam.dir, pattern = "*.bam$")
  dir.raw <- paste0(out.dir, "raw/")
  cmd <- paste0("mkdir ", dir.raw)
  system(cmd)
  mclapply(X, gatk, ref, bam.dir, dir.raw, res.info, dbsnp, bedfile, mc.cores = 8)
  
  cmd <- paste0("java -Xmx4g -jar ", res.info["gatk"], " -T GenotypeGVCFs ")
  cmd <- paste0(cmd, "-R ", ref)
  if(!is.null(dbsnp)){
    cmd <- paste(cmd, " --dbsnp ", dbsnp, sep = "")
  }
  if(!is.null(bedfile)){
    cmd <- paste(cmd, " -L ", bedfile);
  }
  for(i in 1:length(X)){
    in.vcf <- paste0(dir.raw, substr(X[[i]], 1, nchar(X[[i]]) - 3), "g.vcf")
    cmd <- paste0(cmd, " -V ", in.vcf)
  }
  
  cmd <- paste0(cmd, " -o ", dir.raw, "all.vcf")
  message(cmd)
  system(cmd)
  
  cmd <- paste0("rm -f ", dir.raw, "*.g.vcf")
  message(cmd)
  system(cmd)
  cmd <- paste0("rm -f ", dir.raw, "*.g.vcf.idx")
  message(cmd)
  system(cmd)
  
  for(i in 1:length(X)){
    cmd <- paste0("java -Xmx4g -jar ", res.info["gatk"], " -T SelectVariants ")
    cmd <- paste0(cmd, "-R ", ref, " -V ", dir.raw, "all.vcf ")
    cmd <- paste0(cmd, "-o ", paste0(out.dir, substr(X[[i]], 1, nchar(X[[i]]) - 3), "vcf "))
    cmd <- paste0(cmd, "-sn ", paste0(substr(X[[i]], 1, nchar(X[[i]]) - 4)))
    message(cmd)
    system(cmd)
  }
}



