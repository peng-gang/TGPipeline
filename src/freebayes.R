freebayes <- function(bamfile, ref,  bam.dir, out.dir, res.info, bedfile = NULL){
  # bamfile: bam file name
  # ref: reference file
  # bam.dir: bam file directory
  # out.dir: output directory
  # res.info: resource information
  # bedfile: bed file for target sequence region
  out.vcf <- paste0(out.dir, substr(bamfile, 1, nchar(bamfile) - 3), "vcf")
  cmd <- paste0(res.info["freebayes"], " -f ", ref, " -b ", bam.dir, bamfile)
  cmd <- paste0(cmd, " --min-coverage 20")
  if(!is.null(bedfile)){
    cmd <- paste0(cmd, " -t ", bedfile)
  }
  cmd <- paste0(cmd, " -v ", out.vcf)
  message(cmd)
  system(cmd)
}

freebayesParallel <- function(bam.dir, ref, out.dir, res.info, bedfile = NULL){
  # bam.dir: bam file directory
  # ref: reference file
  # out.dir: output directory
  # res.info: resorce information
  # bedfile: bed file for target seqeunce region
  X <- list.files(path = bam.dir, pattern = "*.bam$")
  mclapply(X, freebayes, ref, bam.dir, out.dir, res.info, bedfile, mc.cores = 8)
}