trimfastq <- function(filename, fastq.dir, out.dir, trim.length, res.info){
  #filename: file name for fastq file
  #fastq.dir: directory of fastq file
  #out.dir: output directory for trimmed fastq file
  #trim.length: number of bps to trim off
  #res.info: resource information
  
  #Q33: illumina
  trim.length <- trim.length + 1
  cmd <- paste("gunzip -c ", fastq.dir, filename, 
               " | ", res.info["fastx_trimmer"], " -Q33 -f ", 
               sep = "")
  cmd <- paste(cmd, trim.length, " -z -o ", sep = "")
  trim.filename <- gsub("fastq.gz", "trim.fastq.gz", filename)
  cmd <- paste(cmd, out.dir, trim.filename, sep = "")
  message(cmd)
  system(cmd)
}

trimfastqParallel <- function(fastq.dir, out.dir, trim.length, res.info){
  X  <- list.files(path = fastq.dir, pattern = "*.fastq.gz")
  mclapply(X, trimfastq, fastq.dir, out.dir, trim.length, res.info, mc.cores = 8)
  return(NULL)
}
