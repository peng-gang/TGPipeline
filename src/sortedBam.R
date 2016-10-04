#sort sam file to bam file

sortedbam <- function(sam.file, sam.dir, out.dir, res.info){
  # sam.file: file name of sam file
  # sam.dir: directory of sam file
  # out.dir: output directory of sorted bam file
  # res.info: resource information
  cmd <- paste0("java -jar ", res.info["picard"], " SortSam I=");
  cmd <- paste(cmd, sam.dir, sam.file, " O=", sep = "")
  bamfile <- paste0(substr(sam.file, 1, nchar(sam.file)-4), ".bam")
  cmd <- paste(cmd, out.dir, bamfile, 
               " SORT_ORDER=coordinate", sep = "")
  message(cmd)
  system(cmd)
  
  # index bam file
  cmd <- paste0(res.info["samtools"], " index ", out.dir, bamfile)
  message(cmd)
  system(cmd)
}

sortedbamParallel <- function(sam.dir, out.dir, res.info){
  X <- list.files(path = sam.dir, pattern = "*.sam$")
  mclapply(X, sortedbam, sam.dir, out.dir, res.info, mc.cores = 8)
}

