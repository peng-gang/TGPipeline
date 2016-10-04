#alignment
pairfastq <- function(fastq.file){
  #get paired fastq files and sample name
  sample.name <- NULL
  file.name <- NULL
  for(i in fastq.file){
    fn.tmp <- strsplit(i, "[.]")[[1]][1]
    sn.tmp <- strsplit(fn.tmp, "_")[[1]]
    sample.name <- c(sample.name, sn.tmp[1])
    fn.tmp <- ""
    for(j in sn.tmp){
      if(j == "R1" || j == "R2"){
        next
      }
      fn.tmp <- paste(fn.tmp, j, sep = "")
    }
    file.name <- c(file.name, fn.tmp)
  }
  fn.unique <- unique(file.name)
  
  pair.file <- NULL
  sp.unique <- NULL
  for(j in fn.unique){
    index = which(file.name == j)
    sp.unique <- c(sp.unique, sample.name[index][1])
    pair.tmp <- fastq.file[index]
    pair.file <- rbind(pair.file, pair.tmp)
  }
  
  return(list(pair.file = pair.file, sample.name = sp.unique))
}

alignment <- function(X, ref, fastq.dir, out.dir, res.info){
  #	X: string vector of sample name, paired fastq file name 1, and paired fastq file name 2
  # ref: reference sequence
  # fastq.dir: directory of fastq files
  # out.dir: output directory of sam files
  # res.info: resource information
  
  readgroup <- shQuote(paste('@RG', paste("ID:",X[1], sep = ""), 
                     paste("SM:", X[1], sep = ""), "PL:Illumina",
                     paste("LB:", X[1], sep = ""), 'PU:unit1', sep = "\\t"))
  cmd <- paste(res.info["bwa"], " mem -M -R", readgroup, ref, 
               paste(fastq.dir, X[2], sep = ""), 
               paste(fastq.dir, X[3], sep = ""),  ">",
               paste(out.dir, X[1], ".sam", sep = ""),
               sep = " ")
  system(cmd)
}

alignmentParallel <- function(fastq.file, ref, fastq.dir, out.dir, res.info){
  #get paired fastq files and sample name
  sn.pf <- pairfastq(fastq.file)
  paired.fastq <- sn.pf$pair.file
  sample.name <- sn.pf$sample.name
  
  X <- list()
  for(i in 1:length(sample.name)){
    X[[i]] <- c(sample.name[i], paired.fastq[i,1], paired.fastq[i,2])
  }
  
  mclapply(X, alignment, ref, fastq.dir, out.dir, res.info, mc.cores = 8)
  
  return(NULL)
}

