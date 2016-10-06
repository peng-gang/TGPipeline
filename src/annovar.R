annovarInput <- function(gatk.vcf.dir, freebayes.vcf.dir, out.dir, res.info, out.label, 
                    gatk.pattern=".vcf$", freebayes.pattern = ".vcf$"){
  #create input file for annovar
  
  gatk.vcf.file <- list.files(gatk.vcf.dir, gatk.pattern)
  pos.all <- NULL
  pos.gatk <- list()
  for(file in gatk.vcf.file){
    cmd <- paste0("perl ", res.info["annovar"], "convert2annovar.pl -format vcf4 ", 
                  gatk.vcf.dir, file)
    message(cmd)
    tmp <- system(cmd, intern = TRUE)
    tmp <- do.call(rbind, strsplit(tmp, "\t"))[,1:5]
    
    pos.gatk[[file]] <- tmp
    pos.all <- rbind(pos.all, tmp)
  }
  
  freebayes.vcf.file <- list.files(freebayes.vcf.dir, freebayes.pattern)
  pos.freebayes <- list()
  for(file in freebayes.vcf.file){
    cmd <- paste0("perl ", res.info["annovar"], "convert2annovar.pl -format vcf4 ", 
                  freebayes.vcf.dir, file)
    message(cmd)
    tmp <- system(cmd, intern = TRUE)
    tmp <- do.call(rbind, strsplit(tmp, "\t"))[,1:5]
    pos.freebayes[[file]] <- tmp
    pos.all <- rbind(pos.all, tmp)
  }
  
  unique.pos.all <- unique(pos.all)
  
  chr <- unique.pos.all[,1]
  chr.int <- rep(0, length(chr))
  for(i in 1:length(chr)){
    tmp <- toupper(chr[i])
    if(substr(tmp, 1, 3) == "CHR"){
      tmp <- substring(tmp, 4)
    }
    
    if(tmp == "X"){
      chr.int[i] <- 23
    } else if(tmp == "Y"){
      chr.int[i] <- 24
    } else if(tmp == "MT" || tmp == "M"){
      chr.int[i] <- 25
    } else {
      chr.int[i] <- as.integer(tmp)
    }
  }
  
  rlt <- unique.pos.all[order(chr.int, unique.pos.all[,2], unique.pos.all[,3], 
                              unique.pos.all[,4], unique.pos.all[,5]),]
  out.file <- paste0(out.dir, out.label, ".avinput")
  write.table(rlt, file = out.file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  return(list(av.input = out.file, pos.gatk = pos.gatk, pos.freebayes = pos.freebayes))
}


annovar <- function(gatk.vcf.dir, freebayes.vcf.dir, out.dir, res.info, out.label, builder,
                    gatk.pattern=".vcf$", freebayes.pattern = ".vcf$"){
  #use annovar to anotate the variants
  tmp <- annovarInput(gatk.vcf.dir, freebayes.vcf.dir, out.dir, res.info, out.label, 
                                gatk.pattern, freebayes.pattern)
  annovar.input <- tmp[["av.input"]]
  pos.gatk <- tmp[["pos.gatk"]]
  pos.freebayes <- tmp[["pos.freebayes"]]
  
  #annovar
  cmd <- paste0("perl ", res.info["annovar"], "table_annovar.pl ", annovar.input, " ", res.info["annovar"], "humandb")
  cmd <- paste0(cmd, " -buildver ", builder, " -out ", out.dir, out.label , " -remove ")
  cmd <- paste0(cmd, "-protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,avsnp147,exac03,ljb26_all,clinvar_20160302 -operation g,r,r,f,f,f,f,f,f -nastring .")
  
  message(cmd)
  system(cmd)
  
  annovar.output <- paste0(out.label, ".", builder, "_multianno.txt")
  
  annovar.rlt <- read.table(file = paste0(out.dir, annovar.output), header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)
  ind.annovar <- paste0(annovar.rlt[,1], annovar.rlt[,2], annovar.rlt[,3], annovar.rlt[,4], annovar.rlt[,5])
  
  gatk.vcf.file <- list.files(gatk.vcf.dir, gatk.pattern)
  freebayes.vcf.file <- list.files(freebayes.vcf.dir, freebayes.pattern)
  
  n.sample <- length(gatk.vcf.file)
  
  for(i in 1:n.sample){
    g.file <- gatk.vcf.file[i]
    f.file <- freebayes.vcf.file[i]
    message(paste("Checking file", g.file, f.file))
    
    out.file.share <- paste0(out.dir, substr(g.file, 1, nchar(g.file) - 4), "_",
                             substr(f.file, 1, nchar(f.file) - 3), "share.txt")
    out.file.gatk <- paste0(out.dir, substr(g.file, 1, nchar(g.file) - 3), "gakt.txt")
    out.file.freebayes <- paste0(out.dir, substr(f.file, 1, nchar(f.file) - 3), "freebayes.txt")
    
    gatk.pos <- pos.gatk[[g.file]]
    freebayes.pos <- pos.freebayes[[f.file]]
    
    if(is.null(gatk.pos)){
      if(is.null(freebayes.pos)){
        write.table(NULL, file = out.file.share, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        write.table(NULL, file = out.file.gatk, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        write.table(NULL, file = out.file.freebayes, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
      } else if(is.vector(freebayes.pos)) {
        ind.f <- paste0(freebayes.pos[1], freebayes.pos[2], freebayes.pos[3], 
                                freebayes.pos[4], freebayes.pos[5])
        
        write.table(NULL, file = out.file.share, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        write.table(NULL, file = out.file.gatk, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        index <- ind.annovar %in% ind.f
        if(sum(index) != length(ind.f)){
          message("Cannot find the following positions in annovar:")
          message(ind.f[!(ind.f %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.freebayes, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
      } else {
        write.table(NULL, file = out.file.share, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        write.table(NULL, file = out.file.gatk, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        ind.f <- paste0(freebayes.pos[,1], freebayes.pos[,2], freebayes.pos[,3], 
                                freebayes.pos[,4], freebayes.pos[,5])
        index <- ind.annovar %in% ind.f
        if(sum(index) != length(ind.f)){
          message("Cannot find the following positions in annovar:")
          message(ind.f[!(ind.f %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.freebayes, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
      }
    } else if(is.vector(gatk.pos)) {
      if(is.null(freebayes.pos)){
        ind.g <- paste0(gatk.pos[1], gatk.pos[2], gatk.pos[3], gatk.pos[4], gatk.pos[5])
        
        write.table(NULL, file = out.file.share, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        index <- ind.annovar %in% ind.g
        if(sum(index) != length(ind.g)){
          message("Cannot find the following positions in annovar:")
          message(ind.g[!(ind.g %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.gatk, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        write.table(NULL, file = out.file.freebayes, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
      } else if(is.vector(freebayes.pos)) {
        ind.gatk <- paste0(gatk.pos[1], gatk.pos[2], gatk.pos[3], gatk.pos[4], gatk.pos[5])
        ind.freebayes <- paste0(freebayes.pos[1], freebayes.pos[2], freebayes.pos[3], 
                                freebayes.pos[4], freebayes.pos[5])
        
        ind.share <- ind.gatk[ind.gatk %in% ind.freebayes]
        ind.g <- ind.gatk[!(ind.gatk %in% ind.freebayes)]
        ind.f <- ind.freebayes[!(ind.freebayes %in% ind.gatk)]
        
        index <- ind.annovar %in% ind.share
        if(sum(index) != length(ind.share)){
          message("Cannot find the following positions in annovar:")
          message(ind.share[!(ind.share %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.share, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        
        index <- ind.annovar %in% ind.g
        if(sum(index) != length(ind.g)){
          message("Cannot find the following positions in annovar:")
          message(ind.g[!(ind.g %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.gatk, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        index <- ind.annovar %in% ind.f
        if(sum(index) != length(ind.f)){
          message("Cannot find the following positions in annovar:")
          message(ind.f[!(ind.f %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.freebayes, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
      } else {
        ind.gatk <- paste0(gatk.pos[1], gatk.pos[2], gatk.pos[3], gatk.pos[4], gatk.pos[5])
        ind.freebayes <- paste0(freebayes.pos[,1], freebayes.pos[,2], freebayes.pos[,3], 
                                freebayes.pos[,4], freebayes.pos[,5])
        
        ind.share <- ind.gatk[ind.gatk %in% ind.freebayes]
        ind.g <- ind.gatk[!(ind.gatk %in% ind.freebayes)]
        ind.f <- ind.freebayes[!(ind.freebayes %in% ind.gatk)]
        
        index <- ind.annovar %in% ind.share
        if(sum(index) != length(ind.share)){
          message("Cannot find the following positions in annovar:")
          message(ind.share[!(ind.share %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.share, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        
        index <- ind.annovar %in% ind.g
        if(sum(index) != length(ind.g)){
          message("Cannot find the following positions in annovar:")
          message(ind.g[!(ind.g %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.gatk, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        index <- ind.annovar %in% ind.f
        if(sum(index) != length(ind.f)){
          message("Cannot find the following positions in annovar:")
          message(ind.f[!(ind.f %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.freebayes, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
      }
    } else {
      if(is.null(freebayes.pos)){
        ind.g <- paste0(gatk.pos[,1], gatk.pos[,2], gatk.pos[,3], gatk.pos[,4], gatk.pos[,5])
        
        write.table(NULL, file = out.file.share, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        index <- ind.annovar %in% ind.g
        if(sum(index) != length(ind.g)){
          message("Cannot find the following positions in annovar:")
          message(ind.g[!(ind.g %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.gatk, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        write.table(NULL, file = out.file.freebayes, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
      } else if(is.vector(freebayes.pos)) {
        ind.gatk <- paste0(gatk.pos[,1], gatk.pos[,2], gatk.pos[,3], gatk.pos[,4], gatk.pos[,5])
        ind.freebayes <- paste0(freebayes.pos[1], freebayes.pos[2], freebayes.pos[3], 
                                freebayes.pos[4], freebayes.pos[5])
        
        ind.share <- ind.gatk[ind.gatk %in% ind.freebayes]
        ind.g <- ind.gatk[!(ind.gatk %in% ind.freebayes)]
        ind.f <- ind.freebayes[!(ind.freebayes %in% ind.gatk)]
        
        index <- ind.annovar %in% ind.share
        if(sum(index) != length(ind.share)){
          message("Cannot find the following positions in annovar:")
          message(ind.share[!(ind.share %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.share, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        
        index <- ind.annovar %in% ind.g
        if(sum(index) != length(ind.g)){
          message("Cannot find the following positions in annovar:")
          message(ind.g[!(ind.g %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.gatk, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        index <- ind.annovar %in% ind.f
        if(sum(index) != length(ind.f)){
          message("Cannot find the following positions in annovar:")
          message(ind.f[!(ind.f %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.freebayes, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
      } else {
        ind.gatk <- paste0(gatk.pos[,1], gatk.pos[,2], gatk.pos[,3], gatk.pos[,4], gatk.pos[,5])
        ind.freebayes <- paste0(freebayes.pos[,1], freebayes.pos[,2], freebayes.pos[,3], 
                                freebayes.pos[,4], freebayes.pos[,5])
        
        ind.share <- ind.gatk[ind.gatk %in% ind.freebayes]
        ind.g <- ind.gatk[!(ind.gatk %in% ind.freebayes)]
        ind.f <- ind.freebayes[!(ind.freebayes %in% ind.gatk)]
        
        index <- ind.annovar %in% ind.share
        if(sum(index) != length(ind.share)){
          message("Cannot find the following positions in annovar:")
          message(ind.share[!(ind.share %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.share, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        
        index <- ind.annovar %in% ind.g
        if(sum(index) != length(ind.g)){
          message("Cannot find the following positions in annovar:")
          message(ind.g[!(ind.g %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.gatk, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
        
        index <- ind.annovar %in% ind.f
        if(sum(index) != length(ind.f)){
          message("Cannot find the following positions in annovar:")
          message(ind.f[!(ind.f %in% ind.annovar)])
        }
        write.table(annovar.rlt[index,], file = out.file.freebayes, quote = FALSE, sep = "\t", 
                    row.names = FALSE)
      }
    }
  }
}

#args <- commandArgs(TRUE)

#res.info <- c("/home/gp332/soft/annovar/", "lala")
#names(res.info) <- c("annovar", "mm")
  
#annovar(args[1], args[2], args[3], res.info, "annovar", "hg38",
#                    gatk.pattern="vcf$", freebayes.pattern = "vcf$")


