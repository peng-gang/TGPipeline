percent <- function(x, digits = 2, format = "f", ...) {
  if (length(x) == 0) 
    return(character())
  
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

# base coverage
basecoverage <- function(bam.dir, seq.bed, out.dir, sample.name, res.info){
  # bam.dir: bam file directory
  # seq.bed: bed file of target sequence region
  # out.dir: output directory
  # sample.name: name of sequenced samples
  # res.info: resource information
  message(" ")
  message("base coverage")
  cmd <- paste(res.info["samtools"], " depth -a -b ", seq.bed, " ", bam.dir, "*.bam", sep = "")
  rdall.outfile <- paste(out.dir, "rd.all.txt", sep = "")
  cmd <- paste(cmd, " > ", rdall.outfile, sep = "")
  message(cmd)
  system(cmd)
  rd.all <- read.table(rdall.outfile, sep = "\t", stringsAsFactors = FALSE)
  
  num.sample <- length(sample.name)
  min.coverage <- rep(0, num.sample)
  over100 <- rep("", num.sample)
  plot100 <- rep(0, num.sample)
  over200 <- rep("", num.sample)
  over500 <- rep("", num.sample)
  over1000 <- rep("", num.sample)
  num.exon <- nrow(rd.all)
  for(i in 1:num.sample){
    min.coverage[i] <- min(rd.all[,i+2])
    over100[i] <- percent(sum(rd.all[,i+2] > 100) / num.exon)
    plot100[i] <- sum(rd.all[,i+2] > 100) / num.exon
    over200[i] <- percent(sum(rd.all[,i+2] > 200) / num.exon)
    over500[i] <- percent(sum(rd.all[,i+2] > 500) / num.exon)
    over1000[i] <- percent(sum(rd.all[,i+2] > 1000) / num.exon)
  }
  
  rlt <- rbind(min.coverage, over100, over200, over500, over1000)
  colnames(rlt) <- sample.name
  outfile <- paste(out.dir, "basecoverage.txt", sep = "")
  write.table(rlt, file = outfile, sep = "\t", quote = FALSE,
              row.names = TRUE, col.names = TRUE)
  
  
  bc <- plot100 * 100
  bc.1 <- data.frame(y = bc, x = 1:num.sample)
  gp <- ggplot(data = bc.1, aes(x = x, y = y))
  gp <- gp + geom_point(colour="#FF9999")
  gp <- gp + labs(title = "Base Pair Coverage", x = "Samples", y = "% base pairs with 100x coverage")
  #gp <- gp + scale_x_continuous(breaks = 1:num.sample, labels = sample.name)
  #gp <- gp + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  gp <- gp + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  gp <- gp + theme(axis.title = element_text(size = rel(2.2)))
  gp <- gp + theme(axis.text = element_text(size = rel(2.2)))
  gp <- gp + theme(plot.title = element_text(size = rel(2.5)))
  index <- which(bc < 75)
  if(length(index) > 0){
    text <- data.frame(x = index , y = bc[index] - 3, text = sample.name[index])
    gp <- gp + geom_text(data = text, aes(x = x, y = y, label = text), size = 7)
  }
  
  wd <- num.sample * 0.2
  if(wd < 7)
    wd <- 7
  if(wd > 14)
    wd <- 14
  png(paste(out.dir, "BasePairCoverage.png", sep = ""), width=wd, height=7, res=120, units = "in")
  #pdf(paste(path, "BasePairCoverage.pdf", sep = ""), width = 14)
  print(gp)
  dev.off()
  
  #cmd <- paste("rm ", rdall.outfile)
  #system(cmd)
}


# summary
mapsummary <- function(bam.dir, out.dir, sample.name, res.info){
  # bam.dir: bam file directory
  # out.dir: output directory
  # sample.name: sample.name: name of sequenced samples
  # res.info: resource information
  message(" ")
  message("mapping summary")
  bamfile <- list.files(path = bam.dir, pattern = "*.bam$")
  num.sample <- length(sample.name)
  
  QCPassed <- rep(0, num.sample)
  Mapped <- rep("", num.sample)
  ProperlyPaired <- rep("", num.sample)
  for(i in 1:num.sample){
    #total reads with secondary
    cmd <- paste(res.info["samtools"], " view -c ", paste(bam.dir, bamfile[i], sep = ""))
    message(cmd)
    num.total <- as.numeric(system(cmd, intern = TRUE))
    
    #mapped reads
    cmd <- paste(res.info["samtools"], " view -c -F 4 ", paste(bam.dir, bamfile[i], sep = ""))
    message(cmd)
    num.mapped <- as.numeric(system(cmd, intern = TRUE))
    
    #total reads without secondary
    cmd <- paste(res.info["samtools"], " view -c -F 256 ", paste(bam.dir, bamfile[i], sep = ""))
    message(cmd)
    num.mapped.unique <- as.numeric(system(cmd, intern = TRUE))
    
    #properly paired
    cmd <- paste(res.info["samtools"], " view -c -f 2 -F 256 ", paste(bam.dir, bamfile[i], sep = ""))
    message(cmd)
    num.paried <- as.numeric(system(cmd, intern = TRUE))
    
    QCPassed[i] <- num.total
    Mapped[i] <- percent(num.mapped/num.total)
    ProperlyPaired[i] <- percent(num.paried/num.mapped.unique)
  }
  
  rlt <- rbind(QCPassed, Mapped, ProperlyPaired)
  colnames(rlt) <- sample.name
  outfile <- paste(out.dir, "summary.txt", sep = "")
  write.table(rlt, file = outfile, sep = "\t", quote = FALSE,
              row.names = TRUE, col.names = TRUE)
  
  
  dp <- QCPassed
  #print(dp)
  dp.1 <- data.frame(y = dp, x = 1:num.sample)
  gp <- ggplot(data = dp.1, aes(x = x, y = y))
  gp <- gp + geom_point(colour="#FF9999")
  gp <- gp + labs(title = "Sample coverage", x = "Samples", y = "Read Count/sample")
  #gp <- gp + scale_x_continuous(breaks = 1:num.sample, labels = sample.name)
  #gp <- gp + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #  scale_y_log10(breaks = c(1e3, 1e4, 1e5, 1e6, 1e7))
  gp <- gp + scale_y_log10(breaks = c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7), limits = c(1e2, 1e7))
  gp <- gp + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  gp <- gp + theme(axis.title = element_text(size = rel(2.2)))
  gp <- gp + theme(axis.text = element_text(size = rel(2.2)))
  gp <- gp + theme(plot.title = element_text(size = rel(2.5)))
  
  index <- which(dp < 1e4)
  if(length(index) > 0){
    text <- data.frame(x = index , y = dp[index] *0.75, text = sample.name[index])
    gp <- gp + geom_text(data = text, aes(x = x, y = y, label = text), size = 7)
  }
  
  wd <- num.sample * 0.2
  if(wd < 7)
    wd <- 7
  if(wd > 14)
    wd <- 14
  
  png(paste(out.dir, "SampleCoverage.png", sep = ""), width=wd, height=7, res=120, units = "in")
  #pdf(paste(path, "SampleCoverage.pdf", sep = ""), width = 14)
  print(gp)
  dev.off()
}


# uniformity
uniformity <- function(bam.dir, amplicon.bed, out.dir, sample.name, res.info){
  # bam.dir: bam file directory
  # amplicon.bed: bed file of amplicon
  # out.dir: output directory
  # sample.name: sample.name: name of sequenced samples
  # res.info: resource information
  
  message(" ")
  message("uniformity")
  bamfile <- list.files(path = bam.dir, pattern = "*.bam$")
  num.sample <- length(sample.name)
  
  bed.info <- read.table(amplicon.bed, sep = "\t", stringsAsFactors = FALSE)
  num.amp <- nrow(bed.info)
  read.count.all <- NULL
  rd.max <- rep(0, num.sample)
  rd.min <- rep(0, num.sample)
  rd.mean <- rep(0, num.sample)
  rd.mean.2 <- rep(0, num.sample)
  over.mean.2 <- rep("", num.sample)
  plot.om5 <- rep(0, num.sample)
  rd.mean.5 <- rep(0, num.sample)
  over.mean.5 <- rep("", num.sample)
  
  for(i in 1:num.sample){
    cmd <- paste0(res.info["bedtools"], " coverage -a ", amplicon.bed, " -b ", 
                  paste0(bam.dir, bamfile[i]), 
                  " -counts -F 0.9 | awk '{print $NF}'")
    message(cmd)
    read.count <- as.numeric(system(cmd, intern = TRUE))
    read.count.all <- cbind(read.count.all, read.count)
    rd.max[i] <- max(read.count)
    rd.min[i] <- min(read.count)
    rd.mean[i] <- mean(read.count)
    rd.mean.2[i] <- rd.mean[i] * 0.2
    rd.mean.5[i] <- rd.mean[i] * 0.5
    plot.om5[i] <- sum(read.count > rd.mean.5[i]) / num.amp
    over.mean.2[i] <- percent(sum(read.count > rd.mean.2[i]) / num.amp)
    over.mean.5[i] <- percent(sum(read.count > rd.mean.5[i]) / num.amp)
  }
  ratioMM <- rd.max /rd.min
  read.count.all <- rbind(read.count.all, ratioMM, rd.mean, rd.mean.2,
                          over.mean.2, rd.mean.5, over.mean.5)
  
  
  colnames(read.count.all) <- sample.name
  rownames(read.count.all) <- c(bed.info$V4, "max/min", "mean", "0.2X mean",
                                "%>0.2X mean", "0.5X mean", "%>0.5X mean")
  
  outfile <- paste0(out.dir, "uniformity.txt")
  write.table(read.count.all, file = outfile, sep = "\t", quote = FALSE,
              row.names = TRUE, col.names = TRUE)
  
  
  ap <- plot.om5 * 100
  #print(ap)
  #print(num.sample)
  sd.ap <- sd(ap)
  mean.ap <- mean(ap)
  ap.1 <- data.frame(y = ap, x = 1:length(ap))
  gp <- ggplot(data = ap.1, aes(x = x, y = y))
  gp <- gp + geom_point(colour="#FF9999")
  gp <- gp + labs(title = "Amplicon Uniformity", x = "Samples", y = "% amplicons > 0.5x mean")
  #gp <- gp + scale_x_continuous(breaks = 1:num.sample, labels = sample.name)
  #gp <- gp + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  #gp <- gp + geom_abline(slope = 0, intercept = mean.ap - sd.ap, colour = "#663333")
  gp <- gp + geom_abline(slope = 0, intercept = mean.ap - 2*sd.ap, colour = "#996666",
                         linetype = 2)
  gp <- gp + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  gp <- gp + theme(axis.title = element_text(size = rel(2.2)))
  gp <- gp + theme(axis.text = element_text(size = rel(2.2)))
  gp <- gp + theme(plot.title = element_text(size = rel(2.5)))
  
  index <- which(ap < mean.ap - 2*sd.ap)
  if(length(index) > 0){
    text <- data.frame(x = index , y = ap[index] - 1, text = sample.name[index])
    gp <- gp + geom_text(data = text, aes(x = x, y = y, label = text), size = 7)
  }
  
  text <- data.frame(x = length(ap) - 4, y = mean.ap - 2*sd.ap - 1.5, text = "2xSD")
  gp <- gp + geom_text(data = text, aes(x = x, y = y, label = text), size = 7)
  
  wd <- num.sample * 0.2
  if(wd < 7)
    wd <- 7
  if(wd > 14)
    wd <- 14
  
  png(paste(out.dir, "Uniformity.png", sep = ""), width=wd, height=7, res=120, units = "in")
  #pdf(paste(path, "Uniformity.pdf", sep = ""), width = 14)
  print(gp)
  dev.off()
}

