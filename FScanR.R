if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("FScanR")
library(FScanR)
FScanR <- function(blastx_output, 
                   mismatch_cutoff  = 5,
                   evalue_cutoff    = 1e-5,
                   frameDist_cutoff = 10
) {
  
  blastx <- blastx_output
  colnames(blastx) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
  blastx <- blastx[complete.cases(blastx) & blastx$evalue <= evalue_cutoff & blastx$mismatch <= mismatch_cutoff,,drop=FALSE]
  
  blastx_freq <- table(blastx$qseqid)
  blastx_freq_cutoff <- blastx_freq[blastx_freq > 1]
  
  blastx_cutoff <- blastx[blastx$qseqid %in% names(blastx_freq_cutoff),,drop=FALSE]
  
  blastx_cutoff_sort <- blastx_cutoff[order(blastx_cutoff$qseqid, blastx_cutoff$sseqid, blastx_cutoff$qstart),,drop=FALSE]
  
  prf <- data.frame()
  if (nrow(blastx_cutoff_sort) > 0) {
    for (i in 2:nrow(blastx_cutoff_sort)) {
      mismatch <- blastx_cutoff_sort[i,5]
      qseqid <- blastx_cutoff_sort[i,1]
      sseqid <- blastx_cutoff_sort[i,2]
      qstart <- blastx_cutoff_sort[i,7]
      qend <- blastx_cutoff_sort[i,8]
      sstart <- blastx_cutoff_sort[i,9]
      send <- blastx_cutoff_sort[i,10]
      qframe <- blastx_cutoff_sort[i,13]
      qseqid_last <- blastx_cutoff_sort[i-1,1]
      sseqid_last <- blastx_cutoff_sort[i-1,2]
      qstart_last <- blastx_cutoff_sort[i-1,7]
      qend_last <- blastx_cutoff_sort[i-1,8]
      sstart_last <- blastx_cutoff_sort[i-1,9]
      send_last <- blastx_cutoff_sort[i-1,10]
      qframe_last <- blastx_cutoff_sort[i-1,13]
      strand <- ""
      
      if (qseqid == qseqid_last & sseqid == sseqid_last & qframe != qframe_last & qframe * qframe_last > 0) {
        if (qframe > 0 & qframe_last > 0) {
          frameStart <- qend_last
          frameEnd <- qstart
          pepStart <- send_last
          pepEnd <- sstart
          strand <- "+"
        } else if (qframe < 0 & qframe_last < 0) {
          frameStart <- qstart_last
          frameEnd <- qend
          pepStart <- send
          pepEnd <- sstart_last
          strand <- "-"
        }
        qDist <- frameEnd - frameStart - 1
        sDist <- pepEnd - pepStart
        FS_type <- qDist + (1 - sDist) * 3
        if (abs(qDist) <= frameDist_cutoff & abs(sDist) <= floor(frameDist_cutoff/3)) {
          prf_sub <- data.frame(as.character(qseqid), frameStart, frameEnd, as.character(sseqid), send_last + 1, sstart, FS_type, strand)
          prf <- rbind(prf, prf_sub)
        }
        prf <- prf[prf$FS_type < 3 & prf$FS_type > -3,,drop=FALSE]
      }
    }
    if (nrow(prf) > 0) {
      colnames(prf) <- c("DNA_seqid", "FS_start", "FS_end", "Pep_seqid", "Pep_FS_start", "Pep_FS_end", "FS_type", "Strand")
      prf$loci1 = paste(prf$DNA_seqid, prf$FS_start, sep="_")
      prf$loci2 = paste(prf$DNA_seqid, prf$FS_end, sep="_")
      prf$loci3 = paste(prf$Pep_seqid, prf$Pep_FS_start, sep="_")
      prf$loci4 = paste(prf$Pep_seqid, prf$Pep_FS_end, sep="_")
      prf = prf[!duplicated(prf$loci1),,drop=FALSE]
      prf = prf[!duplicated(prf$loci2),,drop=FALSE]
      prf = prf[!duplicated(prf$loci3),,drop=FALSE]
      prf = prf[!duplicated(prf$loci4),,drop=FALSE]
      prf = prf[,!colnames(prf) %in% c("loci1", "loci2","loci3", "loci4"),drop=FALSE]
    }
  } else {
    message("No PRF events detected!")
  }
  
  return(prf)
}
