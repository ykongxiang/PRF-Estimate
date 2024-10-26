library(openxlsx)
setwd("C:/Users/31598/Desktop/移码数据")
va_mRNA <- read.fasta('Euplotes_vannus_transcripts.fasta', seqtype = "DNA", as.string = TRUE)
va_mRNA_df <- data.frame(
  DNA_seqid = names(va_mRNA),
  Sequence = sapply(va_mRNA, function(x) paste(x, collapse = ""))
)
oct_mRNA <- read.fasta('Euplotes_octocarinatus_transcripts.fasta', seqtype = "DNA", as.string = TRUE)
oct_mRNA_df <- data.frame(
  DNA_seqid = names(oct_mRNA),
  Sequence = sapply(oct_mRNA, function(x) paste(x, collapse = ""))
)
all_mrna <- rbind(va_mRNA_df,oct_mRNA_df)
write.xlsx(all_mrna,'all_mrna.xlsx')
va_blastx <- read.table('blastx_Euplotes_vannus.tab', header = TRUE, sep = "\t")
oct_blastx <- read.table('blastx_Euplotes_octocarinatus.tab', header = TRUE, sep = "\t")
colnames(oct_blastx) <- c("DNA_seqid", "Pep_seqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
colnames(va_blastx) <- c("DNA_seqid", "Pep_seqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
all_blastx <- rbind(oct_blastx,va_blastx)
all_blastx$frameshift <- with(all_blastx, 
                                  mapply(function(dna, pep) {
                                    any(PRF$DNA_seqid == dna & PRF$Pep_seqid == pep)
                                  }, DNA_seqid, Pep_seqid))
# 假设 all_mrna 数据框中包含 DNA_seqid 和 Sequence 列
all_blastx <- merge(all_blastx, all_mrna[, c("DNA_seqid", "Sequence")], 
                    by = "DNA_seqid", all.x = TRUE)

write.xlsx(all_blastx,'all_blastx.xlsx')

va_PRF <- read.table('Euplotes_vannus_PRF.tab',header = T)
oct_PRF <- read.table('Euplotes_octocarinatus_PRF.tab',header = T)
PRF <- rbind(oct_PRF,va_PRF)
write.xlsx(PRF,'PRF.xlsx')

FScanR_results <- read_xlsx('FScanR_results.xlsx')
PRF
train_data <- 