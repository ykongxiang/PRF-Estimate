library(readxl)
library(tidyverse)
library(dplyr)
library(openxlsx)

oct_blastx <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/oct_blastx.xlsx")
colnames(oct_blastx) <- c("DNA_seqid", "Pep_seqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
oct_mrna <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/oct_mrna.xlsx")
oct_prf <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/oct_prf.xlsx")
oct_mrna$DNA_seqid <- gsub("Contig.*", "", oct_mrna$DNA_seqid)
va_blastx <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/va_blastx.xlsx")
colnames(va_blastx) <- c("DNA_seqid", "Pep_seqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
va_mrna <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/va_mrna.xlsx")
va_prf <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/va_prf.xlsx")
va_prf_mrna <- merge(va_prf, va_mrna, by = "DNA_seqid", all = TRUE)
oct_prf_mrna <- merge(oct_prf, oct_mrna, by ='DNA_seqid', all = TRUE)
#删除contig标注 使DNA信息纯净

na_indices <- which(is.na(oct_prf_mrna$Sequence))
for (i in seq_along(na_indices)) {
  na_index <- na_indices[i]
  for (j in seq(na_index, nrow(oct_prf_mrna))) {
    if (!is.na(oct_prf_mrna$Sequence[j])) {
      oct_prf_mrna$Sequence[na_index] <- oct_prf_mrna$Sequence[j]
      break 
    }
  }
}
na_rows <- which(is.na(oct_prf_mrna$Strand))
oct_prf_mrna <- oct_prf_mrna %>%
  group_by(Sequence) %>%
  mutate(HasNonNA = any(!is.na(Strand))) %>%
  ungroup()
oct_prf_mrna <- oct_prf_mrna %>%
  filter(!(is.na(Strand) & oct_prf_mrna$HasNonNA))
write.xlsx(va_prf_mrna,'/Volumes/Samsung_T5/移码突变/移码数据/va_prf_mrna.xlsx')
write.xlsx(oct_mrna,'/Volumes/Samsung_T5/移码突变/移码数据/oct_mrna.xlsx')



va_blastx <- read.table('H:/移码突变/移码数据/blastx_Euplotes_vannus.tab', header = TRUE, sep = "\t")
oct_blastx <- read.table('H:/移码突变/移码数据/blastx_Euplotes_octocarinatus.tab', header = TRUE, sep = "\t")
colnames(oct_blastx) <- c("DNA_seqid", "Pep_seqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
colnames(va_blastx) <- c("DNA_seqid", "Pep_seqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
all_blastx <- rbind(oct_blastx,va_blastx)
write.xlsx(all_blastx,'H:/移码突变/移码数据/all_blastx.xlsx')

#——————————————————————————————————————————————————————————————————————————————————————————
if (!require("Biostrings")) install.packages("Biostrings", dependencies = TRUE)
library(Biostrings)
#####调用示例：
a <- set_data_frame(extract_sequences('/Volumes/Samsung_T5/移码突变/移码数据/all_prf_data.xlsx'))
# 调用主函数
a <- extract_sequences('/Volumes/Samsung_T5/移码突变/移码数据/all_prf_data.xlsx')
mrna_file <- 'H:/移码突变/移码数据/all_prf_data.xlsx'
mrna_data1 <- read_excel(mrna_file)
mrna_data1 <- na.omit(mrna_data1)

#————————————————————————————————————————————————————————————————————————————————————————————————————
#  x 是 BLAST 输出的 DataFrame, `mrna_seq` 是包含 mRNA 序列的函数, `all_prf_data` 是包含 PRF 数据的 DataFrame
scored_results <- score_frameshift_mutation(x, all_mrna,all_prf_data)
view(scored_results)
all_mrna <- read_excel("H:/移码突变/移码数据/all_mrna.xlsx")
FScanR_results <- read_excel("H:/移码突变/移码数据/FScanR_results.xlsx")
all_blastx <- read_excel("H:/移码突变/移码数据/all_blastx.xlsx")
scored_results <- score_frameshift_mutation(FScanR_results, all_mrna,all_prf_data)
new_scored <-score_frameshift_mutation(all_prf, all_mrna,all_prf_data)
# 使用bind_rows()合并数据框
#------------------------------------------------------------------------------------------------------
# 使用示例
# 假设df是你的数据框，total_score是列名，阈值为threshold_low和threshold_high
evaluate_scores(scored_results, "total_score")
evaluate_scores(scored_results, "homology_score")
evaluate_scores(scored_results, "distance_score")
evaluate_scores(scored_results, "Alignment_Score")

na.omit(x)
all_prf <- na.omit(all_prf_data)
new_scored <-score_frameshift_mutation(all_prf, all_mrna,all_prf_data)

library(ggplot2)


write.xlsx(new_scored,'new_scored.xlsx')
write.xlsx(scored_results,'scored_results.xlsx')
evaluate_scores(new_scored, "new_distance")
evaluate_scores(scored_results, "new_distance")
