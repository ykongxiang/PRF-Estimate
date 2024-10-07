library(Biostrings)
library(dplyr)
library(ggplot2)
library(pROC)
library(coin)
library(openxlsx)
#base on DNA_seqid to get matched mRNA sequences(when using the blastx data should conain DNA_seqid and the 
#sequence database named all_mrna)
mrna_seq <- function(blastx,all_mrna){
  blastx <- blastx %>%
    mutate(Sequence = all_mrna$Sequence[match(DNA_seqid,all_mrna$DNA_seqid)])
  return(blastx)
}
#find a sequence best match AAATAA/AAATAG-------------------------------------------------------------------------------------------------
mergeByMax <- function(a, b) {
  max_score <- data.frame()
  for (i in 1:nrow(a)) {
    if (a$Alignment_Score[i] >= b$Alignment_Score[i]) {
      max_score <- rbind(max_score, a[i, , drop = FALSE])
    } else {
      max_score <- rbind(max_score, b[i, , drop = FALSE])
    }
  }
  return(max_score)
}
#a function used to get needed data from A -------------------------------------------------------------------------------
test <- function(data, prf, col_match = 'DNA_seqid', get = 'PRF') {
  get_col <- sym(get)
  match_col <- sym(col_match)
  data <- data %>%
    mutate(!!get := prf[[get]][match(data[[col_match]], prf[[col_match]])])
  
  return(data)
}
#dispose negative strand-------------------------------------------------------------------------------
extract_sequences <- function(mrna_data) {
  mrna_data <- na.omit(mrna_data)  
  mrna_data$FS_start <- as.numeric(mrna_data$FS_start)
  mrna_data$FS_end <- as.numeric(mrna_data$FS_end)
  
  extract_and_reverse_substring <- function(text, start, end, strand) {
    if (strand == '-') {
      text <- as.character(reverseComplement(DNAString(text)))
    }
    # 提取子串
    substring <- substr(text, start - 2, end + 1)
    if (start > end){
      substring <- substr(text,end - 2, start + 1)
    }
    return(substring)
  }
  
  seqs <- sapply(1:nrow(mrna_data), function(i) {
    extract_and_reverse_substring(
      text = mrna_data$Sequence[i],
      start = mrna_data$FS_start[i],
      end = mrna_data$FS_end[i],
      strand = mrna_data$Strand[i]
    )
  })
  return(seqs)
}
#总打分代码-------------------------------------------------------------------------------
score_frameshift_mutation <- function(blast_df, mrna_sequences, all_prf_data, all_blastx) {
  blast_df <- test(blast_df,all_blastx,col_match = 'DNA_seqid' ,get = 'pident')
  blast_df <- test(blast_df,all_blastx,col_match = 'DNA_seqid' ,get = 'bitscore')
  blast_df <- test(blast_df,all_blastx,col_match = 'DNA_seqid' ,get = 'evalue')
  blast_df<- mrna_seq(blast_df,mrna_sequences)
  blast_df$FS_Period <- extract_sequences(blast_df)
  # 1. 计算两段序列的同源性得分
  blast_df$homology_score <- with(blast_df, {
    evalue_score <- -log10(evalue + 1e-10)
    #homology_score <- (pident / 100) * bitscore * evalue_score
    return(evalue_score)
  })
  
  # 2. 计算间隔距离得分
  blast_df$distance_score <- with(blast_df, {
    mean_distance <- 8.5  # 假设最有可能发生移码的平均距离
    sigma <- 1.5  # 标准差
    max_score <- dnorm(mean_distance, mean = mean_distance, sd = sigma)
    
    distance <- abs(FS_end - FS_start) + 1
    distance_score <- ifelse(distance >= 7 & distance <= 10, max_score, -log10(dnorm(distance, mean = mean_distance, sd = sigma) + 1e-10))
    return(distance_score)
  })
  blast_df$distance <- with(blast_df, {
    distance <- abs(FS_end - FS_start) + 1
    return(distance)
  })
  
  #blast_df$new_distance <- sapply(blast_df$distance, function(x) {
    #match_value <- prf_distance$Frequency[which(prf_distance$String == x)]
    #if (length(match_value) == 0) {
     # return(NA) # 如果没有找到匹配项，返回NA或其他默认值
    #} else {
    #  return(match_value)
    #}
  #})
  
  # 3. 计算序列模式匹配得分，同时保留alignment信息
  align_seq <- function(seq1, standard = 'AAATAA', match_score = 1, mismatch_score = -2, gap_open = -2, gap_extend = -2) {
    nucleotide_subs_matrix <- matrix(
      c(match_score, mismatch_score, mismatch_score, mismatch_score,
        mismatch_score, match_score, mismatch_score, mismatch_score,
        mismatch_score, mismatch_score, match_score, mismatch_score,
        mismatch_score, mismatch_score, mismatch_score, match_score),
      nrow = 4, dimnames = list(c("A", "C", "G", "T"), c("A", "C", "G", "T"))
    )
    
    alignment <- pairwiseAlignment(
      seq1, standard, 
      type = "local", 
      substitutionMatrix = nucleotide_subs_matrix, 
      gapOpening = gap_open, 
      gapExtension = gap_extend
    )
    
    result <- list(
      score = score(alignment),
      pattern = as.character(aligned(pattern(alignment))),
      subject = as.character(aligned(subject(alignment))),
      start_pattern = start(pattern(alignment)),
      end_pattern = end(pattern(alignment)),
      start_subject = start(subject(alignment)),
      end_subject = end(subject(alignment))
    )
    return(result)
  }
  
  score_seq <- function(df, standard = 'AAATAA', match_score = 1, mismatch_score = -2, gap_open = -2, gap_extend = -2){
    alignment_results <- sapply(df$FS_Period, align_seq, 
                                standard = standard, 
                                match_score = match_score, 
                                mismatch_score = mismatch_score, 
                                gap_open = gap_open, 
                                gap_extend = gap_extend, simplify = FALSE)
    
    df$Alignment_Score <- sapply(alignment_results, function(x) x$score)
    df$Aligned_Pattern <- sapply(alignment_results, function(x) x$pattern)
    df$Aligned_Subject <- sapply(alignment_results, function(x) x$subject)
    df$Pattern_Start <- sapply(alignment_results, function(x) x$start_pattern)
    df$Pattern_End <- sapply(alignment_results, function(x) x$end_pattern)
    df$Subject_Start <- sapply(alignment_results, function(x) x$start_subject)
    df$Subject_End <- sapply(alignment_results, function(x) x$end_subject)
    
    return(df)
  }
  
  df1 <- score_seq(blast_df, 'AAATAA', 1, -2, -2, -2)
  df2 <- score_seq(blast_df, 'AAATAG', 1, -2, -2, -2)
  
  df <- mergeByMax(df1,df2)
  df <- test(df, all_prf_data)
  df$total_score <- with(df, {
    total_score <- homology_score * distance_score * (1 + Alignment_Score)
    return(total_score)
  })
  return(df)
}
#评估打分机制-------------------------------------------------------------------------------
evaluate_scores <- function(df, total_score_col) {
  
  # 提取 PRF 和非 PRF 分数
  PRF_score <- df[[total_score_col]][df$PRF == FALSE]
  nPRF_score <- df[[total_score_col]][df$PRF == TRUE]
  
  # t-test
  t_test_result <- t.test(PRF_score, nPRF_score)
  print(t_test_result)
  
  # Mann-Whitney U 检验
  wilcox_test_result <- wilcox.test(PRF_score, nPRF_score)
  print(wilcox_test_result)
  
  # AUC 值
  labels <- ifelse(df$PRF == FALSE, 1, 0)
  roc_obj <- roc(labels, df[[total_score_col]])
  plot(roc_obj)
  auc_value <- auc(roc_obj)
  print(auc_value)
  
  # Spearman Rank Correlation
  spearman_correlation <- cor(df[[total_score_col]], as.numeric(df$PRF), method = "spearman")
  print(paste('Spearman Rank Correlation:', spearman_correlation))
  
  # Pearson Correlation
  pearson_correlation <- cor(df[[total_score_col]], as.numeric(df$PRF), method = "pearson")
  print(paste('Pearson Correlation:', pearson_correlation))
}
#using example:
path <- "C:/Users/31598/Desktop/移码数据"
setwd(path)
all_mrna <- read.xlsx("all_mrna.xlsx")
FScanR_results <- read.xlsx("FScanR_results.xlsx")
all_blastx <- read.xlsx("all_blastx.xlsx")
all_prf_data <- read.xlsx("all_prf_data.xlsx")
scored_results <- score_frameshift_mutation(FScanR_results, all_mrna,all_prf_data,all_blastx)
# total_score是列名，阈值为threshold_low和threshold_high(表示可能发生移码的概率，默认为0.3和0.8)
evaluate_scores(scored_results, "total_score")
evaluate_scores(scored_results, "homology_score")
evaluate_scores(scored_results, "distance_score")
evaluate_scores(scored_results, "Alignment_Score")
combined_df <- read.xlsx('C:/Users/31598/Documents/combined_df.xlsx')