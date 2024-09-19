##get all detected PRF period 
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
#------------------------------------------------------------------------------------------
#calculate all kinds of seqs'frequence and order it
set_data_frame <- function(df, column_name) {
  # 检查列名是否存在于数据框中
  if (!column_name %in% names(df)) {
    stop("Column name does not exist in the data frame.")
  }
  
  # 提取指定列
  x <- df[[column_name]]
  
  # 计算元素出现次数
  string_counts <- table(x)
  
  # 计算总数
  total_count <- sum(string_counts)
  
  # 计算频率
  string_frequencies <- string_counts / total_count
  
  # 创建数据框
  data <- data.frame(
    String = names(string_frequencies),
    Frequency = as.numeric(string_frequencies)
  )
  
  # 按频率降序排序
  data <- data[order(-data$Frequency), ]
  
  # 返回结果
  return(data)
}
#------------------------------------------------------------------------------------------------
#base on data to get first six strings and order it (data should at least contain cols named String and Frequence )
first_six <-function(data){
  first_six$String <- substr(data$String, start = 1, stop = 6)
  first_six <- first_six %>%
    group_by(String) %>%
    summarise(Sum = sum(Frequency))
  return(first_six)
}
#---------------------------------------------------------------------------------------------------------
#base on DNA_seqid to get matched mRNA sequences(when using the blastx data should conain DNA_seqid and the 
#sequence database named all_mrna)
mrna_seq <- function(blastx,all_mrna){
  blastx <- blastx %>%
    mutate(Sequence = all_mrna$Sequence[match(DNA_seqid,all_mrna$DNA_seqid)])
  return(blastx)
}
library(dplyr)

test <- function(data, prf, col_match = 'DNA_seqid', get = 'PRF') {
  get_col <- sym(get)
  match_col <- sym(col_match)
  data <- data %>%
    mutate(!!get := prf[[get]][match(data[[col_match]], prf[[col_match]])])
  
  return(data)
}

#-------------------------------------------------------------------------------------------------
mergeByMax <- function(a, b) {
  
  # 初始化一个空数据框用于存储最大分数的行
  max_score <- data.frame()
  
  # 遍历数据框的每一行
  for (i in 1:nrow(a)) {
    # 比较两个数据框中相同行的 Alignment_Score 值
    if (a$Alignment_Score[i] >= b$Alignment_Score[i]) {
      # 如果 a 中的分数较大或相等，添加 a 中的行
      max_score <- rbind(max_score, a[i, , drop = FALSE])
    } else {
      # 否则，添加 b 中的行
      max_score <- rbind(max_score, b[i, , drop = FALSE])
    }
  }
  
  # 返回合并后的数据框
  return(max_score)
}

#-----------------------------------------------------------------------------------------------------------------
library(Biostrings)
library(dplyr)

score_frameshift_mutation <- function(blast_df, mrna_sequences,all_prf_data) {
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
  
  blast_df$new_distance <- sapply(blast_df$distance, function(x) {
    match_value <- prf_distance$Frequency[which(prf_distance$String == x)]
    if (length(match_value) == 0) {
      return(NA) # 如果没有找到匹配项，返回NA或其他默认值
    } else {
      return(match_value)
    }
  })
  
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
  
  score_seq <- function(df, standard = 'AAATAA', match_score = 1, mismatch_score = -2, gap_open = -2, gap_extend = -2) {
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

#-------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(pROC)
library(coin)
evaluate_scores <- function(df, total_score_col, low=0.30, high=0.80) {
  
  # 提取 PRF 和非 PRF 分数
  PRF_score <- df[[total_score_col]][df$PRF == FALSE]
  nPRF_score <- df[[total_score_col]][df$PRF == TRUE]
  
  threshold_low <- quantile(PRF_score, low)
  threshold_high <- quantile(PRF_score, high)
  # t-test
  t_test_result <- t.test(PRF_score, nPRF_score)
  print(t_test_result)
  
  # Mann-Whitney U 检验
  wilcox_test_result <- wilcox.test(PRF_score, nPRF_score)
  print(wilcox_test_result)
  
  #  AUC 值
  labels <- ifelse(df$PRF == FALSE, 1, 0)
  roc_obj <- roc(labels, df[[total_score_col]])
  plot(roc_obj)
  auc_value <- auc(roc_obj)
  print(auc_value)
  # 分类数据
  df$Category_Pred <- cut(
    df[[total_score_col]],
    breaks = c(-Inf, threshold_low, threshold_high, Inf),
    labels = c("Low Probability A", "Possible A", "High Probability A")
  )
  
  # 计算混淆矩阵
  confusion_matrix <- table(Predicted = df$Category_Pred, Actual = df$PRF)
  print(confusion_matrix)
  
  # 计算准确率
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  print(paste("Accuracy:", accuracy))
  
  # 计算精确度
  precision <- confusion_matrix["High Probability A", "TRUE"] / 
    (confusion_matrix["High Probability A", "TRUE"] + confusion_matrix["High Probability A", "FALSE"])
  print(paste("Precision:", precision))
  
  # 计算假阳性率
  false_positive_rate <- confusion_matrix["High Probability A", "FALSE"] / 
    (confusion_matrix["High Probability A", "FALSE"] + confusion_matrix["Low Probability A", "FALSE"])
  print(paste("False Positive Rate:", false_positive_rate))
  
  # 计算假阴性率
  false_negative_rate <- confusion_matrix["Low Probability A", "TRUE"] / 
    (confusion_matrix["Low Probability A", "TRUE"] + confusion_matrix["High Probability A", "TRUE"])
  print(paste("False Negative Rate:", false_negative_rate))
  
  # 计算召回率
  recall <- confusion_matrix["High Probability A", "TRUE"] / 
    (confusion_matrix["High Probability A", "TRUE"] + confusion_matrix["Low Probability A", "TRUE"])
  print(paste("Recall:", recall))
  
  #df$PRF <- as.factor(df$PRF)
  #Permutation Test
  #perm_test_result <- independence_test(PRF ~ total_score, data = df)
  #print(paste('Permutation Test:',perm_test_result))
  
  #Spearman Rank Correlation:
  spearman_correlation <- cor(df[[total_score_col]], as.numeric(df$PRF), method = "spearman")
  print(paste('Spearman Rank Correlation:',spearman_correlation))
}
