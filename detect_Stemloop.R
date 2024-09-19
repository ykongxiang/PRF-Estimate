detect_Stemloop <- function(frameEnd, mrnaSeq) {
    # 提取序列从移码结束位点到末尾的部分
    mrnaSeq <- mrnaSeq[frameEnd:length(mrnaSeq)]
    
    # 将序列转换为大写，以统一碱基表示
    mrnaSeq <- toupper(mrnaSeq)
    
    # 如果提取的序列长度小于4，由于茎环至少需要4个核苷酸，直接返回FALSE
    if (nchar(mrnaSeq) < 4) {
        return(FALSE) # 序列太短，无法形成茎环
    }
    
    # 定义碱基互补关系
    complement <- c(A = "U", U = "A", C = "G", G = "C")
    
    # 初始化最长配对长度变量
    maxPairLength <- 0
    
    # 遍历序列，寻找可能的配对起始点
    for (i in 1:(length(mrnaSeq) - 1)) {
        # 从可能的配对起始点开始，寻找配对终点
        for (j in (i + 2):(length(mrnaSeq))) {
            # 检查从i到j是否与从j到i+1完全互补
            if (all(mrnaSeq[i:(j - 1)] == complement[mrnaSeq[j:(i + 1)]])) {
                Stop_position <- detect_Stop_codons(mrnaSeq,i)#该函数返回逻辑值（未编写完善）
                #搜索停止密码子位置是否在寻找的茎环之前
                }
            }
        }
    return(Stop_positin)
}