library(DESeq2)
library(tidyverse)
library(BiocParallel)

deseq_ana <- function(matrix, sample_info, logFC_t, ncores = 4) {
  # 1. 设置并行后端
  register(MulticoreParam(ncores, progressbar = TRUE))
  
  # 2. 创建 DESeq2 对象
  dds <- DESeqDataSetFromMatrix(
    countData = matrix,      
    colData = sample_info,                
    design = ~ Group                     
  )
  
  # 3. 过滤低表达基因
  keep <- rowSums(counts(dds) >= 10) >= 3 
  dds <- dds[keep, ]
  
  # 4. 运行 DESeq2（开启并行）
  dds <- DESeq(dds, parallel = TRUE)
  
  # 5. 提取结果（开启并行）
  res <- data.frame(
    results(dds, parallel = TRUE), 
    stringsAsFactors = FALSE, 
    check.names = FALSE
  )
  
  # 6. 处理结果
  res_na <- res %>%
    filter(!is.na(pvalue)) %>% 
    arrange(desc(log2FoldChange)) %>% 
    mutate(
      change = case_when(
        pvalue > 0.05 ~ "Stable",
        abs(log2FoldChange) < logFC_t ~ "Stable",
        log2FoldChange >= logFC_t ~ "Up",
        log2FoldChange <= -logFC_t ~ "Down",
        TRUE ~ "Stable"
      )
    )
  
  return(res_na)
}
