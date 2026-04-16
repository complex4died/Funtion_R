library(DESeq2)
library(tidyverse)


deseq_ana <- function(matrix, 
                      sample_info，
                     logFC_t){
  dds <- DESeqDataSetFromMatrix(
    countData = matrix,      
    colData = sample_info,                
    design = ~ Group                     
  )
  
  keep <- rowSums(counts(dds) >= 10) >= 3 
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  
  res <- data.frame(results(dds), 
                    stringsAsFactors = FALSE, 
                    check.names = FALSE)
  
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



