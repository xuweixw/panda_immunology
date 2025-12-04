rm(list=ls())

setwd("/Users/apple/2022-日常工作/大熊猫繁育基地/2024-2025年/论文写作/免疫球蛋白抗体轻链/supplement/")

library("gggenomes")
library("dplyr")
library("tibble")

## 过滤重复的比对记录
filt <- function(large, small, ...) {
  for (i in 1:nrow(small)) {
    Sregion <- c(as.numeric(small[i,7]),as.numeric(small[i,8]))
    contain <- FALSE
    for (j in 1:nrow(large)) {
      Lregion <- c(as.numeric(large[j,7]),as.numeric(large[j,8]))
      if (Lregion[1] < Sregion[1] && Lregion[2]> Sregion[2]) {
        contain <- TRUE
        break
      }
    }
    if (!contain) {
      large <- bind_rows(large, small[i,])
    }
  }
  return(large)
}

seq <- tibble(seq_id=c("GPv2.IgL", "GPv1.IgL"), length=c(2579357, 2514939))

GPv1_results <- read_blast("data/GPv2_vs_GPv1.blastn") %>% arrange(-length)
GPlink100 <- GPv1_results %>% filter(length > 100000)
GPlink100_50 <- GPv1_results %>% filter(length < 100000, length >= 50000 )
GPlink50_25 <- GPv1_results %>% filter(length < 50000, length >= 25000)
GPlink25_10 <- GPv1_results %>% filter(length < 25000, length >= 10000, pident > 95)
GPlink10_5 <- GPv1_results %>% filter(length < 10000, length >= 5000, pident > 98)
GPlink <- filt(GPlink100, GPlink100_50)
GPlink <- filt(GPlink, GPlink50_25)
GPlink <- filt(GPlink, GPlink25_10)
GPlink <- filt(GPlink, GPlink10_5)


p <- gggenomes(seqs = seq, links = GPlink, spacing = 1 ) +
  geom_seq(linewidth = 0.5) +
  geom_seq_label(vjust = -3.5, hjust = 1.2, size = 3) +
  geom_link(offset = 0.1, fill="coral", color="white") +
  scale_x_continuous(
    limits=c(-0.1e+6, 2.6e+6),
    breaks = c(0, 0.5e+6, 1e+6, 1.5e+6, 2e+6, 2.5e+6),
    labels = c("0", "0.5", "1", "1.5", "2M", "2.5M")) +
  theme(legend.position = "NULL")
p

ggsave(filename="result/GPv2.IgL-vs-GPv1.IgL.2025-12-01.pdf", width=10, height=4)
