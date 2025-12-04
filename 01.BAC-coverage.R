rm(list=ls())

setwd("2022-日常工作/大熊猫繁育基地/2024-2025年/论文写作/免疫球蛋白抗体轻链/supplement/")

library("xlsx")
library("ggplot2")
library("dplyr")
library("tidyr")

data <- read.xlsx("data/BAC-summary.xlsx", sheetIndex = 1)
data$Chrom <- factor(data$Chrom, levels=data$Chrom)

head(data)
tail(data)

new_data <- data %>% pivot_longer(c(Coverage_Len, Full_Len), names_to = "type", values_to = "Len" )
head(new_data)

p <- ggplot(new_data, aes(Chrom, Len, fill=type)) + 
  geom_bar(stat="identity", position="dodge", width=0.6) +
  scale_y_continuous(breaks=c(2.5e7, 5e7, 7.5e7, 1e8, 1.25e8, 1.5e8, 1.75e8, 2e8),
                     labels = c("25", "50", "75", "100", "125", "150", "175", "200"), 
                     expand = c(0, 0)) +
  scale_fill_discrete(labels=c("BAC Coverage", "Chromosome Length")) +
  labs(x="Chromosome", y="Length (Mb)", fill="") + 
  theme_classic() +
  theme(legend.position = "top") 
p

ggsave("result/BAC-coverage.pdf", p, width=8, height=6)
