rm(list = ls())

setwd("2022-日常工作/大熊猫繁育基地/2024-2025年/论文写作/免疫球蛋白抗体轻链/supplement/")

library("ggplot2")

BAC_small <- read.csv("data/BAC.small.size.csv", header = F)
head(BAC_small)
dim(BAC_small)

BAC_large <- read.csv("data/BAC.large.size.csv", header = F)
head(BAC_large)
dim(BAC_large)

p <- ggplot(BAC_large, aes(x=V1)) + 
  geom_histogram(binwidth = 2e3, color="white", alpha=0.5, fill="#4d7192") +
  geom_histogram(aes(V1), data=BAC_small, binwidth = 2e3, color="white", alpha=0.5, fill="darkgreen") +
  scale_x_continuous(limits = c(0.6e5, 2e5),
                     breaks=c(0.6e5, 0.8e5, 1e5, 1.2e5, 1.4e5, 1.6e5, 1.8e5, 2e5), 
                     labels = c("60", "80", "100", "120", "140", "160", "180", "200")) +
  labs(x="BAC insert sizes (Kb)", y="Number of BAC clones") +
  theme_classic()
p

ggsave("result/BAC-size-distribution.pdf", p, width=8, height = 6, dpi=300)
