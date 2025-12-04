rm(list=ls())

setwd("/Users/apple/2022-日常工作/大熊猫繁育基地/2024-2025年/论文写作/免疫球蛋白抗体轻链/supplement/")

library(ggplot2)
library(ggrepel)
library(dplyr)


BAC <- read.table("data/IGL.BAC.csv", header = T)
head(BAC)
dim(BAC)

# fill = "steelblue"
BAC <- BAC %>% filter(Level != 0)

p <- ggplot(BAC) +
  geom_segment(
    aes(x = GPv2Start, xend = GPv2End, y = Level, yend = Level, color = Haplotype, linewidth = 1),) +
  geom_text(
    aes(x = GPv2Start+0.5e5, y = Level-0.35, label = Name), size = 3) +
  scale_color_manual(values = c("A"="#00996c", "B"="coral")) +
  #scale_color_manual(values = c("N"="white", "Y"="darkgreen", size=2)) +
  lims(x=c(0, 2.7e6), y=c(-0.05, 5.2)) +
  scale_x_continuous(breaks = c(0, 0.5e6, 1e6, 1.5e6, 2e6, 2.5e6), 
                     labels = c("0", "0.5", "1", "1.5", "2", "2.5")) +
  labs(x = "BAC position in IgL locus (Mb)", y= NULL) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.line.x = element_line(color="black", linewidth =1.5),
        legend.position = "none")
p

ggsave("result/IGL.MTP.2025-12-01.pdf", p, width=10, height=2.5)
