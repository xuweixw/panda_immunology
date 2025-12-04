rm(list=ls())

setwd("2022-日常工作/大熊猫繁育基地/2024-2025年/论文写作/免疫球蛋白抗体轻链/supplement/")

# 加载必要的R包
library(ggplot2)
library(dplyr)
library(slider)
library(ggpubr)

# The number of paired-end reads for whole genome sequencing data
SampleDepth <- c("362"=199496532, "386"=147786839, "598"=165923654, "831"=197165638, "963"=178996617)
factorDepth <- mean(SampleDepth) / SampleDepth
factorDepth

depthDealed <- function(file, ID) {
  depth <- read.table(file, header=FALSE, 
                         col.names = c("Chromosome", "Position", "Depth"))
  depth_smoothed <- depth %>%
    mutate(Mean=slide_dbl(Depth, mean, .before = 5000, .after = 5000, .complete = TRUE)) %>%    # window-size=5001
    filter(!is.na(Mean)) %>%
    filter(row_number() %% 500 ==0) %>%      # step = 500
    mutate(NormalizedMean = Mean*factorDepth[ID])
  return(depth_smoothed)
}

depth362 <- depthDealed("data/362.depth.tab", "362")
depth386 <- depthDealed("data/386.depth.tab", "386")
depth598 <- depthDealed("data/598.depth.tab", "598")
depth831 <- depthDealed("data/831.depth.tab", "831")
depth963 <- depthDealed("data/963.depth.tab", "963")


depthPlot <- function(data, ID) {
  p <- ggplot(filter(data, Position >9.17e7, Position< 9.19e7), aes(x = Position, y =  NormalizedMean)) +
    geom_line(color = "darkgreen", linewidth = 0.5) +
    #geom_hline(yintercept = mean(data$NormalizedMean), linetype="dashed", color="red") +
    labs(x= NULL, y = ID) +
    lims(x= c(9.17e7, 9.19e7), y=c(8, 32)) +
    theme(panel.background = NULL,
          axis.text.x = element_blank())
  return(p)
}
p362 <- depthPlot(depth362, "362")
p362
p386 <- depthPlot(depth386, "386")
p598 <- depthPlot(depth598, "598")
p831 <- depthPlot(depth831, "831")


ggplot(depth362, aes(x = Position, y =  NormalizedMean)) +
  geom_line(color = "darkgreen", size = 0.5) +
  labs(x= NULL, y = "963") +
  lims(y=c(8, 32)) +
  theme(panel.background = NULL)
p963

#  绘制平滑后的深度图
p963 <- ggplot(filter(depth963, Position >9.17e7, Position< 9.19e7), 
               aes(x = Position, y =  NormalizedMean)) +
  geom_line(color = "darkgreen", size = 0.5) +
  labs(x= NULL, y = "963") +
  lims(y=c(8, 32)) +
  scale_x_continuous(breaks=c(9.17e7, 9.175e7, 9.18e7, 9.185e7, 9.19e7),
                   labels = c("91.7", "91.75", "91.8", "91.85", "91.9")) +
  theme(panel.background = NULL)
p963


combined <- ggarrange(p362, p386, p598, p831, p963, ncol = 1, nrow = 5)
# 保存图片
p <- annotate_figure(combined, bottom = text_grob("Genome position (MB)", size=10))
ggsave("result/depth-region-91M800K-2025-10-23.pdf", p, width = 8, height = 6)



## For small region
depthDealed2 <- function(file, ID) {
  depth <- read.table(file, header=FALSE, 
                      col.names = c("Chromosome", "Position", "Depth"))
  depth_smoothed <- depth %>%
    mutate(Mean=slide_dbl(Depth, mean, .before = 1000, .after = 1000, .complete = TRUE)) %>%    # window-size=2001
    filter(!is.na(Mean)) %>%
    filter(row_number() %% 200 ==0) %>%      # step = 200
    mutate(NormalizedMean = Mean*factorDepth[ID])
  return(depth_smoothed)
}
depth362Small <- depthDealed2("data/362.depth.tab", "362")
depth386Small <- depthDealed2("data/386.depth.tab", "386")
depth598Small <- depthDealed2("data/598.depth.tab", "598")
depth831Small <- depthDealed2("data/831.depth.tab", "831")
depth963Small <- depthDealed2("data/963.depth.tab", "963")

depthPlot2 <- function(data, ID) {
  ggplot(filter(data, Position >9.164e7, Position< 9.168e7), aes(x = Position, y =  NormalizedMean)) +
    geom_line(color = "darkgreen", size = 0.5) +
    labs(x= NULL, y = ID) +
    lims( y=c(0, 80)) +
    theme(panel.background = NULL,
          axis.text.x = element_blank())
}

p362Small <- depthPlot2(depth362Small, "362")
p362Small
p386Small <- depthPlot2(depth386Small, "386")
p598Small <- depthPlot2(depth598Small, "598")
p831Small <- depthPlot2(depth831Small, "831")

p963Small <-   ggplot(filter(depth963Small, Position >9.164e7, Position< 9.168e7), aes(x = Position, y =  NormalizedMean)) + 
  geom_line(color = "darkgreen", size = 0.5) +
  labs(x= NULL, y = "963") +
  lims(y=c(0, 80)) +
  scale_x_continuous(breaks=c(91.64e6, 91.65e6, 91.66e6, 91.67e6, 91.68e6),
                     labels = c("91.64", "91.65", "91.66", "91.67", "91.68")) +
  theme(panel.background = NULL)
p963Small

combined2 <- ggarrange(p362Small, p386Small, p598Small, p831Small, p963Small, ncol = 1, nrow = 5)
p2 <- annotate_figure(combined2, bottom = text_grob("Genome position (MB)", size=10))
p2
ggsave("result/depth-region-91M660K-2025-10-23.pdf", p2, width = 8, height = 6)
