setwd(dir = "2022-日常工作/大熊猫繁育基地/2024-2025年/论文写作/免疫球蛋白抗体轻链/supplement/")


library("gggenes")
library("xlsx")
library("ggplot2")
library("dplyr")

data <- read.xlsx("../Supplementary-A-2025-11-07.xlsx", sheetIndex = 3)
head(data)
dim(data)


p <- ggplot(data, aes(xmin = L.part1_Start, xmax = V.EXON_End, y="IGK",
                 fill = Func, label = Gene, color=Func)) +
  geom_gene_arrow(size=0.4, arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(0, "mm")) +
  geom_gene_label() +
  lims(x=c(1.7e5, 4.2e5)) +
  # geom_gene_label(angle=45) +
  #scale_x_continuous(name = "IGK genic positions", breaks=seq(0, 4e5, by=5e4),
  #                   labels = function(x) paste0(x/5e4, " Kb")) +
  labs(y="") +
  theme_genes() +
  theme(legend.position = "NULL")

p


ggsave("result/IgK.gene.position-2025-11-17.pdf", height=1, width=6)
