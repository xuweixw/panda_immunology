setwd(dir = "/Users/apple/2022-日常工作/大熊猫繁育基地/2024-2025年/论文写作/免疫球蛋白抗体轻链/supplement/")


library("gggenes")
library("xlsx")
library("ggplot2")
library("ggrepel")
library("dplyr")

data <- read.xlsx("../Supplementary-A-2025-11-07.xlsx", sheetIndex = 1)
head(data)
dim(data)

data <- tibble(data)
width <- 3e5
data <- data %>% mutate(group = floor(V.EXON_End / width)) %>% 
  mutate(related_Start = L.part1_Start- group*width, related_End = V.EXON_End - group*width)

p <- ggplot(data, aes(xmin = related_Start, xmax = related_End, 
                 y = as.character(group), 
                 fill = Func,  color=Func)) +
  geom_gene_arrow(size=0.6, arrowhead_width = grid::unit(0, "mm"),
                  arrowhead_height = grid::unit(0, "mm")) +
  geom_text(aes(x=related_Start, label = Gene, hjust = -0.4), angle=30, size=2.5, color="black") +
  facet_wrap(~ group, scales = "free_y", ncol = 1)  +
  scale_x_continuous(name = "IGL genic positions", breaks=seq(0, 3e5, by=5e4)
                     ) +
  labs(y="") +
  theme_genes() +
  theme(legend.position = "NULL",
        axis.text.y = element_blank())

p


ggsave("result/IgL.gene.position-2025-12-01.pdf", height=8, width=8)
