rm(list=ls())
setwd("/Users/apple/2022-日常工作/大熊猫繁育基地/2024-2025年/论文写作/免疫球蛋白抗体轻链/supplement/")

library("ggseqlogo")
library("xlsx")
library("dplyr")
library("ggpubr")

dataL <- read.xlsx("../Supplementary-A-2025-11-07.xlsx", sheetIndex = 1)
head(dataL)

LVRS <- dataL %>% mutate(group = substr(Gene, 1, 4)) %>% filter(group =="IGLV") %>% select(HEPTAMER, NONAMER)
pLVheptamer <- ggseqlogo(LVRS$HEPTAMER)  +  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
pLVnonamer <- ggseqlogo(LVRS$NONAMER)  +  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

LJRS <- dataL %>% mutate(group = substr(Gene, 1, 4)) %>% filter(group =="IGLJ") %>% select(HEPTAMER, NONAMER)
pLJheptamer <- ggseqlogo(LJRS$NONAMER)  +  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
pLJnonamer <- ggseqlogo(LJRS$HEPTAMER)  +  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

## For IGK
dataK <- read.xlsx("../Supplementary-A-2025-11-07.xlsx", sheetIndex = 3)
head(dataK)

KVRS <- dataK %>% mutate(group = substr(Gene, 1, 4)) %>% filter(group =="IGKV") %>% select(HEPTAMER, NONAMER)
pKVheptamer <- ggseqlogo(KVRS$HEPTAMER) +  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
pKVnonamer <- ggseqlogo(KVRS$NONAMER)  +  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

KJRS <- dataK %>% mutate(group = substr(Gene, 1, 4)) %>% filter(group =="IGKJ") %>% select(HEPTAMER, NONAMER)
pKJheptamer <- ggseqlogo(KJRS$NONAMER)  +  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
pKJnonamer <- ggseqlogo(KJRS$HEPTAMER)  +  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

p <- ggarrange(pLVheptamer, pLVnonamer, pLJnonamer, pLJheptamer,
          pKVheptamer, pKVnonamer, pKJnonamer, pKJheptamer,
          ncol=4, nrow=2, 
          widths = c(7, 9, 9, 7))
p <- annotate_figure(p, left = text_grob(label = "IGL IGK", rot=90),
                bottom = text_grob(label = "V-HEPTAMER"))
p
ggsave(filename = "result/RS_Logo-2025-11-18.pdf", p, height = 2, width=8)
