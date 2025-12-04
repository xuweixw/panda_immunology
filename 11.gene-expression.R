rm(list = ls())
setwd(dir = "/Users/apple/2022-日常工作/大熊猫繁育基地/2024-2025年/论文写作/免疫球蛋白抗体轻链/supplement/")

library("stringr")
library("tidyr")
library("dplyr")
library("xlsx")
library("ggplot2")
library("ggpubr")

import <- function(ID) {
  data <- read.csv(file = paste0("data/expression/", ID, ".txt"), header = FALSE, sep="\t")
  data <- data %>%  mutate(New=trimws(V1))  %>% 
    separate(New, into = c(ID, "Name"), sep = " ") %>%
    select(all_of(c(ID, "Name"))) %>%
    mutate(across(ID, as.integer))
  return(data)
}

geneName <- read.csv("data/expression/geneList.txt",
                     header = FALSE, col.names = c("Name", "Len"), sep="\t")
head(geneName)

samples <-  read.table("data/expression/samples.txt", header=T)
head(samples)

for (i in samples$Sample){
  geneName <- left_join(geneName, import(i), by="Name")
}

geneName[is.na(geneName[,])] <- 0

write.table(geneName, file="expression.txt", quote = F)

geneName <- geneName %>% mutate(group = substr(Name, 1, 4)) 

## For IG lambda joining
dataIGLJ <- geneName %>% filter(group == "IGLJ") %>%
  mutate(across(`40`:jingjing, ~ round(.x /sum(.x), 2), .names = "{.col}_perc")) %>%
  select(Name, `40_perc`:jingjing_perc) %>%
  pivot_longer(cols=`40_perc`:jingjing_perc, names_to = "group", values_to = "perc" ) %>%
  group_by(Name) %>%
  summarise(Mean = mean(perc),
            SD =sd(perc))
pLJ <- ggplot(dataIGLJ, aes(x= reorder(Name, -Mean),, y=Mean, fill="lightcoral")) + 
  geom_col() +
  geom_errorbar(aes(ymin = Mean-0.01, ymax=Mean+SD), width=0.2, color = "lightcoral") +
  scale_x_discrete(label = c("IGLJ2", "IGLJ1", "IGLJ4", "IGLJ3")) +
  labs(x = "Lambda joining families", y = "Percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=30, vjust = 0.8, hjust = 0.8),
        legend.position = "NULL") 
pLJ

## For lambda variable
dataIGLV <- geneName %>% filter(group == "IGLV") %>%
  mutate(across(`40`:jingjing, ~ round(.x /sum(.x), 2), .names = "{.col}_perc")) %>%
  select(Name, `40_perc`:jingjing_perc) %>%
  pivot_longer(cols=`40_perc`:jingjing_perc, names_to = "group", values_to = "perc" ) %>%
  group_by(Name) %>%
  summarise(Mean = mean(perc),
            SD =sd(perc)) %>%
  arrange(-Mean) %>% head(n = 10)
pLV <- ggplot(dataIGLV, aes(x= reorder(Name, -Mean),, y=Mean, fill="lightcoral")) + 
  geom_col() +
  geom_errorbar(aes(ymin = Mean-0.01, ymax=Mean+SD), width=0.2, color = "lightcoral") +
  labs(x = "Top10 lambda variable genes", y = "Percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=30, vjust = 0.8, hjust = 0.8),
        legend.position = "NULL") 
pLV

## For IG kappa joining
dataIGKJ <- geneName %>% filter(group == "IGKJ") %>%
  mutate(across(`40`:jingjing, ~ round(.x /sum(.x), 2), .names = "{.col}_perc")) %>%
  select(Name, `40_perc`:jingjing_perc) %>%
  pivot_longer(cols=`40_perc`:jingjing_perc, names_to = "group", values_to = "perc" ) %>%
  group_by(Name) %>%
  summarise(Mean = mean(perc),
            SD =sd(perc))
pKJ <- ggplot(dataIGKJ, aes(x= reorder(Name, -Mean),, y=Mean, fill="lightcoral")) + 
  geom_col() +
  geom_errorbar(aes(ymin = Mean-0.01, ymax=Mean+SD), width=0.2, color = "lightcoral") +
  labs(x = "Kappa joining genes", y = "Percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=30, vjust = 0.8, hjust = 0.8),
        legend.position = "NULL") 
pKJ

## For kappa variable
dataIGKV <- geneName %>% filter(group == "IGKV") %>%
  mutate(across(`40`:jingjing, ~ round(.x /sum(.x), 2), .names = "{.col}_perc")) %>%
  select(Name, `40_perc`:jingjing_perc) %>%
  pivot_longer(cols=`40_perc`:jingjing_perc, names_to = "group", values_to = "perc" ) %>%
  group_by(Name) %>%
  summarise(Mean = mean(perc),
            SD =sd(perc)) %>%
  arrange(-Mean) %>% head(n = 10)
pKV <- ggplot(dataIGKV, aes(x= reorder(Name, -Mean),, y=Mean, fill="lightcoral")) + 
  geom_col() +
  geom_errorbar(aes(ymin = Mean-0.01, ymax=Mean+SD), width=0.2, color = "lightcoral") +
  labs(x = "Top10 kappa variable genes", y = "Percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=30, vjust = 0.8, hjust = 0.8),
        legend.position = "NULL") 
pKV

## For light chain constant genes
dataC <- geneName %>% filter(group %in% c("IGKC", "IGLC")) %>%
  mutate(across(`40`:jingjing, ~ round(.x /sum(.x), 2), .names = "{.col}_perc")) %>%
  select(Name, `40_perc`:jingjing_perc) %>%
  pivot_longer(cols=`40_perc`:jingjing_perc, names_to = "group", values_to = "perc" ) %>%
  group_by(Name) %>%
  summarise(Mean = mean(perc),
            SD =sd(perc))
pC <- ggplot(dataC, aes(x= reorder(Name, -Mean),, y=Mean, fill="lightcoral")) + 
  geom_col() +
  geom_errorbar(aes(ymin = Mean-0.01, ymax=Mean+SD), width=0.2, color = "lightcoral") +
  scale_x_discrete(label = c("IGLC", "IGKC")) +
  labs(x = "Light-chain constants", y = "Percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=30, vjust = 0.8, hjust = 0.8),
        legend.position = "NULL") 
pC

IGLC <- read.xlsx("data/expression/IGLC-301bp.xlsx", sheetIndex = 1)
head(IGLC)
IGLCmean <- IGLC %>% 
  group_by(Gene) %>%
  summarise(SRR9678840 = sum(SRR9678840), SRR9678841 = sum(SRR9678841), SRR9678842 = sum(SRR9678842),
            SRR9678843 = sum(SRR9678843), SRR9678844 = sum(SRR9678844), SRR9678845 = sum(SRR9678845),
            SRR9678846 = sum(SRR9678846), SRR9678847 = sum(SRR9678847), SRR9678848 = sum(SRR9678848),
            SRR9678849 = sum(SRR9678849), SRR9678850 = sum(SRR9678850), SRR9678851 = sum(SRR9678851),
            SRR9678852 = sum(SRR9678852), SRR9678853 = sum(SRR9678853), SRR9678854 = sum(SRR9678854),
            SRR9678855 = sum(SRR9678855), SRR9678856 = sum(SRR9678856), SRR9678857 = sum(SRR9678857),
            SRR9678858 = sum(SRR9678858), SRR9678859 = sum(SRR9678859), jingjing = sum(jingjing)) %>%
  pivot_longer(cols=-Gene, names_to = "Sample", values_to = "Read") %>%
  pivot_wider(names_from = Gene, values_from = Read) %>%
 mutate(total = IGLC1 + IGLC2, 
               IGLC1perc = IGLC1 / total,
               IGLC2perc = IGLC2 / total) %>%
  select(all_of(c("Sample", "IGLC1perc", "IGLC2perc"))) %>%
  pivot_longer(cols=c(IGLC1perc, IGLC2perc), names_to = "group", values_to = "perc" ) %>%
  group_by(group) %>%
  summarise(Mean = mean(perc),
            SD =sd(perc))
pLC <- ggplot(IGLCmean, aes(x= group, y=Mean, fill="lightcoral")) + 
  geom_col() +
  geom_errorbar(aes(ymin = Mean-0.01, ymax=Mean+SD), width=0.2, color = "lightcoral") +
  scale_x_discrete(label = c("IGLC1", "IGLC2")) +
  labs(x = "Lambda constant families", y = "Percentage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=30, vjust = 0.8, hjust = 0.8),
        legend.position = "NULL") 
pLC

p <- ggarrange(pC, pLC, pLJ,  pLV, pKJ,  pKV,
               labels = c("A", "B", "C", "D", "E", "F"),
               ncol=3, nrow=2)
p

row1 <- ggarrange(pC, pLC, pLJ, pKJ,
                  labels = c("A", "B", "C", "D"),
                  ncol = 4, nrow = 1)
row1
row2 <- ggarrange(pLV, pKV,
                  labels = c("E", "F"),
                  ncol = 2, nrow = 1)
row2
p2 <- ggarrange(row1, row2, 
          ncol = 1, nrow = 2,
          heights = c(1, 1))
p2

ggsave("result/11.expression-2025-11-17.pdf", plot = p2, height = 6, width = 9)



## 废弃代码
geneName$`40` <- round(geneName$`40`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="40"),2]/1e6), 2)
geneName$`41` <- round(geneName$`41`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="41"),2]/1e6), 2)
geneName$`42` <- round(geneName$`42`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="42"),2]/1e6), 2)
geneName$`43` <- round(geneName$`43`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="43"),2]/1e6), 2)
geneName$`44` <- round(geneName$`44`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="44"),2]/1e6), 2)
geneName$`45` <- round(geneName$`45`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="45"),2]/1e6), 2)
geneName$`46` <- round(geneName$`46`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="46"),2]/1e6), 2)
geneName$`47` <- round(geneName$`47`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="47"),2]/1e6), 2)
geneName$`48` <- round(geneName$`48`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="48"),2]/1e6), 2)
geneName$`49` <- round(geneName$`49`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="49"),2]/1e6), 2)
geneName$`50` <- round(geneName$`50`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="50"),2]/1e6), 2)
geneName$`51` <- round(geneName$`51`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="51"),2]/1e6), 2)
geneName$`52` <- round(geneName$`52`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="52"),2]/1e6), 2)
geneName$`53` <- round(geneName$`53`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="53"),2]/1e6), 2)
geneName$`54` <- round(geneName$`54`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="54"),2]/1e6), 2)
geneName$`55` <- round(geneName$`55`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="55"),2]/1e6), 2)
geneName$`56` <- round(geneName$`56`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="56"),2]/1e6), 2)
geneName$`57` <- round(geneName$`57`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="57"),2]/1e6), 2)
geneName$`58` <- round(geneName$`58`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="58"),2]/1e6), 2)
geneName$`59` <- round(geneName$`59`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="59"),2]/1e6), 2)
geneName$`jingjing` <- round(geneName$`jingjing`/ (geneName$Len/1e3) / (samples[which(samples$Sample=="jingjing"),2]/1e6), 2)
geneMean <- geneName %>% 
  rowwise() %>%
  mutate(AVG = mean(c_across(3:last_col())),
         STD = sd(c_across(3:last_col())),
         group = substr(Name, 1, 4)) %>% 
  select(all_of(c("Name", "AVG", "STD", "group"))) %>%
  arrange(-AVG)