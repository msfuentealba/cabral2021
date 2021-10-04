library(tidyverse)
get_genes <- function(comb,d){
  #rank wt_fad
  subset <- comb %>% filter(dir==d)
  subset <- subset %>% arrange(desc(abs(aging)))
  subset$rank.x <- 1:nrow(subset)
  
  #rank fad_fx
  subset <- subset %>% arrange(desc(abs(xbp)))
  subset$rank.y <- 1:nrow(subset)
  
  n <- 100000
  set.seed(10)
  sim_prod <- replicate(n,prod(c(sample(1:nrow(subset),1),sample(1:nrow(subset),1))))
  subset$prod <- sapply(1:nrow(subset), function(x) prod(c(subset$rank.x[x],subset$rank.y[x])))
  subset <- subset %>% group_by(symbol) %>% summarise(prod = mean(prod)) %>% ungroup
  subset$pval <- sapply(subset$prod, function(x) sum(x>=sim_prod)/n)
  #subset <- subset %>% arrange(prod)
  #subset$fdr <- p.adjust(subset$pval, method = "BH")
  return(subset)
}

pairs_cols <- rep("black",3)

#comb_old <- read_rds("./output/data.rds")[,c(2,6)] %>% na.omit #aav-xbp1
comb_old <- read_rds("./output/data.rds")[,c(2,5)] %>% na.omit #xbp1
colnames(comb_old) <- c("aging","xbp")
comb_old$uniprot <- rownames(comb_old)
symbol2uniprot <- readxl::read_xlsx("./data/PROTEOMICS MAY 2020 FELIPE_combined-t-test-shared.xlsx", skip = 1)[,c(4,23)] %>% set_names(c("uniprot","symbol"))
symbol2uniprot$symbol <- gsub(";.*", "",symbol2uniprot$symbol)
comb_old <- comb_old %>% left_join(symbol2uniprot)
comb_old <- comb_old %>% group_by(symbol) %>% summarise(xbp = mean(xbp),aging=mean(aging))
comb_old$dir <- ifelse(comb_old$aging<0&comb_old$xbp>0,"dn_up", ifelse(comb_old$aging>0&comb_old$xbp<0,"up_dn","same"))
comb_old$dir <- factor(comb_old$dir, levels = c("dn_up","up_dn","same"))
dn_up_old <- get_genes(comb_old,"dn_up") %>% filter(pval<0.05) %>% select(symbol) %>% unlist %>% as.character()
up_dn_old <- get_genes(comb_old,"up_dn") %>% filter(pval<0.05) %>% select(symbol) %>% unlist %>% as.character()
hl_old <- comb_old %>% filter(symbol%in%c(dn_up_old,up_dn_old))
hl_old <- hl_old %>% left_join(tibble(dir=c("up_dn","dn_up"), col = brewer.pal(10, "Paired")[1:2]))

p_old <- ggplot(comb_old, aes(x = aging, y = xbp, color = dir)) +
  geom_point(size = 1.5, alpha = 0.1) +
  geom_point(data=hl_old, aes(x=aging,y=xbp), color=hl_old$col, size = 1.5)+
  geom_smooth(method = "lm", col = "red") +
  xlab("Aged vs Young")+
  ylab("TgXBP1s (Aged) vs\nNon Tg (Aged)")+
  scale_color_manual(values = pairs_cols[1:3])+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14))


comb_ma <- read_rds("./output/data.rds")[,c(1,4)] %>% na.omit
colnames(comb_ma) <- c("aging","xbp")
comb_ma$uniprot <- rownames(comb_ma)
symbol2uniprot <- readxl::read_xlsx("./data/PROTEOMICS MAY 2020 FELIPE_combined-t-test-shared.xlsx", skip = 1)[,c(4,23)] %>% set_names(c("uniprot","symbol"))
symbol2uniprot$symbol <- gsub(";.*", "",symbol2uniprot$symbol)
comb_ma <- comb_ma %>% left_join(symbol2uniprot)
comb_ma <- comb_ma %>% group_by(symbol) %>% summarise(xbp = mean(xbp),aging=mean(aging))
comb_ma$dir <- ifelse(comb_ma$aging<0&comb_ma$xbp>0,"dn_up", ifelse(comb_ma$aging>0&comb_ma$xbp<0,"up_dn","same"))
comb_ma$dir <- factor(comb_ma$dir, levels = c("dn_up","up_dn","same"))
dn_up_ma <- get_genes(comb_ma,"dn_up") %>% filter(pval<0.05) %>% select(symbol) %>% unlist %>% as.character()
up_dn_ma <- get_genes(comb_ma,"up_dn") %>% filter(pval<0.05) %>% select(symbol) %>% unlist %>% as.character()
hl_ma <- comb_ma %>% filter(symbol%in%c(dn_up_ma,up_dn_ma))
hl_ma <- hl_ma %>% left_join(tibble(dir=c("up_dn","dn_up"), col = brewer.pal(10, "Paired")[3:4]))

p_ma <- ggplot(comb_ma, aes(x = aging, y = xbp, color = dir)) +
  geom_point(size = 1.5, alpha = 0.1) +
  geom_point(data=hl_ma, aes(x=aging,y=xbp), color=hl_ma$col, size = 1.5)+
  geom_smooth(method = "lm", col = "red") +
  xlab("Middle Aged vs Young")+
  ylab("TgXBP1s (Middle Aged) vs\nNon Tg (Middle Aged)")+
  scale_color_manual(values = pairs_cols[1:3])+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14))

library(ggpubr)
#pdf(file = paste0("./output/scatter.pdf"), height = 7, width = 4, useDingbats = FALSE)
#ggarrange(p_ma,p_old, ncol = 1)
#dev.off()

# old <- rbind(get_genes(comb_old,"dn_up") %>% group_by(symbol) %>% summarise(score = (-log10(pval))*1),
#              get_genes(comb_old,"up_dn") %>% group_by(symbol) %>% summarise(score = (-log10(pval))*-1))
# old$Var2 <- "old"
# old <- old[,c(3,1,2)]
# colnames(old)[3] <- "value"

#library(GeneOverlap)
#library(gridExtra)

#testGeneOverlap(newGeneOverlap(dn_up_old,dn_up_ma, spec = c("mm9.gene")))
#intersect(dn_up_old,dn_up_ma)
#pdf(file = paste0("./output/venn_labels_up.pdf"), height = 4, width = 4)
#grid.table(as.data.frame(matrix(sort(intersect(dn_up_old,dn_up_ma)), ncol = 3)), theme = ttheme_minimal(), rows = NULL, cols = NULL)
#dev.off()
#testGeneOverlap(newGeneOverlap(up_dn_old,up_dn_ma, spec = c("mm9.gene")))
#intersect(up_dn_old,up_dn_ma)
#pdf(file = paste0("./output/venn_labels_dn.pdf"), height = 4, width = 4)
#grid.table(as.data.frame(matrix(sort(intersect(up_dn_old,up_dn_ma)), ncol = 3)), theme = ttheme_minimal(), rows = NULL, cols = NULL)
#dev.off()

#library(ggvenn)
#library(RColorBrewer)
#library(circlize)
# up <- ggvenn(list(`Middle Aged`  = dn_up_ma, Aged = dn_up_old),
#        show_percentage = FALSE,
#        stroke_color = c(rep("#0101EE",100),rep("black",100)),
#        fill_color = c("white","white"),
#        stroke_size = 2, set_name_size = 4)
# 
# dn <- ggvenn(list(`Middle Age`  = up_dn_ma, Aged = up_dn_old),
#              show_percentage = FALSE,
#              stroke_color = c(rep("#0101EE",100),rep("black",100)),
#              fill_color = c("white","white"),
#              stroke_size = 2, set_name_size = 4)

#library(ggpubr)
#pdf(file = paste0("./output/overlaps.pdf"), height = 5, width = 6)
#ggarrange(up,dn,ncol = 1, labels = c("Up-regulated","Down-regulated"), label.x = 0.25)
#dev.off()
#cormat <- read_rds("./data/cormat.rds")[,c(16,8,24)]
#cor(cormat, method = "s", use = "pairwise.complete.obs")

expr <- read_rds("./data/log2qn/Kang2011_HIP_age.rds")
age <-  read_rds("./data/ages/Kang2011_HIP_age.rds")
expr <- expr[,names(age)]
cors <- tibble(gene = rownames(expr), r = sapply(1:nrow(expr), function(x) cor(expr[x,],age, method = "s")))
head(cors)
orto <- read_tsv("./data/ortholog-results-table.tab")
orto <- orto %>% select(3,5) %>% set_names(c("gene","symbol")) %>% unique
cors <- cors %>% left_join(orto) %>% na.omit
cors <- cors %>% group_by(symbol) %>% summarise(r = mean(r))

#aav <- read_rds("./output/data.rds")
#aav$uniprot <- rownames(aav)
#aav <- aav[,c(7,6)] %>% left_join(symbol2uniprot) %>% left_join(cors)
#cor.test(aav$`AAV-XBP1s (Old)`,aav$r) # no global correlation with aging

all <- tibble(group = c("up_dn_old","up_dn_ma","dn_up_old","dn_up_ma"), genes = list(up_dn_old,up_dn_ma,dn_up_old,dn_up_ma))
all$r <- lapply(all$genes, function(x) cors$r[cors$symbol%in%x])

all2 <- all %>% select(group,r) %>% unnest
all2$age <- ifelse(grepl("old",all2$group),"Old","Young")
all2$age <- factor(all2$age, levels = c("Young","Old"))

all2$group <- factor(all2$group, levels = c("up_dn_old","dn_up_old","up_dn_ma","dn_up_ma"))
library(RColorBrewer)
p1 <- ggplot(all2[all2$age=="Young",], aes(x=group, y=r, color = group)) + 
  geom_boxplot(alpha = 1) +
  coord_flip() +
  ylab("Spearman's rho")+
  xlab("")+
  geom_hline(yintercept = 0, linetype="dotted", color = "black", size = 0.5) +
  geom_jitter(position=position_jitter(0.2)) +
  theme_bw() +
  ylim(-0.7,0.7)+
  scale_color_manual(values = brewer.pal(10, "Paired")[3:4]) +
  #coord_capped_cart(bottom='both', left='both') +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

p2 <- ggplot(all2[all2$age=="Old",], aes(x=group, y=r, color = group)) + 
  geom_boxplot(alpha = 1) +
  coord_flip() +
  ylab("Spearman's rho")+
  xlab("")+
  geom_hline(yintercept = 0, linetype="dotted", color = "black", size = 0.5) +
  geom_jitter(position=position_jitter(0.2)) +
  theme_bw() +
  ylim(-0.7,0.7)+
  scale_color_manual(values = brewer.pal(10, "Paired")[1:2]) +
  #coord_capped_cart(bottom='both', left='both') +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

pdf(file = paste0("./output/boxplot.pdf"), height = 5, width = 5, useDingbats = FALSE) 
ggarrange(p1,p2,ncol = 1)
#ggarrange(p2,ncol = 1)
dev.off()

t.test(all2$r[all2$group=="dn_up_ma"], alternative = "less")
t.test(all2$r[all2$group=="up_dn_ma"], alternative = "greater")

t.test(all2$r[all2$group=="dn_up_old"], alternative = "less")
t.test(all2$r[all2$group=="up_dn_old"], alternative = "greater")


cors_human <- cors %>% left_join(comb_old[,c(1:3)]) %>% data.frame %>% na.omit
rownames(cors_human) <- cors_human$symbol
cors_human$symbol <- NULL
cor(cors_human, method = "s", use = "pairwise.complete.obs")

library(GeneOverlap)
testGeneOverlap(newGeneOverlap(rownames(cors_human %>% arrange(r))[1:100],rownames(cors_human %>% arrange(aging))[1:100], spec = "mm9.gene"))
testGeneOverlap(newGeneOverlap(rownames(cors_human %>% arrange(desc(r)))[1:100],rownames(cors_human %>% arrange(desc(aging)))[1:100], spec = "mm9.gene"))


# cogn1 <- readxl::read_xlsx("~/Downloads/41467_2019_9613_MOESM12_ESM.xlsx", sheet = 2, col_names = FALSE) %>% set_names(c("symbol","dir"))
# cogn1$symbol <- gsub("\\|.*", "",cogn1$symbol)
# cogn2 <- readxl::read_xlsx("~/Downloads/41467_2019_9613_MOESM12_ESM.xlsx", sheet = 4, col_names = FALSE) %>% set_names(c("symbol","dir"))
# cogn2$symbol <- gsub("\\|.*", "",cogn2$symbol)
# #write.table(unique(c(cogn1$symbol,cogn2$symbol)), file = "./data/cogn_orto.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# cogn_orto <- read_tsv("./data/cogn_orto.tab")[,c(1,5)] %>% set_names(c("symbol","symbol2"))
# cogn1 <- cogn1 %>% left_join(cogn_orto) %>% na.omit
# cogn2 <- cogn2 %>% left_join(cogn_orto) %>% na.omit
# 
# #ma
# testGeneOverlap(newGeneOverlap(cogn1$symbol2[cogn1$dir=="up"], dn_up_ma, spec = "mm9.gene"))
# testGeneOverlap(newGeneOverlap(cogn1$symbol2[cogn1$dir=="down"], up_dn_ma, spec = "mm9.gene"))
# testGeneOverlap(newGeneOverlap(cogn2$symbol2[cogn2$dir=="up"], dn_up_ma, spec = "mm9.gene"))
# testGeneOverlap(newGeneOverlap(cogn2$symbol2[cogn2$dir=="down"], up_dn_ma, spec = "mm9.gene"))
# 
# #old
# testGeneOverlap(newGeneOverlap(cogn1$symbol2[cogn1$dir=="up"], dn_up_old, spec = "mm9.gene"))
# testGeneOverlap(newGeneOverlap(cogn1$symbol2[cogn1$dir=="down"], up_dn_old, spec = "mm9.gene"))
# testGeneOverlap(newGeneOverlap(cogn2$symbol2[cogn2$dir=="up"], dn_up_old, spec = "mm9.gene"))
# testGeneOverlap(newGeneOverlap(cogn2$symbol2[cogn2$dir=="down"], up_dn_old, spec = "mm9.gene"))

# gsets <- list(dn_up_ma=dn_up_ma,
#               up_dn_ma=up_dn_ma,
#               dn_up_old=dn_up_old,
#               up_dn_old=up_dn_old)
# library(fgsea)
# 
# aging <- readxl::read_xls("./data/xuetal.xls")
# aging$Accession <- gsub("-.*", "",aging$Accession)
# #write.table(aging$Accession, file = "./output/uniprot_aging.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
# orto_aging <- read_tsv("./data/aging_orto.tab")[,c(1,5)] %>% unique %>% set_names(c("uniprot","symbol"))
# aging <- aging %>% left_join(orto_aging, by = c("Accession"="uniprot")) %>% na.omit
# aging <- aging %>% group_by(symbol) %>% summarise(young = mean(log2(`Group B/A`)), ma = mean(log2(`Group C/A`)), old = mean(log2(`Group D/A`)))
# 
# aging %>% filter(symbol%in%dn_up_old) %>% summarise(young = median(young), ma = median(ma), old = median(old))
# aging %>% filter(symbol%in%up_dn_old) %>% summarise(young = median(young), ma = median(ma), old = median(old))
# aging %>% filter(symbol%in%dn_up_ma) %>% summarise(young = median(young), ma = median(ma), old = median(old))
# aging %>% filter(symbol%in%up_dn_ma) %>% summarise(young = median(young), ma = median(ma), old = median(old))
# 
# young <- aging$young
# names(young) <- aging$symbol
# fgsea(gsets,young)
# 
# ma <- aging$ma
# names(ma) <- aging$symbol
# fgsea(gsets,ma)
# 
# old <- aging$old
# names(old) <- aging$symbol
# fgsea(gsets,old)

# g_up <- all2 %>% filter(group=="Up-regulated")
# g_up$pathway <- factor(g_up$pathway, levels = c("dn_up_ma","dn_up_old") %>% rev)
# g_dn <- all2 %>% filter(group=="Down-regulated")
# g_dn$pathway <- factor(g_dn$pathway, levels = c("up_dn_ma","up_dn_old") %>% rev)
# 
# library(ggridges)
# 
# plot_dn <- ggplot(g_dn, aes(x = order, y = pathway, color = group, point_color = group, fill = group)) +
#   geom_density_ridges(jittered_points = TRUE, scale = .95, rel_min_height = .01,
#                       point_shape = "|", point_size = 3, size = 0.25,
#                       position = position_points_jitter(height = 0, width = 0.1)) +
#   scale_y_discrete(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0), name = "") +
#   scale_fill_manual(values = c("#D55E0050"), labels = c("Up-regulated"), guide = "none") +
#   scale_color_manual(values = c("#D55E00"), guide = "none") +
#   scale_discrete_manual("point_color", values = c("#D55E00"), guide = "none") +
#   coord_cartesian(clip = "off") +
#   theme_ridges(center = TRUE) +
#   theme(axis.title.y=element_blank(),
#         axis.text.x=element_blank())
# 
# plot_up <- ggplot(g_up, aes(x = order, y = pathway, color = group, point_color = group, fill = group)) +
#   geom_density_ridges(jittered_points = TRUE, scale = .95, rel_min_height = .01,
#                       point_shape = "|", point_size = 3, size = 0.25,
#                       position = position_points_jitter(height = 0, width = 0.1)) +
#   scale_y_discrete(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0), name = "") +
#   scale_fill_manual(values = c("#0072B250"), labels = c("Up-regulated"), guide = "none") +
#   scale_color_manual(values = c("#0072B2"), guide = "none") +
#   scale_discrete_manual("point_color", values = c("#0072B2"), guide = "none") +
#   coord_cartesian(clip = "off") +
#   theme_ridges(center = TRUE) +
#   theme(axis.title.y=element_blank(),
#         axis.text.x=element_blank())
# 
# 
# library(ggpubr)
# pdf(file = paste0("./output/gsea.pdf"), height = 5, width = 4) 
# ggarrange(plot_dn,plot_up,ncol = 1)
# dev.off()

# 
# 
# exp <- read_rds("./output/data.rds")
# exp$uniprot <- rownames(exp)
# symbol2uniprot <- readxl::read_xlsx("./data/PROTEOMICS MAY 2020 FELIPE_combined-t-test-shared.xlsx", skip = 1)[,c(4,23)] %>% set_names("uniprot","symbol")
# symbol2uniprot$symbol <- gsub(";.*", "",symbol2uniprot$symbol)
# exp <- exp %>% left_join(symbol2uniprot)
# exp <- reshape2::melt(exp)[,c(3,2,4)] %>% set_names(colnames(age))
# exp <- exp %>% group_by(Var2,symbol) %>% summarise(value = mean(value)) %>% ungroup
# age <- rbind(age,exp)
#jittered_points = TRUE, 
#point_shape = "|", point_size = 3, size = 0.25,
#position = position_points_jitter(height = 0),
#bandwidth = 0.1