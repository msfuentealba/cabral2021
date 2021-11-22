library(tidyverse)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggpubr)

get_genes <- function(comb,d){
  subset <- comb %>% filter(dir==d)
  subset <- subset %>% arrange(desc(abs(aging)))
  subset$rank.x <- 1:nrow(subset)
  subset <- subset %>% arrange(desc(abs(xbp)))
  subset$rank.y <- 1:nrow(subset)
  subset$prod <- sapply(1:nrow(subset), function(x) prod(c(subset$rank.x[x],subset$rank.y[x])))
  subset <- subset %>% group_by(symbol) %>% summarise(prod = mean(prod)) %>% ungroup
  subset <- subset %>% arrange(prod)
  return(subset$symbol[1:50])
}

pairs_cols <- rep("black",3)
comb_old <- read_rds("./data/differential_abundance.rds")[,c(2,5)] %>% na.omit #xbp1 old
colnames(comb_old) <- c("aging","xbp")
comb_old$symbol <- rownames(comb_old)
comb_old <- comb_old %>% group_by(symbol) %>% summarise(xbp = mean(xbp),aging=mean(aging))
comb_old$dir <- ifelse(comb_old$aging<0&comb_old$xbp>0,"dn_up", ifelse(comb_old$aging>0&comb_old$xbp<0,"up_dn","same"))
comb_old$dir <- factor(comb_old$dir, levels = c("dn_up","up_dn","same"))
dn_up_old <- get_genes(comb_old,"dn_up") 
up_dn_old <- get_genes(comb_old,"up_dn") 
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


comb_ma <- read_rds("./data/differential_abundance.rds")[,c(1,4)] %>% na.omit #xbp1 middle age
colnames(comb_ma) <- c("aging","xbp")
comb_ma$symbol <- rownames(comb_ma)
comb_ma <- comb_ma %>% group_by(symbol) %>% summarise(xbp = mean(xbp),aging=mean(aging))
comb_ma$dir <- ifelse(comb_ma$aging<0&comb_ma$xbp>0,"dn_up", ifelse(comb_ma$aging>0&comb_ma$xbp<0,"up_dn","same"))
comb_ma$dir <- factor(comb_ma$dir, levels = c("dn_up","up_dn","same"))
dn_up_ma <- get_genes(comb_ma,"dn_up") 
up_dn_ma <- get_genes(comb_ma,"up_dn") 
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

pdf(file = paste0("./output/scatter.pdf"), height = 7, width = 4, useDingbats = FALSE)
ggarrange(p_ma,p_old, ncol = 1)
dev.off()

## kang
expr <- read_rds("./data/human_data/log2qn/Kang2011_HIP_age.rds")
age <-  read_rds("./data/human_data/ages/Kang2011_HIP_age.rds")
expr <- expr[,names(age)]
cors <- tibble(gene = rownames(expr), r = sapply(1:nrow(expr), function(x) cor(expr[x,],age, method = "s")))
head(cors)
orto <- read_tsv("./data/human_mouse_ortologues.tab")
orto <- orto %>% select(3,5) %>% set_names(c("gene","symbol")) %>% unique
cors <- cors %>% left_join(orto) %>% na.omit
cors <- cors %>% group_by(symbol) %>% summarise(r = mean(r))

all <- tibble(group = c("up_dn_old","up_dn_ma","dn_up_old","dn_up_ma"), genes = list(up_dn_old,up_dn_ma,dn_up_old,dn_up_ma))
all$r <- lapply(all$genes, function(x) cors$r[cors$symbol%in%x])
all2 <- all %>% select(group,r) %>% unnest
all2$age <- ifelse(grepl("old",all2$group),"Old","Young")
all2$age <- factor(all2$age, levels = c("Young","Old"))
all2$group <- factor(all2$group, levels = c("up_dn_old","dn_up_old","up_dn_ma","dn_up_ma"))

p1 <- ggplot(all2[all2$age=="Young",], aes(x=group, y=r, color = group)) + 
  coord_flip() +
  ylab("Spearman's rho")+
  xlab("")+
  geom_hline(yintercept = 0, linetype="dotted", color = "black", size = 0.5) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5) +
  geom_jitter(position=position_jitter(0.25), shape = 16, size = 2) +
  theme_bw() +
  ylim(-0.7,0.7)+
  scale_color_manual(values = brewer.pal(10, "Paired")[3:4]) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

p2 <- ggplot(all2[all2$age=="Old",], aes(x=group, y=r, color = group)) + 
  coord_flip() +
  ylab("Spearman's rho")+
  xlab("")+
  geom_hline(yintercept = 0, linetype="dotted", color = "black", size = 0.5) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5) +
  geom_jitter(position=position_jitter(0.25), shape = 16, size = 2) +
  theme_bw() +
  ylim(-0.7,0.7)+
  scale_color_manual(values = brewer.pal(10, "Paired")[1:2]) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

pdf(file = paste0("./output/boxplot_kang.pdf"), height = 5, width = 5, useDingbats = FALSE) 
ggarrange(p1,p2,ncol = 1)
dev.off()

t.test(all2$r[all2$group=="dn_up_ma"], alternative = "less")$p.value
t.test(all2$r[all2$group=="up_dn_ma"], alternative = "greater")$p.value
t.test(all2$r[all2$group=="dn_up_old"], alternative = "less")$p.value
t.test(all2$r[all2$group=="up_dn_old"], alternative = "greater")$p.value


