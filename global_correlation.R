library(tidyverse)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

data <- readxl::read_xlsx("./data/PROTEOMICS MAY 2020 FELIPE_combined-t-test-shared.xlsx", skip = 1)[,c(4,30:44)] %>% data.frame
rownames(data) <- data$Accession
data$Accession <- NULL
data <- data[,c(4:5,7,10,11,3)]
colSums(data==100, na.rm = TRUE)
colSums(data==0.01, na.rm = TRUE)
data[data==100] <- NA
data[data==0.01] <- NA
data <- log2(data)

colnames(data) <- c("Middle-age","Old","XBP1s-Tg (Young)","XBP1s-Tg (Middle-age)","XBP1s-Tg (Old)","AAV-XBP1s (Old)")
#write_rds(data, file = "./output/data.rds")

plot(data$`XBP1s-Tg (Old)`,data$`AAV-XBP1s (Old)`)
sum((data$`XBP1s-Tg (Old)`>0&data$`AAV-XBP1s (Old)`>0)|(data$`XBP1s-Tg (Old)`<0&data$`AAV-XBP1s (Old)`<0), na.rm = TRUE)/(sum(!is.na(data$`XBP1s-Tg (Old)`)&!is.na(data$`AAV-XBP1s (Old)`)))
plot(data$Old,data$`AAV-XBP1s (Old)`)
sum((data$Old>0&data$`AAV-XBP1s (Old)`<0)|(data$Old<0&data$`AAV-XBP1s (Old)`>0), na.rm = TRUE)/(sum(!is.na(data$Old)&!is.na(data$`AAV-XBP1s (Old)`)))

g1 <- ggplot(data, aes(y = `XBP1s-Tg (Old)`, x = `AAV-XBP1s (Old)`)) +
        geom_point(size = 1.5, alpha = 0.1) +
        #geom_point(data=hl_ma, aes(x=aging,y=xbp), color=hl_ma$col, size = 1.5)+
        geom_smooth(method = "lm", col = "red") +
        ylab("TgXBP1s (Aged) vs\nNon Tg (Aged)")+
        xlab("AAV-XBP1s (Aged) vs\nNon Tg (Aged)")+
        #scale_color_manual(values = pairs_cols[1:3])+
        theme_bw()+
        theme(legend.position = "none",
              axis.title = element_text(size = 16), 
              axis.text = element_text(size = 14))
g2 <- ggplot(data, aes(y = Old, x = `AAV-XBP1s (Old)`)) +
        geom_point(size = 1.5, alpha = 0.1) +
        #geom_point(data=hl_ma, aes(x=aging,y=xbp), color=hl_ma$col, size = 1.5)+
        geom_smooth(method = "lm", col = "red") +
        ylab("Aged vs Young")+
        xlab("AAV-XBP1s (Aged) vs\nNon Tg (Aged)")+
        #scale_color_manual(values = pairs_cols[1:3])+
        theme_bw()+
        theme(legend.position = "none",
              axis.title = element_text(size = 16), 
              axis.text = element_text(size = 14))

library(ggpubr)
pdf(file = paste0("./output/aav.pdf"), height = 5, width = 9) 
ggarrange(g1,g2,nrow = 1)
dev.off()

cor <- cor(data, method = "s", use = "pairwise.complete.obs")
cor

col_fun <- colorRamp2(seq(-1,1,2/10),rev(brewer.pal(11, "RdBu")))
perc <- function(v1,v2){
        sum(sign(v1)!=sign(v2), na.rm = TRUE)/sum(!is.na(v1)&!is.na(v2), na.rm = TRUE)
}
perc(data$`(non-Tg, 12) / (non-Tg, 3)`,data$`(XBP1s-Tg, 12) / (non-Tg, 12)`)
perc(data$`(non-Tg, 18) / (non-Tg, 3)`,data$`(XBP1s-Tg, 12) / (non-Tg, 12)`)
perc(data$`(non-Tg, 12) / (non-Tg, 3)`,data$`(XBP1s-Tg, 18) / (non-Tg, 18)`)
perc(data$`(non-Tg, 18) / (non-Tg, 3)`,data$`(XBP1s-Tg, 18) / (non-Tg, 18)`)
perc(data$`(non-Tg, 12) / (non-Tg, 3)`,data$`(AAV-XBP1s, 18) / (non-Tg, 18)`)
perc(data$`(non-Tg, 18) / (non-Tg, 3)`,data$`(AAV-XBP1s, 18) / (non-Tg, 18)`)

pdf(file = paste0("./output/correlations.pdf"), height = 5, width = 6) 

Heatmap(cor,
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_order = colnames(cor),
        column_order = colnames(cor),
        row_split = c("1 Aging", "1 Aging","2 XBP1s-Tg","2 XBP1s-Tg","2 XBP1s-Tg","3 AAV-XBP1s"), 
        column_split = c("1 Aging", "1 Aging","2 XBP1s-Tg","2 XBP1s-Tg","2 XBP1s-Tg","3 AAV-XBP1s"), 
        #cluster_row_slices = FALSE,
        #cluster_column_slices = FALSE,
        
        #column_split = rep(c("C", "D"), 3)
        
        rect_gp = gpar(col = "white", lwd = 1),
        #rect_gp = gpar(type = "none")
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(round(cor[i,j],2), x = x, y = y, gp=gpar(fontsize=10, fontface = "bold" ))
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 1))
        },
        heatmap_legend_param = list(direction = "vertical", title = "Spearman's\nrho", title_position = "topcenter"),
        column_names_rot = 45,
        #top_annotation = ano_col,
        #left_annotation = ano_row,
        #column_dend_side = "bottom",
        #row_dend_side = "left",
        #show_row_dend = FALSE,
        #show_column_dend = FALSE,
        #show_row_names=FALSE,
        #show_column_names=FALSE,
        border = TRUE)
dev.off()


# age <- readxl::read_xlsx("../natalia/data/datasets.xlsx", sheet = 4)[,c(1,6:8)] %>% set_names("uniprot","b","c","d")
# age$uniprot <- gsub("-.*", "",age$uniprot)
# mouse2human <- read_tsv("../natalia/data/mouse2human.txt") %>% set_names("mouse","human") %>% na.omit
# uniprot2genes <- read_tsv("../natalia/data/uniprot2genes.txt") %>% set_names("uniprot","human") %>% na.omit
# age <- age %>% left_join(uniprot2genes) %>% left_join(mouse2human)
# age_62 <- age[,c(6,2)] %>% set_names("symbol","logfc") %>% na.omit
# age_62$group <- "age_62"
# age_84 <- age[,c(6,3)] %>% set_names("symbol","logfc") %>% na.omit
# age_84$group <- "age_84"
# age_95 <- age[,c(6,4)] %>% set_names("symbol","logfc") %>% na.omit
# age_95$group <- "age_95"
# 
# symbol2uniprot <- readxl::read_xlsx("./data/PROTEOMICS MAY 2020 FELIPE_combined-t-test-shared.xlsx", skip = 1)[,c(4,23)]
# symbol2uniprot$`Gene Symbol` <- gsub(";.*", "",symbol2uniprot$`Gene Symbol`)
# data2 <- data[,3:6] %>% set_names(c("young","ma","old","old2")) 
# data2$Accession <- rownames(data2)
# data2 <- reshape2::melt(data2)
# data2 <- data2 %>% left_join(symbol2uniprot)
# head(data2)
# data2 <- data2[,c(4,3,2)] %>% set_names("symbol","logfc","group")
# 
# all <- rbind(age_62,age_84,age_95,data2)
# mat <- reshape2::acast(all, symbol~group, value.var = "logfc", fun.aggregate = mean) %>% na.omit
# mat[is.nan(mat)] <- NA
# cor(mat, method = "s", use = "pairwise.complete.obs")
