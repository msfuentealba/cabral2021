library(tidyverse)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggpubr)

data <- readxl::read_xlsx("./data/proteomics_data.xlsx", skip = 2)[,c(23,75:79,81:83)] %>% data.frame
colnames(data) <- c("symbol","aavmock_old","aavxbp1_old","nontg_young","nontg_ma","nontg_old","tgxbp1_young","tgxbp1_ma","tgxbp1_old")
data$symbol <- gsub(";.*", "",data$symbol)
data <- data[!grepl("Krt",data$symbol),] # remove kerating
data <- data %>% group_by(symbol) %>% summarise_all("mean")
hist(as.matrix(data[,2:ncol(data)]))
data2 <- data.frame(aging_ma = data$nontg_ma/data$nontg_young,
                aging_old = data$nontg_old/data$nontg_young,
                xbp1_young = data$tgxbp1_young/data$nontg_young,
                xbp1_ma = data$tgxbp1_ma/data$nontg_ma,
                xbp1_old = data$tgxbp1_old/data$nontg_old)
rownames(data2) <- data$symbol
data <- log2(data2)
#write_rds(data, "./data/differential_abundance.rds")
#hist(as.matrix(data))
#apply(data, 2, function(x) length(unique(x)))

colnames(data) <- c("Middle-age","Old","XBP1s-Tg (Young)","XBP1s-Tg (Middle-age)","XBP1s-Tg (Old)")
cor <- cor(data, method = "s", use = "pairwise.complete.obs")
cor <- cor[1:5,1:5]
write_rds(data, "./data/differential_abundance.rds")

col_fun <- colorRamp2(seq(-1,1,2/10),rev(brewer.pal(11, "RdBu")))
perc <- function(v1,v2){
  sum(sign(v1)!=sign(v2), na.rm = TRUE)/sum(!is.na(v1)&!is.na(v2), na.rm = TRUE)
}

#perc(data$`Middle-age`,data$`XBP1s-Tg (Middle-age)`)
#perc(data$Old,data$`XBP1s-Tg (Old)`)

#loadfonts()
pdf(file = paste0("./output/correlations.pdf"), height = 5, width = 6) 
Heatmap(cor,
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_order = colnames(cor),
        column_order = colnames(cor),
        row_split = c("1 Aging", "1 Aging","2 XBP1s-Tg","2 XBP1s-Tg","2 XBP1s-Tg"), 
        column_split = c("1 Aging", "1 Aging","2 XBP1s-Tg","2 XBP1s-Tg","2 XBP1s-Tg"), 
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(round(cor[i,j],2), x = x, y = y, gp=gpar(fontsize=10, fontface = "bold" ))
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 1))
        },
        heatmap_legend_param = list(direction = "vertical", title = "Spearman's\nrho", title_position = "topcenter"),
        column_names_rot = 90,
        border = TRUE)
dev.off()

