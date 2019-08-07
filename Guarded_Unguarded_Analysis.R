library(biomaRt)
library(ggplot2)
library(ChIPseeker)

load("~/Downloads/peakAnno_1D4D7D.rda")

genes <- (as.character(peakAnno_1D4D7D@anno@elementMetadata@listData$SYMBOL))
genes2 <- data.frame(table(genes))
genes2 <- genes2[order(genes2$Freq, decreasing = TRUE),]

#initializing a mart
ensembl <- useMart(host="www.ensembl.org","ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

a <- listAttributes(ensembl)
a[regexpr(pattern="gene",text=a$description,ignore.case=TRUE)==1,]

#Use human gene symbol 
filter <- c("hgnc_symbol") #TRY: listFilters(ensembl)
g  <- as.character(unique(genes2$genes))
attrib <- c("hgnc_symbol","start_position","end_position")
res = getBM(attributes=attrib,filters=filter,values=g,mart=ensembl)
res$length <- res$end_position - res$start_position

for(i in genes2$genes){
  genes2$length[genes2$genes == i] <- max(res$length[res$hgnc_symbol == i])
}

topdamage <- genes2
topdamage <- topdamage[is.finite(topdamage$length),]


model <- lm(Freq ~ length , topdamage)
topdamage$predicted <- predict(model,topdamage)
topdamage$pred_dif <- ifelse(topdamage$Freq/topdamage$predicted > 2, "Guarded","Normal")
topdamage[topdamage$Freq/topdamage$predicted < 0.5,"pred_dif"] <- "Unguarded"

ggplot(topdamage, aes(length, Freq, colour = pred_dif))+
  geom_point(alpha = 0.3)+
  #geom_point(data = topdamage[topdamage$pred_dif != "random chance",], aes(log(length), Freq, colour = pred_dif))+
  #geom_smooth(data= topdamage , aes(log(length),predicted), colour = "red")+
  theme_bw()+
  scale_colour_manual(values= c("darkred", "grey","skyblue"))+
  xlab("length")+
  ylab("Peak Count")+
  labs(title = "damaged genes (1d,4d,7d overlap)") + scale_x_log10()


hypo <- as.character(topdamage[topdamage$pred_dif == "Unguarded","genes"])
hyper <- as.character(topdamage[topdamage$pred_dif == "Guarded","genes"])

rownames(topdamage) <- as.character(topdamage$genes)

ggplot(head(topdamage[hypo,],n=20), aes(reorder(genes,-Freq), Freq))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("gene")+
  ylab("Peak count, 1d,4d,7d overlap")+
  labs(title = "Top 20 hypo-sensitive genes, corrected for length")

