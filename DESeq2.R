library(DESeq2)
library(dplyr)
mycounts <- read.csv("readscount.xls",sep="\t",check.names = F,header=T)
mycounts<-as.data.frame(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
colData<- read_file("group.xls",fileheader=TRUE)
colData<-as.data.frame(colData)
colData<-colData[,1:2]
rownames(colData)<-colData[,1]
colnames(colData)<-c("Samples","Group")
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ Group)
dds <- DESeq(dds)
sizefactor <- DESeq2::sizeFactors(dds)
write.table(sizefactor, file="sizefacor.xls", quote=FALSE, sep="\t")
normcounts<-counts(dds, normalized=TRUE)
res <- results(dds, contrast = c('Group', "TreatGroup", "ControlGroup"))
Control_name<-rownames(colData%>%dplyr::filter(Group=="ControlGroup"))
Control_values<- as.data.frame(normcounts) %>% dplyr::select(all_of(Control_name)) %>% mutate_if(is.numeric, round, digits = 1)
Treat_name<-rownames(colData%>%dplyr::filter(Group=="TreatGroup"))
Treat_values<- as.data.frame(normcounts) %>% dplyr::select(all_of(Treat_name)) %>% mutate_if(is.numeric, round, digits = 1)
Treat_values1 <- apply(Treat_values, 1, function(row) paste(row, collapse = ";"))
Control_values1 <- apply(Control_values, 1, function(row) paste(row, collapse = ";"))
AB_df <- data.frame(Control_NormalizedReadsCount = Control_values1, Treat_NormalizedReadsCount = Treat_values1)
mydf<-cbind(rownames(res),AB_df,res)
gene_diff_exp<-mydf%>%
dplyr::filter(log2FoldChange>0 | log2FoldChange<0)%>%
dplyr::mutate(up_down=ifelse(abs(log2FoldChange)>=1 & padj<=0.05,ifelse(log2FoldChange>0, 'Up', 'Down'), '*'))%>%
dplyr::mutate(up_down=ifelse(is.na(up_down), "*", up_down))%>%
dplyr::select(`rownames(res)`,Control_NormalizedReadsCount,Treat_NormalizedReadsCount,log2FoldChange, pvalue, padj, up_down)
gene_diff_filter<-filter(gene_diff_exp, up_down!='*')
names(gene_diff_exp)<-c("GeneID", "ControlGroup_NormalizedReadsCount", "TreatGroup_NormalizedReadsCount", "log2FoldChange(TreatGroup/ControlGroup)", "Pvalue", "Padj", "Up/Down-Regulation(TreatGroup/ControlGroup)")
names(gene_diff_filter)<-names(gene_diff_exp)
ExpFull<-sprintf("ControlGroup-VS-TreatGroup.DESeq2.GeneDiffExp.xls")
write.table(gene_diff_exp, file=ExpFull, quote=F, sep="\t", row.names=FALSE)
ExpFilter<-sprintf("ControlGroup-VS-TreatGroup.DESeq2.GeneDiffExpFilter.xls")
write.table(gene_diff_filter, file=ExpFilter, quote=F, sep="\t", row.names=FALSE)

