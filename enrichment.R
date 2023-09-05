deg<-read.csv("ControlGroup-VS-TreatGroup.DESeq2.GeneDiffExpFilter.xls",sep="\t",check.names = F,header=T)
df<-as.data.frame(deg)

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("KEGG.db"))
suppressPackageStartupMessages(library("stringr"))
outdir<-"./"
vs<-"mmu"
repFC<-rep(1,length(deg[1]))
FC<-data.frame(deg[1],repFC)
colnames(FC)<-c("GeneID","log2FC")

gene.df=data.frame(FC$GeneID,FC$GeneID)

colnames(gene.df)=c("GeneID","ENTREZID")
gene.df1<-left_join(gene.df,FC)
FC.trans<-gene.df1[[3]]
names(FC.trans)<-gene.df1[[2]]

kegg_ORA<-enrichKEGG(gene         = gene.df$ENTREZID,
                     organism     = "mmu",
                     keyType = "kegg" ,
                     pvalueCutoff = 1 ,
                     qvalueCutoff = 1 ,
                     minGSSize = 1 ,
                     maxGSSize = 1000 ,
                     use_internal_data = T)
mylist<-str_split(kegg_ORA@result$Description," - Mus musculus")

kegg_ORA@result$Description<-unlist(lapply(mylist, function(x) x[1]))
kegg_ORA<-setReadable(kegg_ORA, OrgDb="org.Mm.eg.db",keyType ="ENTREZID")
#- Mus musculus
write.table(kegg_ORA@result,paste0(outdir,"/",vs,".Pathway.ORA.Enrich_all.xls"),sep = "\t",quote = F, row.names = F, col.names = T)
kegg_ORAsig <- kegg_ORA[kegg_ORA$qvalue<0.05, asis=T]
write.table(kegg_ORAsig@result,paste0(outdir,"/",vs,".Pathway.ORA.Enrich_filter.xls"),sep = "\t",quote = F, row.names = F, col.names = T)


BP_ORAgroup <- groupGO(gene     = as.character(gene.df$ENTREZID),
                 keyType  = "ENTREZID",
                 OrgDb    = "org.Mm.eg.db",
                 ont      = "BP",
                 level    = 2,
                 readable = TRUE)
write.table(BP_ORAgroup, file=paste0(vs, ".BP.ORA.group.xls"), sep="\t", quote=F, row.names=F, col.names=T)

BP_ORA <- enrichGO(gene       = gene.df$ENTREZID,
                OrgDb         = "org.Mm.eg.db",
                keyType       = "ENTREZID",
                minGSSize     = 1,
                maxGSSize     = 1000,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)
write.table(BP_ORA@result,paste0(outdir,"/",vs,".BP.ORA.Enrich_all.xls"),sep = "\t",quote = F, row.names = F, col.names = T)
BP_ORAsig <- BP_ORA[BP_ORA$qvalue<0.05, asis=T]
write.table(BP_ORAsig@result,paste0(outdir,"/",vs,".BP.ORA.Enrich_filter.xls"),sep = "\t",quote = F, row.names = F, col.names = T)

CC_ORAgroup <- groupGO(gene     = as.character(gene.df$ENTREZID),
                 keyType  = "ENTREZID",
                 OrgDb    = "org.Mm.eg.db",
                 ont      = "CC",
                 level    = 2,
                 readable = TRUE)
write.table(CC_ORAgroup, file=paste0(vs, ".CC.ORA.group.xls"), sep="\t", quote=F, row.names=F, col.names=T)

CC_ORA <- enrichGO(gene       = gene.df$ENTREZID,
                OrgDb         = "org.Mm.eg.db",
                keyType       = "ENTREZID",
                minGSSize     = 1,
                maxGSSize     = 1000,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)

write.table(CC_ORA@result,paste0(outdir,"/",vs,".CC.ORA.Enrich_all.xls"),sep = "\t",quote = F, row.names = F, col.names = T)
CC_ORAsig <- CC_ORA[CC_ORA$qvalue<0.05, asis=T]
write.table(CC_ORAsig@result,paste0(outdir,"/",vs,".CC.ORA.Enrich_filter.xls"),sep = "\t",quote = F, row.names = F, col.names = T)


MF_ORAgroup <- groupGO(gene     = as.character(gene.df$ENTREZID),
                 keyType  = "ENTREZID",
                 OrgDb    = "org.Mm.eg.db",
                 ont      = "MF",
                 level    = 2,
                 readable = TRUE)
write.table(MF_ORAgroup, file=paste0(vs, ".MF.ORA.group.xls"), sep="\t", quote=F, row.names=F, col.names=T)

MF_ORA <- enrichGO(gene       = gene.df$ENTREZID,
                OrgDb         = "org.Mm.eg.db",
                keyType       = "ENTREZID",
                minGSSize     = 1,
                maxGSSize     = 1000,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)

write.table(MF_ORA@result,paste0(outdir,"/",vs,".MF.ORA.Enrich_all.xls"),sep = "\t",quote = F, row.names = F, col.names = T)
MF_ORAsig <- MF_ORA[MF_ORA$qvalue<0.05, asis=T]
write.table(MF_ORAsig@result,paste0(outdir,"/",vs,".MF.ORA.Enrich_filter.xls"),sep = "\t",quote = F, row.names = F, col.names = T)


