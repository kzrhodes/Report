#loads indicated packages into RStudio
require(DESeq2)
require(RCurl)
require(biomaRt)
require(stringr)
#assigns variable to the data file read in
geneCounts<-read.delim("RNAseqCounts.txt",
                       head=T,sep="\t",skip=1, row.names=1)

#counts any expressions that are not zero

nonZeroCounts<-geneCounts[rowSums(geneCounts[,6:28])>0,6:28]

# samples are either plus or minus. 
treatments <- as.factor(str_sub(colnames(nonZeroCounts),-9,-5))
#starts DESeq analysis of matrix
dds <- DESeqDataSetFromMatrix(as.matrix(nonZeroCounts),
                              as.data.frame(treatments),
                              design=~treatments)

dds$treatments <- relevel(dds$treatments,"minus")
#plot results and write to output file
dds <- DESeq(dds)
plotMA(dds, main="DESeq2")
DEresults <- results(dds)
write.csv(as.data.frame(DEresults), file="treatments_treated_results.csv")

# get annotations
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                host="www.ensembl.org", 
                path="/biomart/martservice")

gg4 <- useDataset("ggallus_gene_ensembl",mart=mart)

annot <- getBM(attributes=c("ensembl_gene_id",
                            "external_gene_id",
                            "affy_chicken"),
               filter="ensembl_gene_id",
               values=row.names(DEresults),
               mart=gg4)

annotResults <- merge(annot,DEresults,by.x="ensembl_gene_id",by.y="row.names")
#write results to indicated output file
head(annotResults)
write.table(annotResults, "annotResults.txt", sep="\t", quote=FALSE)
return(annotResults)
#sort results and display the most differentially expressed genes
sortLimma<-head(annotResults[order(annotResults$padj),])
