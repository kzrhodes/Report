# Example of MicroArray Analysis taken from 
# adapted from
#		analysing-microarray-data-in-bioconductor
# on
#		http://bioinformatics.knowledgeblog.org
#
# download the BioC installation routines
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("GEOquery")
#biocLite("hgu133plus2.db")
library(ggplot2)
library(reshape2)
library(GEOquery)
library(affy)
library(affyPLM)
library(arrayQualityMetrics)
library(annotate)
library(chicken.db)
library(chickenprobe)
library(genefilter)
library(limma)

.ggboxIt <- function(dfin,pdfname)
{  # draw a ggplot boxplot of data frame produced by .summarise.it()
  gp <- ggplot(dfin, aes(x = name, ymin = lwisk, lower = lbox,
                         middle = mid, upper = ubox,
                         ymax = uwisk, fill = name))
  gp <- gp + geom_boxplot(stat="identity",show_guide=F) + coord_flip() 
  gp <- gp + theme(panel.background = element_blank())
  gp <- gp + labs(y = "Log2(Intensity)", x = "")
  ggsave(pdfname) 
}

.summariseIt <- function(dfraw,phenoData)
{ # create quantile data from a dataframe of depths
  # useful for munging into ggplot, however probably a plyr function for this too
  dfin <- melt(dfraw)
  colnames(dfin) <- c("probeId","sample","value")
  dfin$group<-factor(dfin$sample)
  stats <- as.data.frame(t(sapply(levels(dfin$sample),
                                  function(x){
                                    quantile(log2(as.numeric(dfin$value))[which(dfin$group==x)],
                                             prob=c(0.05,0.25,0.5,0.75,0.95),na.rm=T)
                                  })))
  colnames(stats)<-c("lwisk","lbox","mid","ubox","uwisk")
  stats$name <- rownames(stats)
  rownames(stats) <- NULL
  return(stats)
}

workingDir = "~/BS32010/assignments/Report_data/GCOS_Cel"
.getData <- function()
{
	filenames <- c("ROS1-_9.CEL", "ROS1+_10.CEL", "ROS2-_11.CEL",
							"ROS2+_12.CEL", "ROS3-_13.CEL", "ROS3+_14.CEL",
							"ROS4-_15.CEL", "ROS4+_16.CEL")
	samplenames <- c("ROS1-_9", "ROS1+_10", "ROS2-_11",
								"ROS2+_12", "ROS3-_13", "ROS3+_14",
								"ROS4-_15", "ROS4+_16")
	targets <- c("chicken", "chicken", "chicken", "chicken", "chicken", "chicken",
	             "chicken", "chicken")

	phenodata<-as.data.frame(cbind(filenames,samplenames,targets))
	write.table(phenodata,paste(workingDir,"phenodata.txt",sep="/")
							,quote=F,row.name=F)
	celRAW <- ReadAffy(celfile.path=workingDir,compress=T,
										 phenoData=phenodata)
	#arrayQualityMetrics(expressionset=celRAW,
	#										outdir=paste0(baseDir,"celRAW_AQM"),
	#										force=T,do.logtransform=T)
}

.plotDensity <- function(exps,filename)
{
	pdf(filename)
	# Plot a density vs log intensity histogram for the unnormalised data
	d <- apply(exps,2,
				function(x){
					density(x)
				})
	xmax <- max(sapply(d,function(x)max(x$x)))
	xmin <- min(sapply(d,function(x)min(x$x)))
	ymax <- max(sapply(d,function(x)max(x$y)))
	plot(0,pch='',ylab='',xlab='',
			 xlim=c(xmin,round(xmax+1)),ylim=c(0,ymax))
	lapply(1:length(d),function(x) lines(d[[x]],col=x))
	dev.off()
}

# Perform probe-level metric calculations on the CEL files:
.doPLM <- function(celRAW)
{
	pdf("celRAWqc.pdf")
	#affyPLM is required to interrogate celRMA
	celRAWqc <- fitPLM(celRAW)

	# Create an image of GSM24662.CEL:
	image(celRAWqc, which=1, add.legend=TRUE)

	# Create an image of GSM524665.CEL
	# There is a spatial artifact present
	image(celRAWqc, which=4, add.legend=TRUE)

	# affyPLM also provides more informative boxplots
	# RLE (Relative Log Expression) plots should have
	# values close to zero. GSM524665.CEL is an outlier
	RLE(celRAWqc, main="RLE")

	# We can also use NUSE (Normalised Unscaled Standard Errors).
	# The median standard error should be 1 for most genes.
	# GSM524665.CEL appears to be an outlier on this plot too
	NUSE(celRAWqc, main="NUSE")
	dev.off()
}

.doCluster <- function(celRMA)
{
#Quick and dirty clustering much better use pvclust
#Clustering like normalization is a bit of a black art
	eset <- exprs(celRMA)
	distance <- dist(t(eset),method="maximum")
	clusters <- hclust(distance)
	plot(clusters)
}

.doFilter <- function(celRMA)
{
	celfiles.filtered <- nsFilter(celRMA, 
															require.entrez=FALSE, 
															remove.dupEntrez=FALSE)
}


.doDE <- function(eset)
{
	samples <- eset$targets
# check the results of this
# convert into factors
	samples <- as.factor(samples)
# set up the experimental design

	design <- model.matrix(~0 + samples)
	colnames(design) <- c("choroid", "huvec", "iris", "retina")

# fit the linear model to the filtered expression set
	fit <- lmFit(exprs(eset), design)

# set up a contrast matrix to compare tissues v cell line
	contrast.matrix <- makeContrasts(huvec_choroid = huvec - choroid,
																 huvec_retina = huvec - retina, 
																 huvec_iris <- huvec - iris, 
																 levels=design)

# check the contrast matrix
	contrast.matrix

# Now the contrast matrix is combined with the per-probeset linear model fit.
	huvec_fits <- contrasts.fit(fit, contrast.matrix)
	huvec_ebFit <- eBayes(huvec_fits)
# return the top 10 results for any given contrast
# coef=1 is huvec_choroid, coef=2 is huvec_retina
	ttab <- topTable(huvec_ebFit, number=10000, coef=1)

	nrow(topTable(huvec_ebFit, coef=1, number=10000, lfc=5))
	nrow(topTable(huvec_ebFit, coef=1, number=10000, lfc=4))
	nrow(topTable(huvec_ebFit, coef=1, number=10000, lfc=3))
	nrow(topTable(huvec_ebFit, coef=1, number=10000, lfc=2))
# Get a list for probesets with a four fold change or more
	tTable <- topTable(huvec_ebFit, coef=1, number=10000, lfc=4)

	annotation <- as.data.frame(select(hgu133plus2.db,
																		 rownames(tTable), 
																		 c("ENSEMBL","SYMBOL")))
	colnames(annotation) <- c("probeId","ensemblId","geneSymbol")
	results <- merge(annotation, tTable,by.x="probeId",by.y="row.names")

	head(results)
	write.table(results, "results.txt", sep="\t", quote=FALSE)
	return(results)
}

if(!exists("celResults"))
{
	celRAW <- .getData()
	eset<-exprs(celRAW)
	celGCRMA <- gcrma(celRAW)
	celRMA <- rma(celRAW)
	.ggboxIt(.summariseIt(log2(exprs(celRAW))),"sumRAW.pdf")
	.ggboxIt(.summariseIt(exprs(celRMA)),"sumRMA.pdf")
	.plotDensity(log2(exprs(celRAW)),"densityRAW.pdf")
	.plotDensity(log2(exprs(celRMA)),"densityRMA.pdf")
	celFilt <- .doFilter(celRMA)
	celResults <- .doDE(celFilt$eset)

}

