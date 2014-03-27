#install+load required packages into R studio for script to function
install.packages("ape")
install.packages("phangorn")
install.packages("picante")
require(ape)
require(phangorn)
require(picante)
#reads in the selected .fasta file (MMEL1 in this case), assigns data to variable
x<-read.dna("MMEL1first10aligned.fasta",format="fasta")
#assigns variable to store calcuated distances for sequences in .fasta file
d<-dist.dna(x)
#stores results as matrix and exports it to output file
write.table(as.matrix(d),"distances.csv")
#use a mutation model and recalculate distances
d2<-dist.dna(x,model="K80")
#use neighbour join algorithim to work out tree positions
tr.nj<-nj(d)
#plots tree
plot(tr.nj)
#use cophenetic analysis to calculate tree distortion
dt.nj<-cophenetic(tr.nj)
dmat<-as.matrix(d)
nms<-rownames(dmat)
dt.nj<-dt.nj[nms, nms]
dt.nj<-as.dist(dt.nj)
plot(dt.nj-d,ylab="residuals", cex=0.5,main="NJ")
#tr.upgma<-upgma(d)
abline(h=0,lty=3) 
fit<-pml(tr.nj,as.phyDat(x))
fit=optim.pml(fit,T)
plot(fit)
set.seed(8)
#bpptstrapping of data 100 times
bs<-bootstrap.pml(fit,bs=100,optNni=T)
treeBS<-plotBS(fit$tree, type="p", bs)
mt<-modelTest(as.phyDat(x),G=F,I=F)
