key <- read.delim("/Users/justinvaughn/Desktop/historicUSRice/masterKey.txt")
#library("gplots")
KI <- read.delim(file="/Users/justinvaughn/Desktop/historicUSRice/popStructure/midsouthCenterIBSKin2AlleleMax.txt",header=FALSE)
#KI <- read.delim(file="/Users/justinvaughn/Desktop/historicUSRice/popStructure/rawInverseDistance.txt",header=FALSE)
theKin=as.matrix(KI[,-1])
rownames(theKin)=KI[,1]
colnames(theKin)=KI[,1]
key <- key[which(key$TruncLocation == "MidSouth"),]
key <- key[which(key$ID != "US_0085"),]
#nrow(key)
#length(rownames(theKin))

out = paste(c("crossGroupAve","g1Ave","g1Var","g2Ave","g2Var"),collapse="\t")
cat(paste(out,"\n"))

all <- as.character(key$ID[order(key$ReleaseDate,decreasing = F)])
#print(all)
#totalMarkers <- 84097
library(mclust)
window = 25
for (i in 1:(length(all)-window)) {
	#print(all[i:(i+window)])
	slice = all[i:(i+window)]
	kin <- theKin[slice,slice] 
	#print(kin)
	#kin[lower.tri(kin)] <- NA
	#diag(kin) <- NA
	#v = as.vector(kin)
	#v = v[!is.na(v)]
	#print(v)
	x.gmm = summary(Mclust(kin,G=2))
	
	kin[lower.tri(kin)] <- NA
	diag(kin) <- NA
	v = as.vector(kin)
	v = v[!is.na(v)]
	Ave = mean(v)
	var = sd(v)
	#print(str(x.gmm))
	group1 <- slice[which(x.gmm$classification == 1)]
	group2 <- slice[which(x.gmm$classification == 2)]
	
	#print(group1)
	
	kin <- theKin[group1,group2] 
	#print(kin)
	#kin[lower.tri(kin)] <- NA
	#diag(kin) <- NA
	v = as.vector(kin)
	v = v[!is.na(v)]
	crossGroupAve = mean(v)
	
	kin <- theKin[group1,group1] 
	kin[lower.tri(kin)] <- NA
	diag(kin) <- NA
	v = as.vector(kin)
	v = v[!is.na(v)]
	g1Ave = mean(v)
	g1Var = sd(v)
	
	kin <- theKin[group2,group2] 
	kin[lower.tri(kin)] <- NA
	diag(kin) <- NA
	v = as.vector(kin)
	v = v[!is.na(v)]
	g2Ave = mean(v)
	g2Var = sd(v)
	
	
	if (g1Ave < g2Ave) {
		temp = g1Ave
		g1Ave = g2Ave
		g2Ave = temp
		temp = g1Var
		g1Var = g2Var
		g2Var = temp
	}
	#intraGroupAve = mean(c(v1,v2))
	#intraGroupVar = var(c(v1,v2))
	
	out = paste(c(crossGroupAve,g1Ave,g1Var,g2Ave,g2Var),collapse="\t")
	cat(paste(out,"\n"))
	
# 	print(paste(as.vector(summary(x.gmm)$mean),collapse="_"))
# 	x.gmm$classification
# 		#print(x.gmm)
# 	l = length(summary(x.gmm)$mean)
# 	print(l)
		# means = paste(as.vector(summary(x.gmm)$mean),collapse="_")
# 		cat(paste(rownames(countData)[i],l,means,sep="\t"),"\n") #for some reason cycles through each summary value
# 		allele = paste(c(rownames(countData)[i],l,x.gmm$classification),collapse="\t")
		#cat(paste(allele,"\n"))
		#}
}



