library(data.table)
library(igraph)

### 0. Function of Clustering Coefficient
ClustBCG_C_star <- function(W) {
  p <- nrow(W)
  A <- (W != 0) * 1
  W_plus_Wt <- W + t(W)
  A_plus_At <- A + t(A)
  s_tot <- rowSums(W_plus_Wt)
  d_tot <- rowSums(A_plus_At)
  s_bilat <- 0.5 * rowSums(W %*% A + A %*% W)
 
  A_plus_At_sq <- A_plus_At %*% A_plus_At
  Num <- 0.5 * diag(W_plus_Wt %*% A_plus_At_sq)
  Den <- s_tot * (d_tot - 1) - 2 * s_bilat
  C_star <- ifelse(Den > 0, Num / Den, 0)
  
  if (!is.null(rownames(W))) {
    names(C_star) <- rownames(W)
  } else {
    names(C_star) <- paste0("Gene_", 1:p)
  }
  return(C_star)
}
############ 

### 1. Data infile
# Expression levels of genes: 
EXP<-read.table("ToyDATA_CellGNEA\\EXP.csv",sep=",")
# Genes involved in a specific pathway:
PW_genes<-read.table("ToyDATA_CellGNEA\\PathwayGENES.csv",sep=",")
# Modulator values in cell lines:
Modulator<-read.table("ToyDATA_CellGNEA\\Modulator.csv",sep=",")

T_col<-ncol(EXP)
p<-T_col
NoCELL<-nrow(Modulator)
gene_names<-colnames(EXP)

SCORE<-matrix(numeric(T_col*NoCELL),ncol=T_col)
colnames(SCORE)<-gene_names; rownames(SCORE)<-rownames(Modulator)

CORR<- matrix(numeric(T_col*1),ncol=1)
rownames(CORR)<-gene_names

### 2. Computing gene score
for (i in 1:NoCELL){

# 2.1 Infile edge weight matrix of i-th cell line:
EdgeW<-read.table(paste("ToyDATA_CellGNEA\\BETA_Sample",i,".csv",sep=""),sep=",")

# 2.2 Clustering Coefficient in (5):
W<-EdgeW
W[]<-0
W[upper.tri(W)] <- EdgeW[upper.tri(EdgeW)]
W[lower.tri(W)] <- EdgeW[upper.tri(EdgeW)]
C_star_coefficients <- ClustBCG_C_star(as.matrix(W))
if ((max(C_star_coefficients)-min(C_star_coefficients))!=0){
nC_star_coefficients<-(C_star_coefficients-min(C_star_coefficients))/(max(C_star_coefficients)-min(C_star_coefficients))
}else if ((max(C_star_coefficients)-min(C_star_coefficients))==0){
nC_star_coefficients<-(C_star_coefficients-min(C_star_coefficients))
}

# 2.3 PageRank in (6):
g_w <- graph_from_adjacency_matrix(as.matrix(abs(W)), mode = "directed", weighted = TRUE)
pr_w <- page.rank(g_w, damping = 0.85, weights = E(g_w)$weight)
PageRank<-pr_w$vector
nPageRank<-(PageRank-min(PageRank))/(max(PageRank)-min(PageRank))
names(nPageRank)<- gsub("\\.", "-", names(nPageRank))
quantile(nPageRank)

# 2.4 Regulatory Effect in (8):
estNP<-c()
for (c in 1:ncol(EdgeW)){
if (sum(EdgeW[,c]!=0)>0){
estNP<-rbind(estNP,cbind(rownames(EdgeW)[EdgeW[,c]!=0],colnames(EdgeW)[c],EdgeW[,c][EdgeW[,c]!=0]))
}
}
estNP<-cbind(estNP,abs(as.numeric(estNP[,3])))
colnames(estNP)<-c("RG","TG","COEF","absCOEF")
rownames(estNP)<-paste(estNP[,1],"_",estNP[,2],sep="")
estNP<-cbind(estNP,t(as.numeric(estNP[,"absCOEF"])*EXP[rownames(Modulator)[i],estNP[,1]]))
colnames(estNP)[5]<-"RE"
dt <- as.data.table(estNP)
result <- dt[, .(RE_sum = sum(as.numeric(RE), na.rm = TRUE)), by = RG]
mat_result <- as.matrix(result)
RE<-matrix(numeric(p*1),ncol=1)
rownames(RE)<-gene_names
RE[mat_result[,1],1]<-as.matrix(as.numeric(mat_result[,2]),ncol=1)

# 2.5 Compute gene score in (9):
SCORE[i,rownames(RE)]<-RE*(nC_star_coefficients[rownames(RE)]+nPageRank[rownames(RE)])*0.5
} # close i in 1:NoCELL

### 3. Correlation coecfficient between gene's score and modulator value in (10):
CORR[apply(SCORE!=0,2,sum)!=0,1] <- apply(SCORE[,apply(SCORE!=0,2,sum)!=0], 2, function(x) cor(x, Modulator))

### 3. Enrichment score 
# 3.1 Compute enrichment score
SC_T<-matrix(numeric(1*(NOpm)),ncol=1)
for (pm in 1:NOpm){
rdPW_genes<-sample(rownames(CORR),length(PW_genes))
erSCORE<-matrix(numeric(length(CORR[,1])*6),ncol=6)
rownames(erSCORE)<-rownames(CORR)
colnames(erSCORE)<-c("R","IN","NiN","Hit","Miss","SC")
erSCORE[,1:1]<-CORR[,1]
erSCORE<-erSCORE[order(erSCORE[,1],decreasing=TRUE),]
if (pm==1){
erSCORE[is.na(match(rownames(erSCORE),PW_genes))==FALSE,2]<-1
erSCORE[is.na(match(rownames(erSCORE),PW_genes))==TRUE,3]<-1
} else if (pm>1){
erSCORE[is.na(match(rownames(erSCORE),rdPW_genes))==FALSE,2]<-1
erSCORE[is.na(match(rownames(erSCORE),rdPW_genes))==TRUE,3]<-1
}
N<-nrow(erSCORE)
Nh<-sum(erSCORE[,2])
NR<-sum(erSCORE[,2]) 
if (NR!=0){
erSCORE[,4]<-cumsum(erSCORE[,2])/NR
}
erSCORE[,5]<-cumsum(erSCORE[,3])/(N-Nh)
erSCORE[,6]<-erSCORE[,4]-erSCORE[,5]
SC_T[pm ,1]<-erSCORE[,6][order(abs(erSCORE[,6]),decreasing=TRUE)[1]]
}

# 3.2 Normalization of the enrichment score
nSC_T<-SC_T
nSC_T[]<-0
nSC_T[SC_T[,1]>0,1]<-SC_T[SC_T[,1]>0,1]/mean(abs(SC_T[2:NOpm,][SC_T[2:NOpm,1]>0]))
nSC_T[SC_T[,1]<0,1]<-SC_T[SC_T[,1]<0,1]/mean(abs(SC_T[2:NOpm,][SC_T[2:NOpm,1]<0]))

# 3.3 Compute p.value
if (nSC_T[1,1]>0){Pvalue<-sum(nSC_T[1,1]<=nSC_T[-1,1][nSC_T[-1,1]>0])/(NOpm-1)};
if (nSC_T[1,1]<0){Pvalue<-sum(abs(nSC_T[1,1])<=abs(nSC_T[-1,1][nSC_T[-1,1]<0]))/(NOpm-1)};

obs <- nSC_T[1,1]
perm <- nSC_T[-1,1]
Pvalue <- sum(sign(perm) == sign(obs) &  abs(perm) >= abs(obs)) / (NOpm - 1)

print(Pvalue)

