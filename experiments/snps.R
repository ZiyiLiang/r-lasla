#----------------------------------------------------------------------------------------
# Name: snps.R
# Description:
# Select 5000 SNPs from chr6 and use ld matrix as the side information. Implement ANEAT
# and compare to BH.
#----------------------------------------------------------------------------------------

source("~/laws/R/neat_funcs_test.R")
## Real Data Application: SNPs


## Simulation 1: 5000 selected Diabete 2 SNPs from chr6

# # read the ld info file
# lddata<-read.delim("F:/gene/myld.ld",sep="")
# 
# # read the summary stats for the selected 5000 SNPs
# testgenom<-read.delim("F:/gene/S5000.txt")
# 
# #calculate the z-scores
# z_snp<-testgenom$Beta/testgenom$SE
# 
# #calculating the distance matrix
# d_snp<-matrix(rep(1, m*m), m, m)
# x_ld<-rep(0,length(lddata$SNP_A))
# y_ld<-rep(0,length(lddata$SNP_B))
# 
# for(i in 1:length(lddata$CHR_A)){
#   cat(i)
#   x_ld[i]<-which(testgenom$Pos==lddata$BP_A[i])
#   y_ld[i]<-which(testgenom$Pos==lddata$BP_B[i])
#   d_snp[x_ld[i],y_ld[i]]=1-lddata$R2[i]
#   d_snp[y_ld[i],x_ld[i]]=d_snp[x_ld[i],y_ld[i]]
# }


m<-5000
mu0<-rep(0,m)
sd0<-rep(1,m)
fdr_level<-c(0.001,0.01,0.05,0.1)
#fdr_level<-c(0.05)
nrep<-length(fdr_level)

de<-matrix(rep(0, m*nrep), nrep, m)

pv<-2*pnorm(-abs(z_snp), 0, 1)

#bon.nr<-rep(0,nrep)
bh.nr<-rep(0,nrep)
aneat.dd.nr<-rep(0,nrep)

for (i in 1:nrep) {
  cat("iteration:", i, "\n")
  q<-fdr_level[i]
  
  #bh.res<-bh.func(pv, q)
  
  bw=density(z_snp,bw="SJ-ste")$bw
  pis_neat <- pis_1D_neat.func(z_snp, d=d_snp, pval=pv, tau=0.1, h=0.5)
  
  weight<-weights_1D.func(z_snp,d_snp,pis_neat,mu0,sd0,ha=bw)
  aneat.dd.res<-neat.func(pvs=pv, pis=pis_neat, ws=weight, q)
  
  #d_none<-matrix(rep(1,m^2),m,m)
  #pis_none <- pis_1D_neat.func(z_snp, d=d_none, pval=pv, tau=0.1, h=0.1)
  #weight_none<-weights_1D.func(z_snp,d_none,pis_none,mu0,sd0,ha=0.1)
  #none.res<-neat.func(pvs=pv, pis=pis_none, ws=weight_none, q)
   
  #bon.nr[i]<-length(which(pv<=q/m))
  #bh.nr[i]<-bh.res$nr
  aneat.dd.nr[i]<-aneat.dd.res$nr
  de[i,]<-aneat.dd.res$de
}