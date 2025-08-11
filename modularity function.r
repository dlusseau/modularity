modularity<-function(AI)
{

################################################################################
### This function defines clusters or communities in networks using the notion #
### of modularity first develop by Michelle Girvan and Mark Newman #############
### It implements the algorithm described in Newman 2006 PNAS 103:8577 to ######
### define clusters based on a modularity matrix of the network. This algorithm#
### can be used with weighted and directed networks. See Lusseau et al. 2008 ###
### Animal Behaviour 75:1809 and Leicht & Newman 2008 Physics Review Letters ###
### 100: 118703 for further details. ###########################################
################################################################################

################################################################################
###code originally developed by Hal Whitehead and David Lusseau in Matlab####### 
###for socprog (http://myweb.dal.ca/hwhitehe/social.htm) #######################
###code adapted in R by David Lusseau and Marianne Marcoux #####################
### last changes in May 2012 by David Lusseau to use Fast partial SVD ##########
### by implicitly-restarted Lanczos bidiagonalization to implement the #########
### algorithm on large matrices. This requires the library irlba. ##############
################################################################################

################################################################################
### the input, AI, is a  n x n square matrix that represents interactions in ###
### a network. Those can be weighted. The output is the modularity coefficient #
### (ModCoef) and a vector (clust) of size n which represents the clusters. ####
################################################################################

################################################################################
### 20 Feb 2017: change in lines 58-59: when useind is a sole node, R change the 
### nature of AI and bq to vector, now coercing back to matrix and for bq, 
### if loop added to avoid empty matrix
################################################################################

#library(irlba)

n<-dim(AI)[1]

#calculates modularity matrix
strength<-rowSums(AI)

rAI<-matrix(0,n,n)
rAI<-(strength%*%t(strength))/(sum(strength))
Qq<-AI-rAI

bq<-Qq
numclus<-1
clusass<-array(1,n)
clusfin<-array(numclus,1)
clusfin[numclus]<-1

while (1) {
 uq<-which(clusfin>0)

 if (length(uq)==0) {
  break
 }

 useind<-which(clusass==uq[1])
 bq<-Qq[useind,useind]
 kg<-colSums(as.matrix(AI[useind,useind]))-strength[useind]*sum(strength[useind])/sum(strength)
 if (length(useind)>1) {
 bq<-bq-diag(kg)        }

# cln<-irlba(as.matrix(bq), nu=0, nv=1)$v[,1] # new code with irlba for very large networks
 cln<-eigen(as.matrix(bq), symmetric=TRUE)$vectors[,1]
 clna<-(cln>0)
 nw<-length(clna)

 while (1) {
  qqq<-clna*matrix(1,1,nw)
  qqq
  clna
  modu<-sum(bq*(matrix(qqq,length(qqq),length(qqq))==t(matrix(qqq,length(qqq),length(qqq)))))

  sw<-1
  for (i in 1:length(clna)) {
    clnu<-clna
    clnu[i]=!clna[i]
    qqqu<-clnu*matrix(1,1,nw)
     moduu<-sum(bq*(matrix(qqqu,length(qqqu),length(qqqu))==t(matrix(qqqu,length(qqqu),length(qqqu)))))
    if (moduu>modu) {
     sw<-0
     clna<-clnu
     modu<-moduu
    }
  }

  if (sw) break
}

 if (sum(clna)*sum(!clna)) {
    numclus<-numclus+1
    clusass[useind[which(clna!=0)]]<-numclus
    clusfin[numclus]<-1
  }
  else {
    clusfin[uq[1]]<-0
  }
}

newunit<-clusass

qqq<-clusass*array(1,n)
moduf<-Qq*(matrix(qqq,length(qqq),length(qqq))==t(matrix(qqq,length(qqq),length(qqq))))

diag(moduf)<-0
Q<-sum(moduf)/sum(strength)


list(ModCoef=Q,clust=newunit)
}

