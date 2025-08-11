modularity_coefficient<-function(AI,clusass){
#estimate modularity of each topological modules in a network

#clusass is a vector of size n, with cluster membership for each node (1 to k clusters)
#AI is a symmetric matrix n x n

n<-dim(AI)[1]
k<-max(clusass)
#calculates modularity matrix
strength<-rowSums(AI)
rAI<-matrix(0,n,n)
rAI<-(strength%*%t(strength))/(sum(strength))
Qq<-AI-rAI  #modularity matrix
diag(Qq)<-0
	
Q<-sum(Qq)/sum(strength)
mod<-array(0,k)

for (i in 1:k) {
	mod[i]<-sum(Qq[which(clusass==i),which(clusass==i)])/sum(strength[which(clusass==i)])
}
	return(list(mod=mod,Q=Q))
}


