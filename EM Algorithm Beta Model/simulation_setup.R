
#################   Assignment of the dimensionality  ##################

K=4;
G=100;
N=500;


#################  Defining the variables of interest  ###################


alpha_true=matrix(rnorm((K)*G,1,1),nrow=(K)); ### the matrix of fixed effects

Label.Batch=c(rep(1,N/2),rep(2,N/2));

phi_true=abs(rnorm(G,100,1));


B=max(Label.Batch);

sigmab_true=0;

beta_true=matrix(0,B,G);       ###  the matrix of the random effect

for(g in 1:G)
{
	beta_true[,g]=rnorm(B,mean=0,sd=sigmab_true);
}

library(gtools)
T=10;
omega_true=matrix(rbind(rdirichlet(T*10,c(3,4,2,6)),rdirichlet(T*10,c(1,4,6,3)),
			rdirichlet(T*10,c(4,1,2,2)),rdirichlet(T*10,c(2,6,3,2)),
			rdirichlet(T*10,c(3,3,5,4))), nrow=N);



over_dis=0;

noise_true=matrix(0,N,G);

for(n in 1:N)
{
	noise_true[n,]=over_dis*rnorm(G,0,1);
}


read_meth=matrix(0,N,G);
indicator_true=matrix(0,N,G);

for(n in 1:N)
{
	for(g in 1:G)
	{
		index=sample(1:K,1,omega_true[n,],replace=T);
		mean=inv.logit(alpha_true[index,g] +beta_true[Label.Batch[n],g]+noise_true[n,g]);
		read_meth[n,g]=rbeta(1,shape1=mean*phi_true[g],shape2=(1-mean)*phi_true[g],ncp=0);
		indicator_true[n,g]=index;
	}
}

windows()
k=K
barplot(t(omega_true),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


barplot(t(omega_preprocess),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
