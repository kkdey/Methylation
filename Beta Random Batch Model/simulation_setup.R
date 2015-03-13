

##############  Methylation Data Simulation (Random Batch) #######################


K=4;
G=100;
N=100;

alpha_true=matrix(rnorm(K*G,0,1),nrow=(K)); ### the matrix of fixed effects


Label.Batch=c(rep(1,N/2),rep(2,N/2));

B=max(Label.Batch);

sigmab_true=2;

beta_true=matrix(0,B,G);       ###  the matrix of the random effect

for(g in 1:G)
{
	beta_true[,g]=rnorm(B,mean=0,sd=sigmab_true);
}


phi_true=rnorm(G,100,1);

library(gtools)
T=2;
omega_true=matrix(rbind(rdirichlet(T*10,c(3,4,2,6)),rdirichlet(T*10,c(1,4,6,3)),
			rdirichlet(T*10,c(4,1,2,2)),rdirichlet(T*10,c(2,6,3,2)),
			rdirichlet(T*10,c(3,3,5,4))), nrow=N);



over_dis=0.5;

noise_true=matrix(0,N,G);

for(n in 1:N)
{
	noise_true[n,]=over_dis*rnorm(G,0,1);
}



###  generating the table 


beta_values=matrix(0,N,G);

for(n in 1:N)
{
	for(g in 1:G)
	{

		mean=inv.logit(omega_true[n,]%*%(alpha_true[,g]+beta_true[Label.Batch[n],g]+noise_true[n,g]));
		beta_values[n,g]=rbeta(1,shape1=mean*phi_true[g],shape2=(1-mean)*phi_true[g],ncp=0);
	}
}

plot(density(beta_values),col="red",type="l");

windows()
par(mar=c(8,5.25,2.5,2.5))

	# - get rid of space space between leftmost bar and y axis
par(xaxs="i")
k=K
barplot(t(omega_true),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)



######################   The  End   ######################################

