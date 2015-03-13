

###############  The Squarem Implementation (Final Step) ##################


meth=read_meth
y=matrix(meth,1,N*G);

scale=1; K=4; N=500; G=100;

####  The starting value of omega (topic prop weights)  ####################

# Use a preset seed so the example is reproducable.
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
seed=54321))

omega0=omega_initial;
alpha0=alpha_initial;
param0=c(matrix(omega0,1,N*K),matrix(alpha0,1,K*G));


library(betareg)

beta_final.em <- function(param_vec,y)
{
	meth=matrix(y,N,G);
	omega0=matrix(param_vec[1:(N*K)],N,K);
	alpha0=matrix(param_vec[-(1:(N*K))],K,G);
	

	###############   Estimation of omega  ############################

	omega=matrix(0,N,K);
	Z=array(0,c(N,K,G));

	for(n in 1:N)
	{
		for(g in 1:G)
		{
			for(k in 1:K)
			{
				mu=inv.logit(alpha0[k,g]);
				mu[mu<1e-04]=1e-04;
				mu[mu>(1-1e-04)]=(1-1e-04);

				Z[n,k,g]=omega0[n,k]*dbeta(meth[n,g],
						shape1=mu*phi_true[g],
						shape2=(1-mu)*phi_true[g],ncp=0);
			}
			Z[n,,g]=(Z[n,,g]+1e-07)/sum(Z[n,,g]+1e-07)
		}
	}

	for(n in 1:N)
	{
		for(k in 1:K)
		{
			omega[n,k]=sum(Z[n,k,])/G;
		}
	}

	##################  Defining the label matrix ###################

	indicator=matrix(0,N,G);
	for(n in 1:N)
	{
		for(g in 1:G)
		{
			indicator[n,g]=sample(1:K,1,(Z[n,,g]/sum(Z[n,,g])),replace=TRUE);
		}
	}
	

	##################  Estimation of alpha  ########################

	alpha=alpha0;

	for(g in 1:G)
	{
		meth_col=meth[,g];
		indicator_col=indicator[,g];
		fit=betareg(meth_col~as.factor(indicator_col)-1);
		alpha[sort(unique(indicator_col)),g]=as.numeric(fit$coef[1]$mean);
		# phi[g]=as.numeric(fit$coef[2]);
	}

	################## Pooling everything back in vector ###########

	param_vec_omega=matrix(omega,1,N*K);
	param_vec_alpha=matrix(alpha,1,K*G);
	param_new_vec=c(param_vec_omega,param_vec_alpha);
	param_vec=param_new_vec
	return(param_new_vec)
}


options(warn=-1)
system.time(res <- squarem(p=param0,y=y, fixptfn=beta_final.em, control=list(maxiter = 30, trace = FALSE)));

##################  Adjusting for labeling and plotting ####################

omega_final=matrix(res$par[1:(N*K)],N,K);
alpha_final=matrix(res$par[-(1:(N*K))],K,G)

docweights=omega_final;
library(permute);
library("BioPhysConnectoR");
perm_set=rbind(1:K,allPerms(1:K));
diff=array(0,dim(perm_set)[1]);
for (p in 1:dim(perm_set)[1])
{
	temp=docweights[,perm_set[p,]];
	diff[p]=fnorm(temp,omega_true);
}

p_star=which(diff==min(diff));
docweights=docweights[,perm_set[p_star,]];

omega_final=docweights;
alpha_final=alpha_final[perm_set[p_star,],];

windows()
k=K
barplot(t(omega_final),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

windows()
k=K
barplot(t(omega_true),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


#######################   The End   #####################################


		


