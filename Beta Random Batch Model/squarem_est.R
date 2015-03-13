
meth=read_meth
y=matrix(meth,1,N*G);


omega_loglik_beta = function(omega_vec,meth_vec,alpha,phi)
{

	G=length(meth_vec);
	sum=0;
	for(g in 1:G)
	{
		mean=inv.logit(omega_vec%*%alpha[,g]);
		#cat(lambda,"\n")
		sum=sum - dbeta(meth_vec[g],shape1=mean*phi[g],shape2=(1-mean)*phi[g],log=TRUE) ;
		# cat(sum,"\n");
	}

	return(sum);
}


scale=1; K=4; N=500; G=100;

####  The starting value of omega (topic prop weights)  ####################

# Use a preset seed so the example is reproducable.
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
seed=54321))

# omega0=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);
omega0=omega_initial;
alpha0=matrix(rnorm(K*G,0,1),nrow=(K)); ### the matrix of fixed effects
phi0=rnorm(G,10,1);
param0=c(matrix(omega0,1,N*K),matrix(alpha0,1,K*G),phi0);


library(betareg)
library(lme4)
alpha=matrix(0,K,G);

library(SQUAREM)

beta_est.probe <- function(param_vec,y)
{
	meth=matrix(y,N,G);
	omega0=matrix(param_vec[1:(N*K)],N,K);

	############  estimation of alpha (effect) and phi (disp) #############

	alpha=matrix(0,K,G);
	phi=array(0,G);
	for(g in 1:G)
	{
		meth_col=meth[,g];
		fit=betareg(meth_col~omega0[,1]+omega0[,2]+omega0[,3]+omega0[,4]-1);
		alpha[,g]=as.numeric(fit$coef[1]$mean);
		phi[g]=as.numeric(fit$coef[2]);
	}

	omega=matrix(0,N,K);
	for(n in 1:dim(meth)[1])
	{
		omega_vec=omega0[n,];
		meth_vec=meth[n,];
		res=optim(reverse_transform(omega_vec), function(v) omega_loglik_beta(transform(v),meth_vec,alpha,phi) );
		omega[n,]=transform(res$par);
	}

	param_vec_omega=matrix(omega,1,N*K);
	param_vec_alpha=matrix(alpha,1,K*G);
	param_new_vec=c(param_vec_omega,param_vec_alpha,phi);
	param_vec=param_new_vec
	return(param_new_vec)
}


options(warn=-1)
system.time(res <- squarem(p=param0,y=y, fixptfn=beta_est.probe, control=list(maxiter = 10, trace = FALSE)));

omega_final=matrix(res$par[1:(N*K)],N,K);
alpha_final_vec=res$par[-(1:(N*K))];
alpha_final=matrix(alpha_final_vec[1:(K*G)],K,G);
phi_final=tail(alpha_final_vec,G);


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




