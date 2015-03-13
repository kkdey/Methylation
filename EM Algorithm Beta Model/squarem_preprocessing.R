


meth=read_meth

scale=1; K=4; N=500; G=100;

####  The starting value of omega (topic prop weights)  ####################

# Use a preset seed so the example is reproducable.
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
seed=54321))

omega_preprocess=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);

logit_meth=log(meth)-log(1-meth);
y=matrix(logit_meth,1,N*G);


omega0=omega_preprocess
alpha0=matrix(rnorm((K)*G,1,1),nrow=(K));
param0=c(matrix(omega0,1,N*K),matrix(alpha0,1,K*G));

library(SQUAREM)

beta_process.em <- function(param_vec,y)
{
	logit_meth=matrix(y,N,G);
	#N=dim(logit_meth)[1]; G=dim(logit_meth)[2];
	omega0=matrix(param_vec[1:(N*K)],N,K);
	alpha0=matrix(param_vec[-(1:(N*K))],K,G);
	
	svd_omega=svd(omega0);
	temp1=t(svd_omega$v)%*%diag(1/svd_omega$d^2,dim(omega0)[2])%*%svd_omega$v;
	temp2=t(omega0)%*%logit_meth;
	# temp2=t(omega0)%*%pnorm_meth;
	temp1=solve(t(omega0)%*%omega0);
	H = temp1%*%temp2;

	###  Estimation of the matrix W (or omega) 

	omega=matrix(0,N,K);
	for(n in 1:N)
	{
		omega_vec=omega0[n,];
		meth_vec=logit_meth[n,];
		res=optim(reverse_transform(omega_vec), function(v) loglik_norm(transform(v),meth_vec,t(H)) );
		omega[n,]=transform(res$par);
	}

	param_vec_omega=matrix(omega,1,N*K);
	param_vec_alpha=matrix(H,1,K*G);
	param_new_vec=c(param_vec_omega,param_vec_alpha);
	param_vec=param_new_vec
	return(param_new_vec)
}

options(warn=-1)
system.time(res <- squarem(p=param0,y=y, fixptfn=beta_process.em, control=list(maxiter = 100, trace = FALSE)));

######################################################################################
omega_preprocess=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);

logit_meth=log(meth)-log(1-meth);
maxIter=100;

iteration=1;
omega0=omega_preprocess;
alpha0=matrix(rnorm((K)*G,1,1),nrow=(K));
param0=c(matrix(omega0,1,N*K),matrix(alpha0,1,K*G));
y=matrix(logit_meth,1,N*G);


system.time(
while(iteration<=maxIter)
{
	param0=beta_process.em(param0,y);
	iteration=iteration+1;
}   )

omega_initial=matrix(param0[1:(N*K)],N,K);
alpha_initial=matrix(param0[-(1:(N*K))],K,G);

################################################################################################

omega_initial=matrix(res$par[1:(N*K)],N,K);
alpha_initial=matrix(res$par[-(1:(N*K))],K,G)

docweights=omega_initial;
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

omega_initial=docweights;
alpha_initial=alpha_initial[perm_set[p_star,],];

windows()
k=K
barplot(t(omega_initial),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

######################  Using the alpha log likelihood to modify alpha est  #########################

alpha0=alpha_initial;
omega0=omega_initial; 

for(g in 1:G)
{
	res=optim(alpha0[,g], function(v) meth_loglik_EM_alpha(v,meth,omega0,phi_true,g),
		method="Nelder-Mead")
	alpha0[,g]=res$par;
}

alpha_initial=alpha0;



