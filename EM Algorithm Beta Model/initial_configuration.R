
meth=read_meth

mean_temp=array(0,c(B,K,G));
scale=1; K=4;

omega_preprocess=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);

# omega_preprocess=matrix(1/K,N,K);

iteration=1;
MaxIter=1000
diff2=1000; diff1=20;

logit_meth=log(meth)-log(1-meth);



# for(b in 1:B)
# {
	
	meth_batch=meth[which(Label.Batch==b),];
	logit_meth_batch=log(meth_batch)-log(1-meth_batch);
	# pnorm_meth_batch=qnorm(meth_batch);
	scale=1; K=4;

	#omega0=matrix(rdirichlet(dim(counts_batch)[1],
	#		c(scale/K,scale/K,scale/K,scale/K)), nrow=dim(counts_batch)[1]);

	# omega0=omega_true[Label.Batch==b,];

	omega0=omega_preprocess[which(Label.Batch==b),];


	iteration=1;
	MaxIter=100
	diff2=1000; diff1=100;

	while(iteration < MaxIter)
	{

		####   Estimation of the matrix H

		svd_omega=svd(omega0);
		temp1=t(svd_omega$v)%*%diag(1/svd_omega$d^2,dim(omega0)[2])%*%svd_omega$v;
		temp2=t(omega0)%*%logit_meth_batch;
		# temp2=t(omega0)%*%pnorm_meth_batch;
		temp1=solve(t(omega0)%*%omega0);
		H = temp1%*%temp2;

		###  Estimation of the matrix W (or omega) 

		omega=matrix(0,dim(logit_meth_batch)[1],K);
		for(n in 1:dim(logit_meth_batch)[1])
		{
			omega_vec=omega0[n,];
			meth_vec=logit_meth_batch[n,];
			res=optim(reverse_transform(omega_vec), function(v) loglik_norm(transform(v),meth_vec,t(H)) );
			omega[n,]=transform(res$par);
		}
		diff2=diff1;
		diff1=fnorm(logit_meth_batch,omega%*%H);
		cat("The difference is",diff1,"\n");
		omega0=omega;
		iteration=iteration+1;
	}

	mean_temp[b,,]=H;
	omega_preprocess[which(Label.Batch==b),]=omega0;
}

for(k in 1:K)
{
	for(g in 1:G)
	{
		alpha[k,g]=mean(mean_temp[,k,g]);
	}
}

docweights=omega_preprocess;
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

omega_preprocess=docweights;

alpha_initial=alpha[perm_set[p_star,],];




windows()
barplot(t(omega_preprocess),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


omega_initial=omega_preprocess


###  we carry from this code "omega_initial" and alpha_initial" to the "estimation_EM.R"
###   code as the initial values

alpha0=alpha_initial;
omega0=omega_preprocess

for(g in 1:G)
{
	res=optim(alpha0[,g], function(v) meth_loglik_EM_alpha(v,meth,omega0,phi_true,g),
		method="Nelder-Mead")
	alpha0[,g]=res$par;
}

alpha_initial=alpha0;




