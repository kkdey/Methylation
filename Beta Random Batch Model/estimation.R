

#############  Estimation of Random Effect Beta Topic Model #################


omega_loglik_beta_random = function(omega_vec,meth_vec,alpha,beta,noise,phi)
{

	G=length(meth_vec);
	sum=0;
	for(g in 1:G)
	{
		mean=inv.logit(omega_vec%*%(alpha[,g]+beta[Label.Batch[n],g]+noise[n,g]));
		#cat(lambda,"\n")
		sum=sum - dbeta(meth_vec[g],shape1=mean*phi[g],shape2=(1-mean)*phi[g],log=TRUE) ;
		# cat(sum,"\n");
	}

	return(sum);
}


#################  Initialization of the Beta values  #######################

meth=read_meth;
scale=1; K=4; N=dim(meth)[1]; G=dim(meth)[2];
omega0=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);
library(betareg)
library(lme4)
alpha=matrix(0,K,G);


#############  True log likelihood of the beta values  #####################


true_loglik=0;
for(n in 1:dim(meth)[1])
{
	omega_vec=omega_true[n,];
	meth_vec=meth[n,];
	res=optim(reverse_transform(omega_vec), function(v) omega_loglik_beta_random(transform(v),
					meth_vec,alpha_true,beta_true,noise_true,phi_true) );
	true_loglik=true_loglik+res$value+(scale/K)*sum(log(omega_true[n,]));
}

#############  Estimation of the Model ##############################


iteration=1;
MaxIter=500;

phi=array(0,G);

while(iteration <= MaxIter)
{
	for(g in 1:G)
	{
		meth_col=meth[,g];
		fit=betareg(meth_col~omega0[,1]+omega0[,2]+omega0[,3]+omega0[,4]-1);
		alpha[,g]=as.numeric(fit$coef[1]$mean);
		phi[g]=as.numeric(fit$coef[2]);
	}

	omega=matrix(0,N,K);
	curr_loglik=0;
	for(n in 1:dim(meth)[1])
	{
		omega_vec=omega0[n,];
		meth_vec=meth[n,];
		res=optim(reverse_transform(omega_vec), function(v) omega_loglik_beta(transform(v),meth_vec,alpha,phi) );
		omega[n,]=transform(res$par);
		curr_loglik=curr_loglik+res$value+(scale/K)*sum(log(omega[n,]));
	}
	
	omega0=omega; iteration=iteration+1;
	#diff=(curr_loglik-true_loglik)/true_loglik;
	diff=curr_loglik;
	cat("The relative diff is",diff,"\n");
}


docweights=omega0;
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


windows()
par(mar=c(8,5.25,2.5,2.5))

# - get rid of space space between leftmost bar and y axis
par(xaxs="i")

k=K
# Make plot 
# - plot data
barplot(t(docweights),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


abline(v=17)
axis(1, at=1:57, labels=names(Methyl_Data)[2:58],las=2)


windows()
barplot(t(omega_true),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)




