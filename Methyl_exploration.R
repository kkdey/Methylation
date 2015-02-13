
library(data.table)

Methyl_Data=read.table("D:/Matthew_Stephens_Project/Methylation/Methylation Expt 2/betas_file_copy_reduced.txt",header=T);

numData=as.matrix(Methyl_Data[-1,-1]);
plot(density(numData),col="red",xlab="Beta values",ylab="Density",main="");
title(main="Density plot of the beta values");

batch=c(rep(1,18),rep(2,(dim(numData)[2]-18))); batch=as.factor(batch);
log_pval=array(0,dim(numData)[1]);

for(l in 1:dim(numData)[1])
{
	log_pval[l]=-log(ks.test(numData[l,which(batch==1)],numData[l,which(batch==2)])$p.value);
}

plot(1:dim(numData)[1],log_pval,col="red",main="- log p-values of CpGs between tissues",xlab="",ylab="");
title(xlab="CpG sites (indexed)");
title(ylab="-log p-value of KS test");

plot(1:500,log_pval[1:500],col="red",type="l",main="- log p-values of 500 CpGs between tissues",xlab="", ylab="");
title(xlab="CpG sites (indexed)");
title(ylab="-log p-value of KS test");

###############  Sample to sample heat maps #############################

sample_cov= t(numData)%*%numData;
heatmap(sample_cov,col = cm.colors(256));

    library(gplots) # for colour panel of heatmap
    library(heatmap.plus)

    heatmap(sample_cov,
    Rowv=FALSE, Colv=FALSE,
    cexCol=1, cexRow=1,
    key=TRUE, keysize=0.1, # display colour key
    density.info=c("none"), # options of different plots to be drawn in colour key
    trace="none", # character string indicating whether a solid “trace” line should be drawn across ‘row’s or down ‘column’s, ‘both’ or ‘none’.
    margins=c(5,9),
    lmat=rbind( c(0,3), c(2,1), c(0,4) ), lhei=c(0.2, 8.5, 2), # where to display colour key
    col=grey(seq(1,0,-0.01)) # custom colours for colour key
    )

#############  The N and Ts seem to cluster together #####################

############   Model  Based  clustering  #################################

library(betareg)
precision=array(0,dim(numData)[1]); effect_size=array(0,dim(numData)[1]);
sd_res=array(0,dim(numData)[1]);

for(l in 12765:dim(numData)[1])
{
	numData[l,which(numData[l,]==0)]=1e-07;
	numData[l,which(numData[l,]==1)]=(1-1e-07);

	L=betareg(numData[l,]~factor(batch));
	precision[l]=L$coef[2];
	effect_size[l]=as.numeric(L$coef[1]$mean[2])
	sd_res[l]=sd(L$residuals);
}

##  the CpG sites that gave error in above loop are:  6162,6282,12151,12764

bad_CpGs=c(6162,6282,12151,12764);
sd_res=sd_res[-bad_CpGs];
precision=precision[-bad_CpGs];
effect_size=effect_size[-bad_CpGs];

write.table(sd_res,"D:/Matthew_Stephens_Project/Methylation/Methylation Expt 2/Data/sd_residuals.txt");
write.table(precision,"D:/Matthew_Stephens_Project/Methylation/Methylation Expt 2/Data/precision_est.txt");
write.table(effect_size,"D:/Matthew_Stephens_Project/Methylation/Methylation Expt 2/Data/regression_effect_size.txt");

CpGs=1:dim(numData)[1];
pdf("D:/Matthew_Stephens_Project/Methylation/Methylation Expt 2/figures/std_err_residuals.pdf");
plot(CpGs[-bad_CpGs],sd_res,col="red",
		type="l",main="Standard error of residuals for each CpG site",
			xlab="CpG sites",ylab="Std. error of residuals");
dev.off()

pdf("D:/Matthew_Stephens_Project/Methylation/Methylation Expt 2/figures/precision_CpGs.pdf");
plot(CpGs[-bad_CpGs],precision,col="red",
		type="l",main="Precision estimate CpG site",
			xlab="CpG sites",ylab="precision");
dev.off()

pdf("D:/Matthew_Stephens_Project/Methylation/Methylation Expt 2/figures/effect_size_CpGs.pdf");
plot(CpGs[-bad_CpGs],effect_size,col="red",
		type="l",main="Effect size of each CpG site",
			xlab="CpG sites",ylab="effect sizes");
dev.off()

sd_samples=array(0,dim(numData)[2]);

for(lab in 1:dim(numData)[2])
{
	sd_samples[lab]=sd(numData[,lab]);
}

pdf("D:/Matthew_Stephens_Project/Methylation/Methylation Expt 2/figures/sd_error_samples.pdf");
plot(1:dim(numData)[2],sd_samples,col=c(rep("blue",18),rep("red",(57-18))),
	main="Std error of methylation beta values of various samples",
	xlab="sample index", ylab="std error",pch=20,lwd=4,ylim=c(0.3,0.35));
legend("topright",legend=c("N samples","T samples"),fill=c("blue","red"));
dev.off()

#####  This plot implies that there is difference in variance between the 

####  N samples and the T samples 

###  So a Beta Mixed Model would probably make more sense 



mean_samples=array(0,dim(numData)[2]);

for(lab in 1:dim(numData)[2])
{
	mean_samples[lab]=mean(numData[,lab]);
}

pdf("D:/Matthew_Stephens_Project/Methylation/Methylation Expt 2/figures/mean_beta_samples.pdf");
plot(1:dim(numData)[2],mean_samples,col=c(rep("blue",18),rep("red",(57-18))),
	main="Mean val of methylation beta values for various samples",
	xlab="sample index", ylab="means",pch=20,lwd=4);
legend("topright",legend=c("N samples","T samples"),fill=c("blue","red"));
dev.off()


#############  Comparative density plots for N and T samples  ####################
pdf("D:/Matthew_Stephens_Project/Methylation/Methylation Expt 2/figures/comparative_density_plots.pdf");
plot(density(numData),col="red",xlab="Beta values",ylab="Density",main="");
title(main="Density plot of the beta values");
lines(density(numData[,1:18]),col="blue");
lines(density(numData[,19:57]),col="green");
legend("topright",legend=c("ALL samples","N samples","T samples"),fill=c("red","blue","green"));

dev.off()


###################   The  End   #########################################



