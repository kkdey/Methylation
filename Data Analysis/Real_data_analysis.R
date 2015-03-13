

#############  Real  Data Analysis #################################

sample=order(rowSums(numData),decreasing=T)[1:1000];
numData_reduced=numData[sample,];
betas=t(numData_reduced);


###  now run the remaining part of the Methyl.estimation.R file 

###  The number of clusters in the real Data -  by Principal Component Analysis

###  Probably from the PC plots (PC1 vs PC2) or the Circular Projection plot, we

###   can get some intuitive idea


Princomp.analysis=prcomp(betas,center=TRUE);
eigenPrincomp=Princomp.analysis$sdev;
rotPrincomp=Princomp.analysis$rotation;
windows()
barplot(eigenPrincomp,col="red", 
		ylab="eigenvalues of beta value matrix",
		xlab="Labels of clusters",xaxt="n");
#axis(1,at=0.7+1.2*(0:6),c("PC1","PC2","PC3","PC4","PC5","PC6","PC7"));

projection=t(rotPrincomp[,-dim(rotPrincomp)[2]])%*%t(betas);

windows()
plot(projection[1,],projection[2,],main="FirstPC vs SecondPC plot",xlab="First PC",
	ylab="Second PC",col=c(rep("red",17),rep("blue",40)));
legend("bottomright",c("N samples", "T samples"),fill=c("red","blue"));

windows()
plot(projection[1,],projection[3,],main="FirstPC vs ThirdPC plot",xlab="First PC",
	ylab="Third PC",col=c(rep("red",17),rep("blue",40)));
legend("bottomright",c("N samples", "T samples"),fill=c("red","blue"));





