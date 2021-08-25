library(ggplot2)
library(RadialMR)
library(MVMR)
library(RMVMR)

###############################
###   Univariate Analyses   ###
###############################

#Exposure 1

#Select only instruments for exposure 1
X1datuni<-sum.data[sum.data$pgammaX1hat<0.05,]

#Format exposure 1 data
X1f.data<-format_radial(X1datuni[,2], X1datuni[,8], X1datuni[,5], X1datuni[,9], X1datuni[,1])

#Estimate causal effect of exposure 1 on Y
resX1uni<-ivw_radial(X1f.data,0.05,1,0.0001,T)

#Generate univariable radial plot for exposure 1
plotX1<-plotly_radial(resX1uni)


#Exposure 2

#Select only instruments for exposure 2
X2datuni<-sum.data[sum.data$pgammaX2hat<0.05,]

#Format exposure 2 data
X2f.data<-format_radial(X2datuni[,3], X2datuni[,8], X2datuni[,6], X2datuni[,9], X2datuni[,1])

#Estimate causal effect of exposure 2 on Y
resX2uni<-ivw_radial(X2f.data,0.05,1,0.0001,T)

#Generate univariable radial plot for exposure 2
plotX2<-plotly_radial(resX2uni)


#Exposure 3

#Select only instruments for exposure 3
X3datuni<-sum.data[sum.data$pgammaX3hat<0.05,]

#Format exposure 3 data
X3f.data<-format_radial(X3datuni[,4], X3datuni[,8], X3datuni[,7], X3datuni[,9], X3datuni[,1])

#Estimate causal effect of exposure 3 on Y
resX3uni<-ivw_radial(X3f.data,0.05,1,0.0001,T)

#Generate univariable radial plot for exposure 3
plotX3<-plotly_radial(resX3uni)



###########################################
###   Radial MVMR Analyses: X1 and X2   ###
###########################################

# Select instruments for either exposure 1 or exposure 2

X1X2datuni<-sum.data[sum.data$pgammaX1hat<0.05 | sum.data$pgammaX2hat<0.05,]

#Format the X1-X2 data
X1X2f.data<-format_rmvmr(sum.data[,2:3], sum.data[,8], sum.data[,5:6], sum.data[,9], sum.data[,1])

#Estimate conditional instrument strength
strength_X1X2<- strength_rmvmr(X1X2f.data)

#Show conditional F statistic
strength_X1X2$f

#Generate conditional strength radial MVMR plot for exposure 1
strength_X1X2$plot[1]

#Generate conditional strength radial MVMR plot for exposure 2
strength_X1X2$plot[2]

#Estimate causal effects of exposure 1 and exposure 2 using Radial MVMR
X1X2res<-ivw_rmvmr(X1X2f.data,T)

#Estimate pleiotropic effects and detect outliers
pleioX1X2<-pleiotropy_rmvmr(X1X2f.data,X1X2res)

#Show Q-statistics and p-values for each SNP across all exposures
pleioX1X2$qdat

#Show global Q-statistics for all exposures
pleioX1X2$gq

#Generate unadjusted and adjusted radial MVMR plots
X1X2plots<-plot_rmvmr(X1X2f.data,X1X2res)
X1X2plots$p1
X1X2plots$p2



################################################
###   Radial MVMR Analyses: X1, X2, and X3   ###
################################################



# Select instruments for any exposure

X1X2X3datuni<-sum.data

#Format the X1-X2 data
X1X2X3f.data<-format_rmvmr(sum.data[,2:4], sum.data[,8], sum.data[,5:7], sum.data[,9], sum.data[,1])

#Estimate conditional instrument strength
strength_X1X2X3<- strength_rmvmr(X1X2X3f.data)

#Show conditional F statistic
strength_X1X2X3$f

#Generate conditional strength radial MVMR plot for exposure 1
strength_X1X2X3$plot[1]

#Generate conditional strength radial MVMR plot for exposure 2
strength_X1X2X3$plot[2]

#Generate conditional strength radial MVMR plot for exposure 3
strength_X1X2X3$plot[3]

#Estimate causal effects  using Radial MVMR
X1X2X3res<-ivw_rmvmr(X1X2X3f.data,T)

#Estimate pleiotropic effects and detect outliers
pleioX1X2X3<-pleiotropy_rmvmr(X1X2X3f.data,X1X2X3res)

#Show Q-statistics and p-values for each SNP across all exposures
pleioX1X2X3$qdat

#Show global Q-statistics for all exposures
pleioX1X2X3$gq

#Generate unadjusted and adjusted radial MVMR plots
X1X2X3plots<-plot_rmvmr(X1X2X3f.data,X1X2X3res)
X1X2X3plots$p1
X1X2X3plots$p2



#############################################################
###   Radial MVMR Analyses: X1, X2, and X3 with pruning   ###
#############################################################
X1X2X3fp.data<-X1X2X3f.data

outliers<-"start"

while(!is.null(outliers)){
  
  X1X2X3pres<-ivw_rmvmr(X1X2X3fp.data,F)
  pleioX1X2X3p<-pleiotropy_rmvmr(X1X2X3fp.data,X1X2X3pres)
  
  ###Prune###
  
  if(min(pleioX1X2X3p$qdat$qj_p)< 0.05){
    outliers<-pleioX1X2X3p$qdat[pleioX1X2X3p$qdat$qj_p <0.05,]$snp
    outliers<-unique(outliers)
    
    X1X2X3fp.data<-X1X2X3fp.data[!X1X2X3fp.data$SNP %in% outliers,]
    
  }else{
    outliers<-NULL
    
  }
  
}

#Estimate conditional instrument strength
strength_X1X2X3p<- strength_rmvmr(X1X2X3fp.data)

#Show conditional F statistic
strength_X1X2X3p$f

#Generate conditional strength radial MVMR plot for exposure 1
strength_X1X2X3p$plot[1]

#Generate conditional strength radial MVMR plot for exposure 2
strength_X1X2X3p$plot[2]

#Generate conditional strength radial MVMR plot for exposure 3
strength_X1X2X3p$plot[3]

#Estimate causal effects  using Radial MVMR
X1X2X3resp<-ivw_rmvmr(X1X2X3fp.data,T)

#Show Q-statistics and p-values for each SNP across all exposures
pleioX1X2X3p$qdat

#Show global Q-statistics for all exposures
pleioX1X2X3p$gq

#Generate unadjusted and adjusted radial MVMR plots
X1X2X3pplots<-plot_rmvmr(X1X2X3fp.data,X1X2X3resp)
X1X2X3pplots$p1
X1X2X3pplots$p2


##Plot for detecting outliers

plotdat<-pleioX1X2X3$qdat

plotdat<-plotdat[plotdat$qj_p < 0.05/(nrow(plotdat)),]

names(plotdat)[6]<- "Exposure"

plotdat$Exposure<-factor(plotdat$Exposure)

levels(plotdat$Exposure)<-c("Exposure 1","Exposure 2","Exposure 3")

F6<-ggplot(data = plotdat, aes(x = snp, y = -log10(qj_p)))+geom_point(aes(colour = Exposure))+theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"))+ylab(expression(-log10(p)))+xlab("SNP")+ggtitle("Q-statistics for identified outliers")+
    geom_hline(yintercept = -log10(0.05/240))+
    scale_color_manual(values=c("Exposure 1"="#56B4E9","Exposure 2"="#D55E00","Exposure 3"="#009E73"))

png(filename = "Figure 6.png",
    width = 1920 , height = 1080, units = "px", res=200,
    bg = "white")

F6

dev.off()


#RMVMR Plots

png(filename = "Figure 7A.png",
    width = 1920 , height = 1920, units = "px", res=300,
    bg = "white")

X1X2plots$p1

dev.off()

png(filename = "Figure 7B.png",
    width = 1920 , height = 1920, units = "px", res=300,
    bg = "white")

X1X2plots$p2

dev.off()

png(filename = "Figure 7C.png",
    width = 1920 , height = 1920, units = "px", res=300,
    bg = "white")

X1X2X3plots$p1

dev.off()

png(filename = "Figure 7D.png",
    width = 1920 , height = 1920, units = "px", res=300,
    bg = "white")

X1X2X3plots$p2

dev.off()


png(filename = "Figure 7E.png",
    width = 1920 , height = 1920, units = "px", res=300,
    bg = "white")

X1X2X3pplots$p1

dev.off()

png(filename = "Figure 7F.png",
    width = 1920 , height = 1920, units = "px", res=300,
    bg = "white")

X1X2X3pplots$p2

dev.off()




