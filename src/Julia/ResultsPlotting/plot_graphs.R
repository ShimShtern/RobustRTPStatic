library(latex2exp)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel)
library("rstudioapi")                                 # Load rstudioapi package

setwd(dirname(getActiveDocumentContext()$path))       # Set working directory to source file location


Folder="../ResultsFiles" #Relative to this location
results="patient4_visit1_20230906_no_dose_vol_no_zeros.csv"
evaluation="patient4_visit1_20230930_solution_evaluations.csv"
PlotFolder =sprintf("%s/Plots/Patient4/",Folder)


ResultsFile=sprintf("%s/%s",Folder,results)
EvaluationFile=sprintf("%s/%s",Folder,evaluation)


snip_df<-function(df,g_limit,field){
  if(field=="test"){
    temp1=df%>%group_by(Method,gamma,delta,gamma_test,delta_test)%>%arrange(gamma,delta,mu_test,by_group = TRUE)%>%mutate(g_test_next=lag(g_test),mu_test_next=lag(mu_test))%>%ungroup()%>%
      group_by(Method,gamma,delta,gamma_test,delta_test)%>%arrange(desc(mu_test),by_group = TRUE)%>%mutate(g_test_prev=lag(g_test),mu_test_prev=lag(mu_test))%>%
      filter(g_test<g_limit)%>%mutate(mu_test_min=mu_test,g_test_min=g_test)%>%ungroup()
    temp2=temp1%>%filter(g_test_next>=g_limit)%>%mutate(mu_test=mu_test_next,g_test=g_test_next)
    temp4=temp1%>%filter(g_test_prev>=g_limit)%>%mutate(mu_test=mu_test_prev,g_test=g_test_prev)
    temp3=rbind(temp2,temp4)%>%mutate(mu_test=mu_test+(mu_test_min-mu_test)*(g_limit-g_test)/(g_test_min-g_test),g_test=g_limit)%>%select(-g_test_min,-mu_test_min,-mu_test_prev,-mu_test_next,-g_test_prev,-g_test_next)
    dfnew=rbind(df,temp3)
    dfnew=dfnew%>%group_by(Method,gamma,delta)%>%arrange(mu_test,by_group = TRUE)%>%ungroup()
  }
  else if (field=="phys"){
    temp1=df%>%group_by(Method,gamma,delta)%>%arrange(gamma,delta,phys_hom,by_group = TRUE)%>%mutate(g_phys_next=lag(min_phys_dose),mu_phys_next=lag(phys_hom))%>%ungroup()%>%
      group_by(Method,gamma,delta)%>%arrange(desc(phys_hom),by_group = TRUE)%>%mutate(g_phys_prev=lag(min_phys_dose),mu_phys_prev=lag(phys_hom))%>%
      filter(min_phys_dose<g_limit)%>%mutate(mu_phys_min=phys_hom,g_phys_min=min_phys_dose)%>%ungroup()
    temp2=temp1%>%filter(g_phys_next>=g_limit)%>%mutate(phys_hom=mu_phys_next,min_phys_dose=g_phys_next)
    temp4=temp1%>%filter(g_phys_prev>=g_limit)%>%mutate(phys_hom=mu_phys_prev,min_phys_dose=g_phys_prev)
    temp3=rbind(temp2,temp4)%>%mutate(phys_hom=phys_hom+(mu_phys_min-phys_hom)*(g_limit-min_phys_dose)/(g_phys_min-min_phys_dose),min_phys_dose=g_limit)%>%select(-g_phys_min,-mu_phys_min,-mu_phys_prev,-mu_phys_next,-g_phys_prev,-g_phys_next)
    dfnew=rbind(df,temp3)
    dfnew=dfnew%>%group_by(Method,gamma,delta)%>%arrange(phys_hom,by_group = TRUE)%>%ungroup()
    
  }
  else if (field=="EUD"){
    temp1=df%>%group_by(Method,gamma,delta)%>%arrange(gamma,delta,phys_hom,by_group = TRUE)%>%mutate(g_phys_next=lag(EUD),mu_phys_next=lag(phys_hom))%>%ungroup()%>%
      group_by(Method,gamma,delta)%>%arrange(desc(phys_hom),by_group = TRUE)%>%mutate(g_phys_prev=lag(EUD),mu_phys_prev=lag(phys_hom))%>%
      filter(EUD<g_limit)%>%mutate(mu_phys_min=phys_hom,g_phys_min=EUD)%>%ungroup()
    temp2=temp1%>%filter(g_phys_next>=g_limit)%>%mutate(phys_hom=mu_phys_next,EUD=g_phys_next)
    temp4=temp1%>%filter(g_phys_prev>=g_limit)%>%mutate(phys_hom=mu_phys_prev,EUD=g_phys_prev)
    temp3=rbind(temp2,temp4)%>%mutate(phys_hom=phys_hom+(mu_phys_min-phys_hom)*(g_limit-EUD)/(g_phys_min-EUD),EUD=g_limit)%>%select(-g_phys_min,-mu_phys_min,-mu_phys_prev,-mu_phys_next,-g_phys_prev,-g_phys_next)
    dfnew=rbind(df,temp3)
    dfnew=dfnew%>%group_by(Method,gamma,delta)%>%arrange(phys_hom,by_group = TRUE)%>%ungroup()
    
  }
  return(dfnew)
}

add_zero_dose<-function(df){
  df_temp=df%>%group_by(Method,gamma,delta,gamma_test,delta_test)%>%slice(which.min(mu))%>%ungroup()%>%mutate(mu=mu-0.0001,mu_test=1,g_test=0)
  if("phys_hom" %in% colnames(df)){
      df_temp=df_temp%>%mutate(phys_hom=1,min_phys_dose=0,EUD=0,g=0)
  }
  df=rbind(df,df_temp)
  return(df)
}

df_choose_marker<-function(df,field){
  if(field=="test"){
    df_marker=df%>%group_by(gamma,delta,Method,delta_test,gamma_test)%>%arrange(gamma,delta,mu_test)%>%mutate(mu_test_diff1=abs(mu_test-lag(mu_test,default=0)),g_test_diff1=abs(g_test-lag(g_test,default=0)))%>%
    mutate(is_mark1=(mu_test_diff1>=0.01)|(g_test_diff1>=0.5))%>%
    arrange(mu_test,order="decreasing")%>%mutate(mu_test_diff2=abs(mu_test-lag(mu_test,default=0)),g_test_diff2=abs(g_test-lag(g_test,,default=0)))%>%
    mutate(is_mark2=(mu_test_diff2>=0.01)|(g_test_diff2>=0.5))%>%
    filter(is_mark1|is_mark2)
  }
  if (field=="phys"){
    df_marker=df%>%group_by(gamma,delta,Method)%>%arrange(gamma,delta,phys_hom,.by_group=TRUE)%>%mutate(phys_hom_diff1=abs(phys_hom-lag(phys_hom,default=0)),g_phys_diff1=abs(min_phys_dose-lag(min_phys_dose,default=0)))%>%
      mutate(is_mark1=(phys_hom_diff1>=0.01)|(g_phys_diff1>=0.5))%>%
      arrange(desc(phys_hom),.by_group=TRUE)%>%mutate(phys_hom_diff2=abs(phys_hom-lag(phys_hom,default=0)),g_phys_diff2=abs(min_phys_dose-lag(min_phys_dose,,default=0)))%>%
      mutate(is_mark2=(phys_hom_diff2>=0.01)|(g_phys_diff2>=0.5))%>%
      filter(is_mark1|is_mark2)
  }
  if (field=="EUD"){
    df_marker=df%>%group_by(gamma,delta,Method)%>%arrange(gamma,delta,phys_hom,.by_group=TRUE)%>%mutate(phys_hom_diff1=abs(phys_hom-lag(phys_hom,default=0)),g_phys_diff1=abs(EUD-lag(EUD,default=0)))%>%
      mutate(is_mark1=(phys_hom_diff1>=0.01)|(g_phys_diff1>=0.5))%>%
      arrange(desc(phys_hom),.by_group=TRUE)%>%mutate(phys_hom_diff2=abs(phys_hom-lag(phys_hom,default=0)),g_phys_diff2=abs(EUD-lag(EUD,,default=0)))%>%
      mutate(is_mark2=(phys_hom_diff2>=0.01)|(g_phys_diff2>=0.5))%>%
      filter(is_mark1|is_mark2)
  }
  return(df_marker)
}


#create basic data set
df1 <- read.csv(ResultsFile)
df1 <- subset(df1,!(mu %in% c(1.05,1.135)))%>%mutate(delta=round(delta,2),gamma=round(gamma,2))
df1 <- subset(df1,(g != 0))
df1 <- subset(df1, select = c(mu,delta,gamma,beta, lambda1, lambda2,t1,t2,t3,g,reg1,reg2,phys_hom,min_phys_dose,time,delta_test,gamma_test,mu_test,g_test,EUD) )
dd=df1%>%group_by(gamma,delta)%>%summarize(min_mutest=min(mu_test),min_gtest=min(g_test),minmutest_gtest=g_test[which.min(mu_test)],minmutest_mu=mu[which.min(mu_test)],mingtest_mutest=mu_test[which.min(g_test)],mingtest_mu=mu[which.min(g_test)])
df1$g_test <- as.numeric(df1$g_test)
head(df1)

df_box <- subset(df1, gamma==1)
df_box <- subset(df_box, delta>0)
df_box <- subset(df_box,g_test>0)
df_box$mu_test <- as.numeric(df_box$mu_test)
df1a <- df1%>%filter(delta>0, gamma<1)
df1a <- subset(df1a,!is.na(g_test))
df1a$Method=sprintf("Robust $\\gamma=%.2f$",df1a$gamma)
dfnom <- df1%>%filter(delta == 0 , gamma == 1)
for(delta in seq(0.01,0.14,0.01)){
  newdata<-data.frame(dfnom)
  newdata$delta=delta
  newdata$Method=("Nominal")
  df1a=rbind(df1a,newdata)
}
df_box$Method=("Robust Box")
df1a=rbind(df1a,df_box)
df1a=df1a%>%filter(round(delta*100)%in%c(seq(2,14,2)))
df1a$gamma=as.numeric(df1a$gamma)
df1a$delta=round(as.numeric(df1a$delta),2)
df1a=df1a[order(df1a$Method),]
Methods=unique(df1a$Method)
MethodLabels=lapply(Methods,TeX)
df1a$deltaLabel=
  factor(df1a$delta, labels=sprintf("delta~'='~%.2f",as.numeric(unique(sort(df1a$delta)))))


##################################
#figure 13 start here
gamma_val=0.04
df2 = subset(df1,gamma %in% c(gamma_val) & mu<1.18)
df3 <- subset(df2,!is.na(g_test))
df3 <-subset(df3,(mu*10000-875)%%125==0)
g <- ggplot(data=df3,aes(x=delta)) + geom_point(aes(group=mu,y=g,color=as.factor(mu))) + geom_line(aes(group=mu,y=g,color=as.factor(mu),linetype=as.factor(mu))) + ylab(TeX("\\underline{d}")) + xlab(TeX("$\\delta")) #+ xlim(0.01,0.11) + ylim(72.5,80.0) #+ ylab('g-test') + xlab(TeX("$\\delta"))
#g <- g + geom_point(aes(group=mu,y=g_test,color=as.factor(mu)),shape=22) + geom_line(aes(group=mu,y=g_test,color=as.factor(mu),linetype=as.factor(mu))) 
g +  labs(linetype=TeX("$\\mu$"),color=TeX("$\\mu$"))+ 
  scale_x_continuous(breaks=c(seq(0,0.14,0.02)))+
  scale_y_continuous(breaks=c(seq(45,60,1)))+
  theme_bw(base_size = 12, base_family = "Helvetica")  +
  theme(legend.title = element_text(size = 10), 
       legend.text  = element_text(size = 10),
       legend.background = element_rect(color="black"),
       legend.direction="vertical", 
       strip.background = element_blank())
outputFile=sprintf("./%s/%s%.2f%s",PlotFolder,"robustGvsDelta_gamma",gamma_val, ".png")
ggsave(outputFile,height=3,width=5)
#figure 13 ends here



df1evaluate <- read.csv(EvaluationFile)
df1evaluate <- df1evaluate%>%filter(delta_test==0,gamma_test==1)

dfnom <- df1evaluate%>%filter(delta == 0,gamma==1)
dfnom$Method="Nominal"
#dfnom=add_zero_dose(dfnom)
dfnom = rbind(dfnom,dfnom%>%mutate(gamma=0.04),dfnom%>%mutate(gamma=0.02))
df5 <- df1evaluate%>%arrange(gamma,delta,mu)%>%filter(gamma %in% c(1)) %>% filter(round(delta*100)%%2==0, delta>0)
df5$Method=sprintf("Robust $\\delta=%.2f$",df5$delta)
#df5=add_zero_dose(df5)

## nominal performance for different gamma values, showing diffrent delta values
#figure 6 starts here
gtest_limit=25
dfnomnew=snip_df(dfnom,gtest_limit,'test')
df5new=snip_df(df5,gtest_limit,'test') 

df5marks=df_choose_marker(df5,"test")
dfnommarks=df_choose_marker(dfnom,"test")
Methods="Nominal"
Methods=c(Methods,unique(df5$Method))
MethodLabels=lapply(Methods,TeX)
g <- ggplot(data=df5new) 
#df5marks
g <- g + geom_point(data=df5marks, aes(x=mu_test,y=g_test,color=Method,shape=Method)) + geom_line(aes(x=mu_test,y=g_test,color=Method,linetype=Method)) 
#dfnommarks
g <- g + geom_point(data=dfnommarks,aes(x=mu_test,y=g_test,color=Method,shape=Method)) + geom_line(data=dfnomnew,aes(x=mu_test,y=g_test,color=Method,linetype=Method)) 
g <- g + scale_shape_manual(label=MethodLabels,values=c(16, 17, 15,18,6,3,4,2,8,5,7))+
  scale_color_discrete(label=MethodLabels )+
  #scale_size_manual(values=c(2,2,2,3,2,3,3,2), label=MethodLabels, guide=guide_legend(ncol=legendcolnum))+
  scale_linetype_discrete(label=MethodLabels)+
  scale_x_continuous(breaks=seq(1.03,1.2,0.02),lim=c(1.07,1.15))+
  scale_y_continuous(lim=c(gtest_limit,57.5))+
  ylab(TeX("\\hat{\\underline{d}}")) + 
  xlab(TeX("$\\hat{\\mu}")) + 
  labs(linetype="Method",color="Method",shape="Method") + 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 8),
        legend.background = element_rect(color="black"),
        #legend.spacing = unit(-1, "cm"),
        #legend.key.height= unit(0.05, "npc"),
        #legend.margin = margin(),
        legend.direction="vertical", 
        strip.background = element_blank())
g
outputFile=sprintf("./%s/%s",PlotFolder,"Patient4_Visit1_gamma1robloss.png")
ggsave(outputFile,width=3.9,height=2.7,units=c("in"))

for (gamma_val in c(0.02,0.04)){
  df4 <- subset(df1evaluate,gamma %in% c(gamma_val) &  round(delta*100)%%2==0 & delta>0)
  df4 <- subset(df4,!is.na(g_test))
  df4$Method=sprintf("Robust $\\delta=%.2f$",df4$delta)
  #df4=add_zero_dose(df4)
  df4new=snip_df(df4,gtest_limit,"test")
  df4marks=df_choose_marker(df4,"test")
  g <- ggplot(data=df4new) #+ geom_point(aes(group=delta,y=g_g_geq_0,color=as.factor(delta))) + geom_line(aes(group=delta,y=g_g_geq_0,color=as.factor(delta),linetype=as.factor(delta))) + ylab(TeX("\\underline{d}")) + xlab(TeX("$\\mu")) #+ xlim(0.01,0.11) + ylim(72.5,80.0) #+ ylab('g-test') + xlab(TeX("$\\delta"))
  g <- g + geom_point(data=df4marks,aes(group=delta,x=mu_test,y=g_test,color=Method,shape=Method)) + geom_line(aes(group=delta,x=mu_test,y=g_test,color=Method,linetype=Method))
  g <- g + geom_point(data=dfnommarks,aes(x=mu_test,y=g_test,color=Method,shape=Method)) + geom_line(data=dfnomnew,aes(x=mu_test,y=g_test,color=Method,linetype=Method)) 
  g <- g +   scale_x_continuous(breaks=seq(1.03,1.18,0.02),lim=c(1.07,1.15))+
    scale_y_continuous(lim=c(gtest_limit,57.5))+
    scale_shape_manual(label=MethodLabels,values=c(16, 17, 15,18,6,3,4,2,8,5,7))+
    scale_color_discrete(label=MethodLabels )+
    #scale_size_manual(values=c(2,2,2,3,2,3,3,2), label=MethodLabels, guide=guide_legend(ncol=legendcolnum))+
    scale_linetype_discrete(label=MethodLabels)+
    ylab(TeX("\\hat{\\underline{d}}")) + xlab(TeX("$\\hat{\\mu}")) +
    labs(linetype=TeX("$\\delta$"),color=TeX("$\\delta$")) + 
    theme_bw(base_size = 12, base_family = "Helvetica") + #ylim(81,83.3)  + 
    theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.position="none",
        legend.background = element_rect(color="black"),
        legend.direction="vertical", 
        strip.background = element_blank())
  g
  outputFile=sprintf("./%s/%s%.2f%s",PlotFolder,"Patient4_Visit1_gamma",gamma_val,"robloss.png")
  ggsave(outputFile,width=2.5,height=2.7,units=c("in"))
}
#figure 6 ends here


## plot nominal evaluation of nominal vs robust for delta=0.06 and different gamma

df_box <- subset(df1evaluate, gamma==1)
df_box <- subset(df_box, delta>0)
df_box <- subset(df_box,g_test>0)
df_box$mu_test <- as.numeric(df_box$mu_test)

df6 <- df1evaluate%>%filter(delta>0, gamma<1)
df6 <- subset(df6,!is.na(g_test))
df6$Method=sprintf("Robust $\\gamma=%.2f$",df6$gamma)
dfnom <- df1evaluate%>%filter(delta == 0 , gamma == 1)
for(delta in seq(0.01,0.14,0.01)){
  newdata<-data.frame(dfnom)
  newdata$delta=delta
  newdata$Method=("Nominal")
  df6=rbind(df6,newdata)
}
df_box$Method=("Robust Box")
df6=rbind(df6,df_box)
df6=df6%>%filter(round(delta*100)%in%c(seq(2,14,2)))
df6$gamma=as.numeric(df6$gamma)
df6$delta=round(as.numeric(df6$delta),2)
df6=df6[order(df6$Method),]
Methods=unique(df6$Method)
MethodLabels=lapply(Methods,TeX)
df6$deltaLabel=
  factor(df6$delta, labels=sprintf("delta~'='~%.2f",as.numeric(unique(sort(df6$delta)))))

df6_marks=df_choose_marker(df6,"test")
snip_g_test_val=47.5
df6_new=snip_df(df6,snip_g_test_val,"test")

#plot 14 starts here
legendcolnum=2
facetcolnum=3
g <- ggplot() 
g <- g + geom_point(data=df6_marks,aes(x=mu_test,y=g_test,color=Method,shape=Method,size=Method)) + geom_line(data=df6_new,aes(x=mu_test,y=g_test,color=Method,linetype=Method)) 
g +   facet_wrap(~deltaLabel,ncol = facetcolnum,labeller=label_parsed)+
  scale_shape_manual(label=MethodLabels, values=c(16, 17, 15,18,6,3,4,0),guide=guide_legend(ncol=legendcolnum) )+
  scale_color_discrete(label=MethodLabels,  guide=guide_legend(ncol=legendcolnum) )+
  scale_size_manual(label=MethodLabels, values=c(2,2,2,3,2,3,3,2),  guide=guide_legend(ncol=legendcolnum))+
  scale_linetype_discrete(label=MethodLabels, guide=guide_legend(ncol=legendcolnum) )+
  scale_x_continuous(breaks=c(1.075,1.1,1.125,1.15,1.175,1.2,1.225,1.25,1.275,1.3),lim=c(1.06,1.15))+
  ylab(TeX("\\hat{\\underline{d}}")) + xlab(TeX("$\\hat{\\mu}$")) + 
  labs(linetype="Method",color="Method",shape="Method") + 
  theme_bw(base_size = 12, base_family = "Helvetica") + ylim(snip_g_test_val,57.5) +
  theme(legend.title = element_text(size = 10), legend.text  = element_text(size = 10),
        #legend.position = c(.1, 0.05),
        #legend.position = c(.35, .2),
        legend.position = c(.4, .25),
        legend.justification = c(0, 1),
        legend.box = "vertical",
        legend.background = element_rect(color="black"),
        legend.direction="vertical", legend.margin=margin(t = 1, r = 6, b = 1, l = 3, unit = "pt"),
        strip.background = element_blank())
outputFile=sprintf("./%s/%s",PlotFolder,"Patient4_Visit1_nominal_performance.png")
ggsave(outputFile,height=7,width=6)
dev.off()
#plot 14 ends here

#figure 5 starts here
delta_val=0.14
df6b=df6%>%filter(delta==delta_val)
df6b_marks=df6_marks%>%filter(delta==delta_val)
g_limit=35
df6b_new=snip_df(df6b,g_limit,"test")
colnumlegend=1
g <- ggplot(data=df6b_new) 
g <- g + geom_point(data=df6b_marks,aes(x=mu_test,y=g_test,color=Method,shape=Method,size=Method)) + geom_line(aes(x=mu_test,y=g_test,color=Method,linetype=Method)) 
g +  scale_shape_manual(label=MethodLabels, values=c(16, 17, 15,18,6,3,4,0),guide=guide_legend(ncol=colnumlegend) )+
  scale_color_discrete(label=MethodLabels,guide=guide_legend(ncol=colnumlegend) )+
  scale_size_manual(label=MethodLabels, values=c(2,2,2,3,2,3,3,2), guide=guide_legend(ncol=colnumlegend))+
  scale_linetype_discrete(label=MethodLabels, guide=guide_legend(ncol=colnumlegend) )+
  scale_x_continuous(breaks=seq(1.025,1.5,0.025),lim=c(1.075,1.15))+
  ylab(TeX("\\hat{\\underline{d}}")) + xlab(TeX("$\\hat{\\mu}$")) + 
  labs(linetype="Method",color="Method",shape="Method") + 
  theme_bw(base_size = 12, base_family = "Helvetica") + ylim(g_limit,57.5) +
  theme(legend.title = element_text(size = 10), legend.text  = element_text(size = 10),
        legend.box = "vertical",
        legend.background = element_rect(color="black"),
        legend.direction="vertical", legend.margin=margin(t = 1, r = 6, b = 1, l = 3, unit = "pt"),
        strip.background = element_blank())
outputFile=sprintf("./%s/%s%.2f%s",PlotFolder,"Patient4_Visit1_nominal_performance_delta", delta_val,".png")
ggsave(outputFile,height=2.8,width=6)
dev.off()
#figure 5 ends here


## plot physical performance

#start figure 15
min_phys_dose_value=50
df1a$deltaLabel=
  factor(df1a$delta, labels=sprintf("delta~'='~%.2f",as.numeric(unique(sort(df6_new$delta)))))
df6_new=snip_df(df1a,min_phys_dose_value,"phys")
df6_marks=df_choose_marker(df1a,"phys")
legendcolnum=2
facetcolnum=3
g <- ggplot(data=df6_new) 
g <- g + geom_point(data=df6_marks,aes(x=phys_hom,y=min_phys_dose,color=Method,shape=Method,size=Method)) + geom_line(aes(x=phys_hom,y=min_phys_dose,color=Method,linetype=Method)) 
g<-g +   facet_wrap(~deltaLabel,ncol = facetcolnum,labeller=label_parsed)+
  scale_shape_manual(label=MethodLabels, values=c(16, 17, 15,18,6,3,4,0),guide=guide_legend(ncol=legendcolnum) )+
  scale_color_discrete(label=MethodLabels, guide=guide_legend(ncol=legendcolnum) )+
  scale_size_manual(values=c(2,2,2,3,2,3,3,2), label=MethodLabels, guide=guide_legend(ncol=legendcolnum))+
  scale_linetype_discrete(label=MethodLabels, guide=guide_legend(ncol=legendcolnum) )+
  scale_y_continuous(lim=c(min_phys_dose_value,63.5))+
  scale_x_continuous(breaks=c(seq(1,1.18,0.05)))+
  ylab(TeX("Min physical dose")) + xlab(TeX("Physical $\\mu$")) + 
  labs(linetype="Method",color="Method",shape="Method") + 
  theme_bw(base_size = 12, base_family = "Helvetica") +  
  theme(legend.title = element_text(size = 10), legend.text  = element_text(size = 10),
        legend.position = c(.4, .25),
        legend.justification = c(0, 1),
        legend.box = "vertical",
        legend.background = element_rect(color="black"),
        legend.direction="vertical", legend.margin=margin(t = 1, r = 6, b = 1, l = 3, unit = "pt"),
        strip.background = element_blank())
outputFile=sprintf("./%s/%s",PlotFolder,"Patient4_Visit1_physical_performance.png")
ggsave(outputFile,height=7,width=6)
dev.off()
#end figure 15

#start figure 7A
delta_val=0.14
df6b=df1a%>%filter(delta==delta_val)
min_phys_dose_value=35
df6b=snip_df(df6b,min_phys_dose_value,"phys")
df6b_marks=df6_marks%>%filter(delta==delta_val)
legendcolnum=1
g <- ggplot(data=df6b) 
g <- g + geom_point(data=df6b_marks,aes(x=phys_hom,y=min_phys_dose,color=Method,shape=Method,size=Method)) + geom_line(aes(x=phys_hom,y=min_phys_dose,color=Method,linetype=Method)) 
g +   scale_shape_manual(label=MethodLabels, values=c(16, 17, 15,18,6,3,4,0),guide=guide_legend(ncol=legendcolnum) )+
  scale_color_discrete(label=MethodLabels, guide=guide_legend(ncol=legendcolnum) )+
  scale_size_manual(values=c(2,2,2,3,2,3,3,2), label=MethodLabels, guide=guide_legend(ncol=legendcolnum))+
  scale_linetype_discrete(label=MethodLabels, guide=guide_legend(ncol=legendcolnum) )+
  scale_x_continuous(breaks=c(seq(1.01,1.18,0.02)),lim=c(1.01,1.15))+
  scale_y_continuous(lim=c(min_phys_dose_value,64))+
  ylab(TeX("Min physical dose")) + xlab(TeX("Physical $\\mu$")) + 
  labs(linetype="Method",color="Method",shape="Method") + 
  theme_bw(base_size = 12, base_family = "Helvetica") + #ylim(90,94) +  
  theme(legend.position = "none",
        legend.title = element_text(size = 10), legend.text  = element_text(size = 10),
        legend.box = "vertical",
        legend.background = element_rect(color="black"),
        legend.direction="vertical", legend.margin=margin(t = 1, r = 6, b = 1, l = 3, unit = "pt"),
        strip.background = element_blank())
outputFile=sprintf("./%s/%s%.2f%s",PlotFolder,"Patient4_Visit1_physical_performance_delta",delta_val,".png")
ggsave(outputFile,height=2.8,width=4.5)
dev.off()
#end figure 7A

##plot EUD

#not figure in paper
df9=df1a
min_EUD=40
df9_new=snip_df(df1a,min_EUD,"EUD")
df9_marks=df_choose_marker(df9,"EUD")
legendcolnum=2
facetcolnum=3
g <- ggplot(data=df9_new) 
g <- g + geom_point(data=df9_marks,aes(x=phys_hom,y=EUD,color=Method,shape=Method,size=Method)) + geom_line(aes(x=phys_hom,y=EUD,color=Method,linetype=Method)) 
g +   facet_wrap(~deltaLabel,ncol = facetcolnum,labeller=label_parsed)+
  scale_shape_manual(label=MethodLabels, values=c(16, 17, 15,18,6,3,4,0),guide=guide_legend(ncol=legendcolnum) )+
  scale_color_discrete(label=MethodLabels, guide=guide_legend(ncol=legendcolnum) )+
  scale_size_manual(values=c(2,2,2,3,2,3,3,2), label=MethodLabels, guide=guide_legend(ncol=legendcolnum))+
  scale_linetype_discrete(label=MethodLabels, guide=guide_legend(ncol=legendcolnum) )+
  scale_y_continuous(lim=c(min_EUD,66))+
  scale_x_continuous(breaks=c(seq(1,1.18,0.02)))+
  ylab(TeX("EUD")) + xlab(TeX("Physical $\\mu$")) + 
  labs(linetype="Method",color="Method",shape="Method") + 
  theme_bw(base_size = 12, base_family = "Helvetica") +  
  theme(legend.title = element_text(size = 10), legend.text  = element_text(size = 10),
        legend.position = c(.4, .25),
        legend.justification = c(0, 1),
        legend.box = "vertical",
        legend.background = element_rect(color="black"),
        legend.direction="vertical", legend.margin=margin(t = 1, r = 6, b = 1, l = 3, unit = "pt"),
        strip.background = element_blank())
outputFile=sprintf("./%s/%s",PlotFolder,"Patient4_Visit1_EUD_performance.png")
ggsave(outputFile,height=7,width=6)
dev.off()

#start figure 7B
delta_val=0.14
df9b=df9%>%filter(delta==delta_val)
min_EUD=40
df9b=snip_df(df9b,min_EUD,"EUD")
df9b_marks=df9_marks%>%filter(delta==delta_val)
legendcolnum=1
g <- ggplot(data=df9b) 
g <- g + geom_point(data=df9b_marks,aes(x=phys_hom,y=EUD,color=Method,shape=Method,size=Method)) + geom_line(aes(x=phys_hom,y=EUD,color=Method,linetype=Method)) 
g +   scale_shape_manual(label=MethodLabels, values=c(16, 17, 15,18,6,3,4,0),guide=guide_legend(ncol=legendcolnum) )+
  scale_color_discrete(label=MethodLabels, guide=guide_legend(ncol=legendcolnum) )+
  scale_size_manual(values=c(2,2,2,3,2,3,3,2), label=MethodLabels, guide=guide_legend(ncol=legendcolnum))+
  scale_linetype_discrete(label=MethodLabels, guide=guide_legend(ncol=legendcolnum) )+
  scale_x_continuous(breaks=c(seq(1.01,1.18,0.02)),lim=c(1.01,1.13))+
  scale_y_continuous(lim=c(min_EUD,65))+
  ylab(TeX("EUD")) + xlab(TeX("Physical $\\mu$")) + 
  labs(linetype="Method",color="Method",shape="Method") + 
  theme_bw(base_size = 12, base_family = "Helvetica") + #ylim(90,94) +  
  theme(legend.title = element_text(size = 10), legend.text  = element_text(size = 10),
        legend.box = "vertical",
        legend.background = element_rect(color="black"),
        legend.direction="vertical", legend.margin=margin(t = 1, r = 6, b = 1, l = 3, unit = "pt"),
        strip.background = element_blank())
outputFile=sprintf("./%s/%s%.2f%s",Folder,"Patient4_Visit1_EUD_delta",delta_val,".png")
ggsave(outputFile,height=2.8,width=6)
dev.off()
#end figure 7B


## sensitivity to misspesified gamma and delta values

#start figures 10 and 16-17
for (gamma_val in c(0.02,0.04)){
df7 <- read.csv(EvaluationFile)
df7=df7[df7$gamma_test == gamma_val,]
df7=df7[df7$delta_test >0,]
df7_nominal=df7%>%filter(delta==0)
df7_nominal$Method=("Nominal")
df7=df7[df7$delta == df7$delta_test,]
df7_box=df7[df7$gamma == 1,]
df7_box$Method=("Robust Box")
df7=df7[df7$gamma < 1,]
df7$Method=sprintf("Robust $\\gamma=%.2f$",df7$gamma)
df7=rbind(df7_nominal,df7,df7_box)
df7$Method=ordered(df7$Method)
Methods=levels(df7$Method)
MethodLabels=lapply(Methods,TeX)
df7$deltaLabel=
  factor(df7$delta_test, labels=sprintf("delta~'='~%.2f",as.numeric(unique(sort(df7$delta_test)))))
df7=df7%>%filter(round(delta_test*100)%%2==0)


df7a=df7%>%filter(gamma==1 | delta>0)

min_g_test=40
df7_new=snip_df(df7a,min_g_test,"test")
df7_new1=df7_new%>%filter(g_test<min_g_test)%>%mutate(g_test=NaN)
df7_new=df7_new%>%filter(g_test>=min_g_test)
df7_new=rbind(df7_new,df7_new1)%>%filter(mu_test>=1.05)
df7_marks=df_choose_marker(df7a,"test")
df7_marks=df7_marks%>%filter(g_test>=min_g_test)


legendcolnum=2
facetcolnum=3
g <- ggplot(data=df7_new) 
g <- g + geom_point(data=df7_marks,aes(x=mu_test,y=g_test,color=Method,shape=Method,size=Method)) + geom_line(aes(x=mu_test,y=g_test,color=Method,linetype=Method)) 
g +   facet_wrap(~deltaLabel,ncol = facetcolnum,labeller=label_parsed,scales="free")+
  scale_y_continuous(breaks=seq(min_g_test,57,2))+
  scale_x_continuous(breaks=seq(1.05,1.4,0.025))+
  scale_shape_manual(label=MethodLabels, values=c(16, 17, 15,18,6,3,4,0),guide=guide_legend(ncol=legendcolnum) )+
  scale_color_discrete(label=MethodLabels, guide=guide_legend(ncol=legendcolnum) )+
  scale_size_manual(values=c(2,2,2,3,2,3,3,2), label=MethodLabels, guide=guide_legend(ncol=legendcolnum))+
  scale_linetype_discrete(label=MethodLabels, guide=guide_legend(ncol=legendcolnum) )+
  ylab(TeX("\\hat{\\underline{d}}")) + xlab(TeX("$\\hat{\\mu}$")) + 
  labs(linetype="Method",color="Method",shape="Method") + 
  theme_bw(base_size = 12, base_family = "Helvetica") + #ylim(73,83) +  
  theme(legend.title = element_text(size = 10), legend.text  = element_text(size = 10),
        legend.position = c(.4, .25),#
        legend.justification = c(0.1, 1.1),
#        legend.position = "bottom",
        legend.box = "vertical",
        legend.background = element_rect(color="black"),
        #legend.direction="horizontal", legend.margin=margin(t = 1, r = 6, b = 1, l = 3, unit = "pt"),
        legend.direction="vertical", legend.margin=margin(t = 1, r = 6, b = 1, l = 3, unit = "pt"),
        strip.background = element_blank(),
        axis.text.x=element_text(angle = 45,vjust = 1.2, hjust=1)
        )
fig_name=sprintf("./%s/Patient4_Visit1_performance_gamma%.2f.png",PlotFolder,gamma_val)
ggsave(fig_name,height=8,width=7)

## fixing delta=delta test and checking sensitivity to misspesified gamma

g_limit=43
delta_test_value=0.1
df7b=df7a%>%filter((delta==delta_test | delta==0) & (delta_test==delta_test_value))
df7b=snip_df(df7b,g_limit,"test")
df7b_marks=df_choose_marker(df7a,"test")%>%filter(g_test>=g_limit & (delta==delta_test | delta==0)  & delta_test==delta_test_value)

colnumlegend=1
g <- ggplot(data=df7b) 
g <- g + geom_line(aes(x=mu_test,y=g_test,color=Method,linetype=Method)) +
  geom_point(data=df7b_marks,aes(x=mu_test,y=g_test,color=Method,shape=Method,size=Method)) 
g<- g + scale_shape_manual(label=MethodLabels, values=c(16,17,15,18,6,3,4,0),guide=guide_legend(ncol=colnumlegend) )+
  scale_color_discrete(label=MethodLabels,guide=guide_legend(ncol=colnumlegend) )+
  scale_size_manual(label=MethodLabels, values=c(2,2,2,3,2,3,3,2), guide=guide_legend(ncol=colnumlegend))+
  scale_linetype_discrete(label=MethodLabels, guide=guide_legend(ncol=colnumlegend) )+
  scale_y_continuous(breaks=c(seq(20,60,2)),lim=c(g_limit,52))+
  ylab(TeX("\\hat{\\underline{d}}")) + xlab(TeX("$\\hat{\\mu}$")) + 
  labs(linetype="Method",color="Method",shape="Method") + 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.box = "vertical",
        legend.background = element_rect(color="black"),
        legend.direction="vertical", legend.margin=margin(t = 1, r = 6, b = 1, l = 3, unit = "pt"))
fig_name=sprintf("%s/Patient4_Visit1_performance_gamma%.2f_delta%.2f.png",PlotFolder,gamma_val,delta_test_value)
if (gamma_val==0.02){
  g+scale_x_continuous(breaks=seq(1.01,1.4,0.02),lim=c(1.09,1.23))+
    theme(legend.position="none")
  ggsave(fig_name,height=3,width=3.5)
}
else if (gamma_val==0.04){
  g+scale_x_continuous(breaks=seq(1.02,1.4,0.02),lim=c(1.12,1.26))+
    theme(legend.title = element_text(size = 10), legend.text  = element_text(size = 10))
  ggsave(fig_name,height=3,width=5)
}

}
#end figures 10, and 16-17


## sensitivity to misspesified delta values

df8 <- read.csv(EvaluationFile)
df8 <-  subset(df8, round(delta*100)%%2 ==0 & gamma==gamma_test & delta_test%in%c(0.02,0.04,0.06,0.08,0.1,0.12,0.14) & gamma_test %in% c(0.02,0.04))

df8 <-df8[order(df8$delta),]
df8$Method=sprintf("Robust $\\delta=%.2f$",df8$delta)
df8$Method[df8$delta==0]=sprintf("Nominal")
df8$gammaLabel=
  factor(df8$gamma_test, labels=sprintf("gamma~'='~%.2f",as.numeric(unique(sort(df8$gamma_test)))))
df8$deltaLabel=
  factor(df8$delta_test, labels=sprintf("delta~'='~%.2f",as.numeric(unique(sort(df8$delta_test)))))

#fix delta and gamma values and plot performance of misspecified delta

#start figure 11
gamma_val=0.04
for (delta_val in c(0.02,0.14)){
  if (delta_val==0.02){
    g_limit=50
    g_ub =56
    x_lb = 1.1
    x_ub = 1.24
  } else if(delta_val==0.14){
    g_limit=43.5
    g_ub = 48.5
    x_lb = 1.12
    x_ub = 1.26
  }
  df8b=df8[df8$delta_test==delta_val & df8$gamma_test==gamma_val,]
  df8b_new=snip_df(df8b,g_limit,"test")
  df8b_new$Method=ordered(df8b_new$Method)
  MethodLabels=lapply(levels(df8b_new$Method),TeX)
  df8b_marker=df_choose_marker(df8b,"test")%>%filter(g_test>=g_limit)
  colnumlegend=1
  g <- ggplot(data=df8b_new) 
  g <- g + geom_point(data=df8b_marker,aes(x=mu_test,y=g_test,color=Method,shape=Method,size=Method)) + geom_line(aes(x=mu_test,y=g_test,color=Method,linetype=Method)) 
  g<-  g +  scale_x_continuous(breaks=seq(1.1,1.4,0.02),lim=c(x_lb,x_ub))+
    scale_y_continuous(breaks=c(seq(40,60,1)),lim=c(g_limit,g_ub))+
    scale_shape_manual(label=MethodLabels,  values=c(16, 17, 15,18,6,3,4,2,8,5,7,9,0),guide=guide_legend(ncol=colnumlegend) )+
    scale_color_discrete(label=MethodLabels,  guide=guide_legend(ncol=colnumlegend) )+
    scale_size_manual(values=c(2,2,2,3,2,3,3,2,2,2,2), label=MethodLabels, guide=guide_legend(ncol=colnumlegend))+
    scale_linetype_discrete(label=MethodLabels, guide=guide_legend(ncol=colnumlegend) )+
    ylab(TeX("\\hat{\\underline{d}}")) + xlab(TeX("$\\hat{\\mu}$")) +  
    labs(linetype="Method",color="Method",shape="Method") + 
    theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(legend.title = element_text(size = 10), 
      legend.text  = element_text(size = 10),
      legend.box = "vertical",
      legend.background = element_rect(color="black"),
      legend.direction="vertical", legend.margin=margin(t = 1, r = 6, b = 1, l = 3, unit = "pt"),
      strip.background = element_blank())
  fig_name=sprintf("%s/Patient4_Visit1_performance_delta%.2f_gamma%.2f.png",PlotFolder,delta_val,gamma_val)
  if (delta_val==0.02){
    g+theme(legend.position="none")
    ggsave(fig_name,height=3,width=3.5)
  } else if (delta_val==0.14){
    g+theme(legend.title = element_text(size = 10), legend.text  = element_text(size = 10))
    ggsave(fig_name,height=3,width=5)
  }
}
#end figure 11
