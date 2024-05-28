library(latex2exp)
library(ggplot2)
library(dplyr)
library(tidyr)
library("rstudioapi")                                 # Load rstudioapi package


setwd(dirname(getActiveDocumentContext()$path))       # Set working directory to source file location

Brain_voxels=1360182
Folder="../ResultsFiles" #Relative to this location
Nominal_Results = "Patient4_Visit1_ChangingBeta_delta0.00_gamma1.00_mu1.125.csv"
Robust_Results = "Patient4_Visit1_ChangingBeta_delta0.08_gamma0.04_mu1.1875.csv"
PlotFolder =sprintf("%s/Plots/Patient4/",Folder)


df1 <- read.csv(sprintf("%s/%s",Folder,Nominal_Results))
df1a <- read.csv(sprintf("%s/%s",Folder,Robust_Results))
df1=df1%>%mutate(method="Nominal")
df1a=df1a%>%mutate(method="Robust")
df1=rbind(df1,df1a)
df1=df1%>%arrange(beta)%>%mutate(beta=beta/(2*Brain_voxels)) #why 2?
df2=df1%>%filter(l0<Brain_voxels*0.01)%>%group_by(method)%>%slice(which.min(beta))%>%mutate(min_beta=beta)%>%ungroup()

##
ggplot(df1,aes(x=beta+1))+
    geom_line(aes(y=g,linetype=method),color='blue',size=1)+
    geom_line(aes(y=(l0/Brain_voxels)*100+55,linetype=method),color='red',linewidth=1)+
    scale_x_log10(TeX('\\beta+1'),breaks=c(1e-6,1e-5,1e-4,0.001,0.01,0.1,1,10,100,1000,10000,10000))+
scale_y_continuous(
  # Features of the first axis
  name = "g (Gy)",
  limits = c(54,60),
  # Add a second axis and specify its features
  sec.axis = sec_axis(trans=~(.-50),name="% deviation")
) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.title.y = element_text(color = "blue", size=13),
    axis.title.y.right = element_text(color = "red", size=13)
  )+
  geom_vline(xintercept =df2$min_beta+1,linetype="dotted",size=1)+
  annotate("text", x = df2$min_beta+1, y = 55, label = TeX(paste("\\beta=",sprintf("%2.2e",df2$min_beta),sep="")))
#ggsave("violation_and_g_vs_beta_robust.png",scale=1)
#dev.off()

## plot the comparison between nominal and robust
#start figure 12
df1$beta[df1$beta==0]=3/(2*Brain_voxels)
ggplot(df1,aes(x=beta))+
  geom_line(aes(y=g,linetype=method),color='blue',linewidth=1)+
  #geom_line(aes(y=(dev_num/Brain_voxels)*100+80),color='red',linetype="dashed",size=1)+
  scale_x_log10(TeX('\\beta'),
                breaks=c(10^seq(-6,3,1)),
                limits=c(1e-6,3e-2))+
  scale_y_continuous(
    # Features of the first axis
    name = TeX('\\underline{d} (Gy)'),
    limits = c(48,58),
    # Add a second axis and specify its features
    #sec.axis = sec_axis(trans=~(.-80),name="% deviation")
  ) +
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 20),
    legend.title = element_blank(),#element_text(size = 13), 
    legend.text  = element_text(size = 13),
    legend.background = element_rect(color="black"),
    legend.position=c(0.84, 0.9),
    legend.margin = margin(-7, 6, 5, 6),
    legend.direction="vertical"
  #  axis.title.y = element_text(color = "blue", size=13),
  #  axis.title.y.right = element_text(color = "red", size=13)
  )+
  annotate("segment", x =  df2[df2$method=="Nominal",]$min_beta, xend = df2[df2$method=="Nominal",]$min_beta, y = df2[df2$method=="Nominal",]$g-1, yend = df2[df2$method=="Nominal",]$g,
           size=0.75,linetype="solid",arrow = arrow(type = "closed",length = unit(.2,"cm")))+
  annotate("text", x = df2[df2$method=="Nominal",]$min_beta, y = df2[df2$method=="Nominal",]$g-1.5, size=5, label = TeX(paste("\\beta=",sprintf("%2.2e",df2[df2$method=="Nominal",]$min_beta),sep="")))+
  annotate("segment", x =  df2[df2$method=="Robust",]$min_beta, xend = df2[df2$method=="Robust",]$min_beta, y = df2[df2$method=="Robust",]$g+1, yend = df2[df2$method=="Robust",]$g,
           size=0.75,linetype="dashed",arrow = arrow(type = "closed",length = unit(.2,"cm")))+
  annotate("text", x = df2[df2$method=="Robust",]$min_beta, vjust=-0.1,y = df2[df2$method=="Robust",]$g+1, size=5, color='black',label = TeX(paste("\\beta=",sprintf("%2.2e",df2[df2$method=="Robust",]$min_beta),sep="")))
outputfile=sprintf("./%s/%s",PlotFolder,"Patient4_Visit1_udbar_vs_beta_nominal_robust.png")
ggsave(outputfile,height=4.2,width=5)
dev.off()


ggplot(df1,aes(x=beta))+
  geom_line(aes(y=l1,linetype=method),color='blue',size=1)+
  geom_line(aes(y=l0,linetype=method),color='red',size=1)+
  scale_x_log10(TeX('\\beta'),breaks=c(1e-6,1e-5,0.0001,0.001,0.01,0.1,1,10,100,1000,10000,10000),
                limits=c(1e-6,2.5e-2))+
  scale_y_continuous(
    # Features of the first axis
    #name = TeX('$\\sum_{v\\in H_3} \\delta_v$ (Gy)') ,
    name = TeX('$\\|y\\|_1$ (Gy)') ,
    breaks=c(0,10000,20000,30000,40000,50000),
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~(.*100/Brain_voxels),#name=TeX('$\\frac{\\sum_{v\\in H_3} 1(\\delta_v>0)}{|H_3|}$ (\\%)'))
                        name=TeX('$\\frac{\\|y\\|_0}{|H_2|}$ (\\%)'))
  ) +
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 20),
    axis.title.y = element_text(color = "blue", size=13),
    axis.title.y.right = element_text(color = "red", size=13),
    legend.title = element_blank(),#element_text(size = 13), 
    legend.text  = element_text(size = 13),
    legend.background = element_rect(color="black"),
    legend.position=c(0.83, 0.9),
    legend.margin = margin(-7, 6, 5, 6),
    legend.direction="vertical"
  )+
  geom_hline(yintercept = 1*Brain_voxels/100,linetype="dotted",size=1)+
  annotate("segment", x = 1e-4, xend = df2[df2$method=="Nominal",]$min_beta, y = df2[df2$method=="Nominal",]$l1+20000, yend = 1*Brain_voxels/100,
         size=0.75,linetype="solid",arrow = arrow(type = "closed",length = unit(.2,"cm")))+
  annotate("text", x =1e-4, y = df2[df2$method=="Nominal",]$l1+20000, size=5, hjust=0, label = TeX(paste("\\beta=",sprintf("%2.2e",df2[df2$method=="Nominal",]$min_beta),sep="")))+
  annotate("segment", x = df2[df2$method=="Nominal",]$min_beta , xend = df2[df2$method=="Robust",]$min_beta, y = (df2[df2$method=="Robust",]$l1+30000), yend = 1*Brain_voxels/100,
           size=0.75,linetype="dashed",arrow = arrow(type = "closed",length = unit(.2,"cm")))+
  annotate("text", df2[df2$method=="Nominal",]$min_beta, y = (df2[df2$method=="Robust",]$l1+31000), size=5, hjust=0,vjust=0,color='black',label = TeX(paste("\\beta=",sprintf("%2.2e",df2[df2$method=="Robust",]$min_beta),sep="")))
outputfile=sprintf("./%s/%s",PlotFolder,"Patient4_Visit1_violation_sum_and_num_vs_beta_nominal_robust.png")
ggsave(outputfile,height=4.2, width=6)
dev.off()
#end figure 12


#CVaR Results
ggplot(df1,aes(x=l0/Brain_voxels*100))+
  geom_line(aes(y=g,linetype=method),color='blue',size=1)+
  geom_vline(xintercept = 1,linetype="dotted",size=1)+
  #CVAR results
  geom_point(aes(x=0.004096511*100,y=56.46320733),color='red',shape='+',size=5)+  
  geom_point(aes(x=5612/Brain_voxels*100,y=51.3845907499786),color='red',shape='x',size=3)+
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 20),
    axis.title.y = element_text(size=13),
    axis.title.x = element_text(size=13),
    legend.title = element_blank(),#element_text(size = 13), 
    legend.text  = element_text(size = 13),
    legend.background = element_rect(color="black"),
    legend.position=c(0.12, 0.88),
    legend.margin = margin(-7, 6, 5, 6),
    legend.direction="vertical"
  )+
  scale_y_continuous(name = TeX('\\underline{d} (Gy)'))+
  scale_x_log10(
    name = TeX('$\\frac{\\|y\\|_0}{|H_2|}$ (\\%)') ,
    #breaks=c(0,10000,20000,30000,40000,50000),
  )

  