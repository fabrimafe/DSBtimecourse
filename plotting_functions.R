
library(data.table)
library(tidyverse)
library(ggplot2)
library(forcats)

minimum_exp=9
import_flowdata<-function(myfileflow)
    {
    xdataflow<-fread(myfileflow,header=FALSE)
    names(xdataflow)<-c("rate","flow")
    xdataflow$rate[xdataflow$rate=="k11.k11"]<-"k11"
    xdataflow$model<-mymodel
    xdataflow$target<-mytarget
    xdataAIC<-xdataflow %>% filter(rate=="AIC")
    xdataflow<-xdataflow %>% filter(rate!="AIC")
    return(list(xdataAIC,xdataflow))
    }

parse_rates<-function(xdata.all)
{
xdata.all$rate[xdata.all$rate=="k11"]<<-"Cutting"
xdata.all$rate[xdata.all$rate=="r11"]<<-"Precise Repair (from direct DSB)"
xdata.all$rate[xdata.all$rate=="r12"]<<-"Repair-error (from direct DSB)"
xdata.all$rate[xdata.all$rate=="r21"]<<-"Precise Repair (from processed DSB)"
xdata.all$rate[xdata.all$rate=="r22"]<<-"Repair-error (from processed DSB)"
xdata.all$rate[xdata.all$rate=="K"]<<-"u (uncut fraction)"
xdata.all$rate[xdata.all$rate=="rr12"]<<-"Processing"
xdata.all$rate[xdata.all$rate=="er1"]<<-"Error"
xdata.all$rate[xdata.all$rate=="r0"]<<-"r (induction speed)"
xdata.all$rate[xdata.all$rate=="r2"]<<-"d (induction decay)"
xdata.all$target[xdata.all$target=="PhyB2.2.nov2021"]<<-"PhyB2"
xdata.all$target[xdata.all$target=="CRTISO.nov2021.49and50"]<<-"CRTISO"
}

parse_flows<-function(xdataflow.all)
{
xdataflow.all$rate[xdataflow.all$rate=="k11"]<<-"Cutting"
xdataflow.all$rate[xdataflow.all$rate=="r11"]<<-"Precise Repair (from direct DSB)"
xdataflow.all$rate[xdataflow.all$rate=="r12"]<<-"Repair-error (from direct DSB)"
xdataflow.all$rate[xdataflow.all$rate=="r21"]<<-"Precise Repair (from processed DSB)"
xdataflow.all$rate[xdataflow.all$rate=="r22"]<<-"Repair-error (from processed DSB)"
xdataflow.all$rate[xdataflow.all$rate=="K"]<<-"u (uncut fraction)"
xdataflow.all$rate[xdataflow.all$rate=="rr12"]<<-"Processing"
xdataflow.all$rate[xdataflow.all$rate=="er1"]<<-"Error"
xdataflow.all$rate[xdataflow.all$rate=="r0"]<<-"r (induction speed)"
xdataflow.all$rate[xdataflow.all$rate=="r2"]<<-"d (induction decay)"
xdataflow.all$target[xdataflow.all$target=="PhyB2.2.nov2021"]<<-"PhyB2"
xdataflow.all$target[xdataflow.all$target=="CRTISO.nov2021.49and50"]<<-"CRTISO"
}

parse_flows_all<-function(xdata.allall.table)
{
xdata.allall.table$rate[xdata.allall.table$rate=="k11"]<<-"Cutting"
xdata.allall.table$rate[xdata.allall.table$rate=="r11"]<<-"Precise Repair (from direct DSB)"
xdata.allall.table$rate[xdata.allall.table$rate=="r12"]<<-"Repair-error (from direct DSB)"
xdata.allall.table$rate[xdata.allall.table$rate=="r21"]<<-"Precise Repair (from processed DSB)"
xdata.allall.table$rate[xdata.allall.table$rate=="r22"]<<-"Repair-error (from processed DSB)"
xdata.allall.table$rate[xdata.allall.table$rate=="K"]<<-"u (uncut fraction)"
xdata.allall.table$rate[xdata.allall.table$rate=="rr12"]<<-"Processing"
xdata.allall.table$rate[xdata.allall.table$rate=="er1"]<<-"Error"
xdata.allall.table$rate[xdata.allall.table$rate=="r0"]<<-"r (induction speed)"
xdata.allall.table$rate[xdata.allall.table$rate=="r2"]<<-"d (induction decay)"
xdata.allall.table$target[xdata.allall.table$target=="PhyB2.2.nov2021"]<<-"PhyB2"
xdata.allall.table$target[xdata.allall.table$target=="CRTISO.nov2021.49and50"]<<-"CRTISO"
}

parse_rates_mean<-function(xdatamean.all,orderforplotting=orderforplotting,minimum_exp=9)
{
xdatamean.all$rates<-xdatamean.all$max
xdatamean.all$rate<-factor(xdatamean.all$rate,levels=orderforplotting)
xdatamean.all$rates_trans<-log(xdatamean.all$rates+10^(-minimum_exp),10)
xdatamean.all<-xdatamean.all %>% filter(rate!="r0",rate!="r2")
xdatamean.all$target<-fct_inorder(as.factor(as.character(xdatamean.all$target)))
return(xdatamean.all)
}

parse_rates_all<-function(xdata.all,orderforplotting=orderforplotting,mypalette=mypalette,induction_params=induction_params)
{
xdata.all$rate<-factor(xdata.all$rate,levels=orderforplotting)
xdata.all<-xdata.all %>% arrange(rate)
xdata.all$rate<-fct_inorder(as.factor(as.character(xdata.all$rate)))
temp<-xdata.all$rate %>% unique
#mypalette_t<-mypalette[orderforplotting %in% temp]
xdata.all$target<-fct_inorder(as.factor(as.character(xdata.all$target)))
#mypalette_t<-mypalette[orderforplotting %in% (xdata.all$rate %>% unique) & !orderforplotting %in% induction_params]
return(xdata.all)
#return(list(xdata.all,mypalette_t))
}
#parse_rates_all(xdata.all)

set_breaks<-function(minimum_exp=9)
{
xbreaks<-sapply(minimum_exp:(-1),function(x) 10^(-x))
xlabels<-xbreaks
xbreaks<-log(xbreaks+10^(-minimum_exp),10)
xbreaks[1]<-log(10^(-minimum_exp),10)
xlabels[1]<-0 #
return(list(xbreaks,xlabels))
}

plot_violin_params<-function(xdata.all.temp,xdatamean.all.temp,induction_params=induction_params,outputfile="plot1.pdf")
{
p<-ggplot(xdata.all.temp, aes(x=rate, y=rates_trans, fill=rate)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=mypalette_t) +
theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 45, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),legend.position = "bottom",aspect.ratio=1) + ylab("rate ( event/hour )") +
facet_wrap(~target) + scale_y_continuous(breaks=xbreaks,labels=xlabels,limits=c(min(xbreaks)-0.01,log(20,10)))+
geom_point(data=xdatamean.all.temp,aes(rate, rates_trans),col = "black",size=0.8)+geom_jitter(alpha = 0.11, width = 0.15, size=0.2) #, size=1
ggsave(p,file=outputfile)
}

plot_violin_induction<-function(xdata.all.temp,xdatamean.all.temp,induction_params=induction_params,outputfile="plot2.pdf")
{
p<-ggplot(xdata.all.temp, aes(x=rate, y=rates_trans, fill=rate)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=mypalette_t) +
theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 45, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),legend.position = "bottom",aspect.ratio=1) + ylab("rate ( event/hour )") +
facet_wrap(~target) + scale_y_continuous(breaks=xbreaks_ind,labels=xlabels_ind,limits=c(min(xbreaks_ind)-0.01,log(10000,10)))+
geom_point(data=xdatamean.all.temp,aes(rate, rates_trans),col = "black",size=0.8)+geom_jitter(alpha = 0.11, width = 0.15, size=0.2) #, size=1
ggsave(p,file=outputfile)
}

plot_violin_induction_wrapper<-function(){
xdata.all.temp<-xdata.all %>% filter(rate %in% induction_params)
xdatamean.all.temp<-xdatamean.all %>% filter(rate %in% induction_params)
xdata.all.temp$rates_trans[xdata.all.temp$rates_trans<(-6.5) & xdata.all.temp$rate=="d (induction decay)"]<-(-9)
xdatamean.all.temp$rates_trans[xdatamean.all.temp$rates_trans<(-6.5) & xdatamean.all.temp$rate=="d (induction decay)"]<-(-9)
mypalette_t<-mypalette[orderforplotting %in% (xdata.all$rate %>% unique) & orderforplotting %in% induction_params]
plot_violin_induction(xdata.all.temp,xdatamean.all.temp,outputfile=paste0("violininduction_",mymodel,"_indc",ninduction,".pdf"))
}

plot_flow<-function()
{
pflow<-ggplot(xdataflow.all.temp, aes(x=rate, y=flow, fill=rate)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=mypalette_t)+
theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 45, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),legend.position = "bottom",aspect.ratio=1) + ylab("flow ( prop. of initial molecules )") +
scale_y_continuous(trans = "sqrt",breaks=c(0,0.1,0.5,1,5,10,15,20),limits=c(0,22))+facet_wrap("target")+
geom_point(data=xdata.allall.table.temp,aes(rate, flow),col = "black",size=0.8)+
geom_jitter(alpha = 0.11, width = 0.15, size=0.2)
ggsave(pflow,file=paste0("violinflow_",mymodel,"_indc",ninduction,".pdf"))
}

plots_AIC<-function()
{
AICnok12<-xdataAIC.all %>% filter(model=="modelDSBs1i1_nok12") %>% select(flow)
AIC3x4<-xdataAIC.all %>% filter(model=="modelDSBs1i1_3x4") %>% select(flow)
AICtarget<-xdataAIC.all %>% filter(model=="modelDSBs1i1_nok12") %>% select(target)
AIC_plot<-data.frame(cbind(AICtarget,deltaAIC=AICnok12-AIC3x4))
names(AIC_plot)<-c("target","deltaAIC")
AIC_plot<-AIC_plot %>% mutate(relL=exp(0.5*deltaAIC))
AIC_plot$target[AIC_plot$target=="PhyB2.2.nov2021"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="CRTISO.nov2021.49and50"]<-"CRTISO"

pAIC2<-ggplot(AIC_plot, aes(x=relL)) + geom_histogram(color="darkblue", fill="lightblue")+facet_wrap("target")+
theme(strip.background = element_rect(fill = "white"), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),
axis.text.x = element_text(angle = 90),aspect.ratio=1) +
ylab("density") + xlab("relative likelihood in support of 4state vs 3state model")+
scale_x_continuous(breaks=c(0.05,1,2,3))
ggsave(pAIC2,file=paste0("histrelL_indc",ninduction,".pdf"))

pAIC<-ggplot(AIC_plot, aes(x=deltaAIC)) + geom_histogram(color="darkblue", fill="lightblue",breaks=seq(-500,100,25))+facet_wrap("target")+
theme(strip.background = element_rect(fill = "white"), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),aspect.ratio=1) + ylab("density") + xlab("delta AIC (4state - 3state model)")
ggsave(pAIC,file=paste0("histdeltaAIC_indc",ninduction,".pdf"))

res1<-AIC_plot %>% filter(target=="PhyB2") %>% filter(relL>0.00001) #4
res2<-AIC_plot %>% filter(target=="Psy1") %>% filter(relL>0.00001) #0
res3<-AIC_plot %>% filter(target=="CRTISO") %>% filter(relL>0.00001) #19
res4<-AIC_plot %>% filter(target=="PhyB2") %>% filter(deltaAIC>0) #0
res5<-AIC_plot %>% filter(target=="Psy1") %>% filter(deltaAIC>0) #0
res6<-AIC_plot %>% filter(target=="CRTISO") %>% filter(deltaAIC>0) #0
print(c(res1,res2,res3,res4,res5,res6))
}


xbreaks<-set_breaks()[[1]]
xlabels<-set_breaks()[[2]]
xbreaks_ind<-set_breaks()[[1]]
xlabels_ind<-set_breaks()[[2]]

mypalette=c("firebrick","darkorange","dodgerblue2","dodgerblue4","green4","darkgreen","burlywood4","darksalmon","deeppink3","deeppink4")
mypalette=c("Cutting"="firebrick","Processing"="darkorange1","Repair-error (from direct DSB)"="deepskyblue2","Repair-error (from processed DSB)"="dodgerblue4",
"Precise Repair (from direct DSB)"="chartreuse3","Precise Repair (from processed DSB)"="darkgreen","Error"="gray65","u (uncut fraction)"="burlywood4",
"r (induction speed)"="orchid3","d (induction decay)"="purple3")
orderforplotting=c("Cutting","Processing","Repair-error (from direct DSB)","Repair-error (from processed DSB)","Precise Repair (from direct DSB)","Precise Repair (from processed DSB)","Error","u (uncut fraction)","r (induction speed)","d (induction decay)")
induction_params<-c("Error","u (uncut fraction)","r (induction speed)","d (induction decay)")



