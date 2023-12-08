
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

parse_rate_tables<-function(x)
{
y<-x %>% 
mutate( model=case_when(
  model == "modelDSBs1i1_3x4" ~ "3 states",
  model == "modelDSBs1i1_nok12" ~ "4 states",
  .default = as.character(model)  
)) %>%
mutate( target=case_when(
  target == "PhyB2.2.nov2021" ~ "PhyB2",
  target == "CRTISO.nov2021.49and50" ~ "CRTISO",  
  target == "Psy1.20231005.I72h" ~ "Psy1",
  target == "PhyB2.20231005.I72h" ~ "PhyB2",
  target == "CRTISO.20231005.I72h" ~ "CRTISO(-4bp Processed DSB)",
  target == "CRTISO.49and50.20231005.I72h" ~ "CRTISO",
  target == "Psy1.20231005" ~ "Psy1",
  target == "CRTISO.20230903" ~ "CRTISO(-4bp Processed DSB)",
  target == "CRTISO.49and50.20230903" ~ "CRTISO",
  target == "PhyB2.20230903" ~ "PhyB2",
  target == "Psy1.20231026" ~ "Psy1",
  target == "CRTISO.20231026" ~ "CRTISO(-4bp Processed DSB)",
  target == "CRTISO.49and50.20231026" ~ "CRTISO",
  target == "PhyB2.20231026" ~ "PhyB2",
  target == "Psy1.20230910" ~ "Psy1",
  target == "CRTISO.20230910" ~ "CRTISO(-4bp Processed DSB)",
  target == "CRTISO.49and50.20230910" ~ "CRTISO",
  target == "PhyB2.20230910" ~ "PhyB2",
  target == "Psy1.20231101" ~ "Psy1",
  target == "Psy1.20231102" ~ "Psy1",
  target == "Psy1.20230909" ~ "Psy1",
  target == "CRTISO.20231102" ~ "CRTISO(-4bp Processed DSB)",
  target == "CRTISO.49and50.20231102" ~ "CRTISO",
  target == "CRTISO.49and50.20230909" ~ "CRTISO",
  target == "PhyB2.20231102" ~ "PhyB2",
  target == "PhyB2.20231115" ~ "PhyB2",
  target == "PhyB2.20231116" ~ "PhyB2",
  target == "PhyB2.20231117" ~ "PhyB2",
  target == "PhyB2.cas12" ~ "PhyB2",
  target == "PhyB2.20231130" ~ "PhyB2",
  target == "Psy1.20231130" ~ "Psy1",
  target == "CRTISO.49and50.20231130" ~ "CRTISO.49and50",
  target == "PhyB2.20231130.I72h" ~ "PhyB2",
  target == "Psy1.20231130.I72h" ~ "Psy1",
  target == "CRTISO.49and50.20231130.I72h" ~ "CRTISO.49and50",
  .default = as.character(target)
)) %>% 
mutate( rate=case_when(
  rate == "k11" ~ "Cutting",
  rate == "r11" & model == "3 states" ~ "Precise Repair",
  rate == "r12" & model == "3 states" ~ "Repair-error",
  rate == "r11" & model == "4 states" ~ "Precise Repair (from direct DSB)",
  rate == "r12" & model == "4 states" ~ "Repair-error (from direct DSB)",
  rate == "r21" ~ "Precise Repair (from processed DSB)",
  rate == "r22" ~ "Repair-error (from processed DSB)",
  rate == "K" ~ "u (uncut fraction)",
  rate == "rr12" ~ "Processing",
  rate == "r0" ~ "r (induction speed)",
  rate == "r2" ~ "d (induction decay)",
  rate == "er1" ~ "Error",
  .default = as.character(rate)
))
return(y)
}



parse_rates_mean<-function(xdatamean.all,orderforplotting=orderforplotting,minimum_exp=9)
{
xdatamean.all$rates<<-xdatamean.all$max
xdatamean.all$rates_trans<<-log(xdatamean.all$rates+10^(-minimum_exp),10)
xdatamean.all<<-xdatamean.all %>% filter(rate!="r0",rate!="r2")
xdatamean.all$rate<<-factor(xdatamean.all$rate,levels=orderforplotting)
xdatamean.all$target<<-fct_inorder(as.factor(as.character(xdatamean.all$target)))
#return(xdatamean.all)
}

parse_rates_all<-function(x,orderforplotting=orderforplotting,mypalette=mypalette,induction_params=induction_params)
{
y<-x
y$rates[y$rates<0]<-0
y$rates_trans<-log(y$rates+10^(-minimum_exp),10)
y$rate<-factor(y$rate,levels=orderforplotting)
y<-y %>% arrange(rate)
y$rate<-fct_inorder(as.factor(as.character(y$rate)))
temp<-y$rate %>% unique
#mypalette_t<-mypalette[orderforplotting %in% temp]
y$target<-fct_inorder(as.factor(as.character(y$target)))
#mypalette_t<-mypalette[orderforplotting %in% (xdata.all$rate %>% unique) & !orderforplotting %in% induction_params]
return(y)
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

plot_violin_induction<-function(xdata.all,xdatamean.all,induction_params=induction_params,outputfile="plot2.pdf")
{
#xdata.all.temp<-xdata.all %>% filter(rate %in% induction_params)
#xdatamean.all.temp<-xdatamean.all %>% filter(rate %in% induction_params)
#xdata.all.temp$rates_trans[xdata.all.temp$rates_trans<(-6.5) & xdata.all.temp$rate=="d (induction decay)"]<-(-9)
#xdatamean.all.temp$rates_trans[xdatamean.all.temp$rates_trans<(-6.5) & xdatamean.all.temp$rate=="d (induction decay)"]<-(-9)
#mypalette_t<-mypalette[orderforplotting %in% (xdata.all$rate %>% unique) & orderforplotting %in% induction_params]
#mypalette_t<-mypalette[orderforplotting %in% induction_params]
#plot_violin_induction(xdata.all.temp,xdatamean.all.temp,outputfile=paste0("violininduction_",mymodel,"_indc",ninduction,".pdf"))
p<-ggplot(xdata.all.temp, aes(x=rate, y=rates_trans, fill=rate)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=mypalette_t) +
theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 45, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),legend.position = "bottom",aspect.ratio=1) + ylab("rate ( event/hour )") +
facet_wrap(~target) + scale_y_continuous(breaks=xbreaks_ind,labels=xlabels_ind,limits=c(min(xbreaks_ind)-0.01,log(10000,10)))+
geom_point(data=xdatamean.all.temp,aes(rate, rates_trans),col = "black",size=0.8)+geom_jitter(alpha = 0.11, width = 0.15, size=0.2) #, size=1
ggsave(p,file=outputfile)
}

plot_flow<-function(outputfile=paste0("violinflow_",summarytag,"_",mymodel,"_indc",ninduction,".pdf"))
{
pflow<-ggplot(xdataflow.all.temp, aes(x=rate, y=flow, fill=rate)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=mypalette_t)+
theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 45, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),legend.position = "bottom",aspect.ratio=1) + ylab("flow ( prop. of initial molecules )") +
scale_y_continuous(trans = "sqrt",breaks=c(0,0.1,0.5,1,5,10),limits=c(0,11))+facet_wrap("target")+
geom_point(data=xdata.allall.table.temp,aes(rate, flow),col = "black",size=0.8)+
geom_jitter(alpha = 0.11, width = 0.15, size=0.2)
ggsave(pflow,file=outputfile)
}

plots_AIC<-function(outputfile=paste0("histdeltaAIC_",summarytag,"_indc",ninduction,".pdf"))
{
AICnok12<-xdataAIC.all %>% filter(model=="modelDSBs1i1_nok12") %>% select(flow)
AIC3x4<-xdataAIC.all %>% filter(model=="modelDSBs1i1_3x4") %>% select(flow)
AICtarget<-xdataAIC.all %>% filter(model=="modelDSBs1i1_nok12") %>% select(target)
AIC_plot<-data.frame(cbind(AICtarget,deltaAIC=AICnok12-AIC3x4))
names(AIC_plot)<-c("target","deltaAIC")
AIC_plot<-AIC_plot %>% mutate(relL=exp(0.5*deltaAIC))
AIC_plot$target[AIC_plot$target=="PhyB2.2.nov2021"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="CRTISO.nov2021.49and50"]<-"CRTISO"
AIC_plot$target[AIC_plot$target=="PhyB2.2.nov2021"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="CRTISO.20231005.I72h"]<-"CRTISO(-4bp Processed DSB)"
AIC_plot$target[AIC_plot$target=="CRTISO.49and50.20231005.I72h"]<-"CRTISO"
AIC_plot$target[AIC_plot$target=="Psy1.20231005.I72h"]<-"Psy1"
AIC_plot$target[AIC_plot$target=="PhyB2.20231005.I72h"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="CRTISO.20230903"]<-"CRTISO(-4bp Processed DSB)"
AIC_plot$target[AIC_plot$target=="CRTISO.49and50.20230903"]<-"CRTISO"
AIC_plot$target[AIC_plot$target=="Psy1.20231005"]<-"Psy1"
AIC_plot$target[AIC_plot$target=="PhyB2.20230903"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="CRTISO.20231026"]<-"CRTISO(-4bp Processed DSB)"
AIC_plot$target[AIC_plot$target=="CRTISO.49and50.20231026"]<-"CRTISO"
AIC_plot$target[AIC_plot$target=="Psy1.20231101"]<-"Psy1"
AIC_plot$target[AIC_plot$target=="Psy1.20231102"]<-"Psy1"
AIC_plot$target[AIC_plot$target=="PhyB2.20231102"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="CRTISO.20231102"]<-"CRTISO(-4bp Processed DSB)"
AIC_plot$target[AIC_plot$target=="CRTISO.49and50.20231102"]<-"CRTISO"
AIC_plot$target[AIC_plot$target=="Psy1.20231026"]<-"Psy1"
AIC_plot$target[AIC_plot$target=="PhyB2.20231026"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="CRTISO.20230910"]<-"CRTISO(-4bp Processed DSB)"
AIC_plot$target[AIC_plot$target=="CRTISO.49and50.20230910"]<-"CRTISO"
AIC_plot$target[AIC_plot$target=="CRTISO.49and50.20230909"]<-"CRTISO"
AIC_plot$target[AIC_plot$target=="Psy1.20230909"]<-"Psy1"
AIC_plot$target[AIC_plot$target=="Psy1.20230910"]<-"Psy1"
AIC_plot$target[AIC_plot$target=="PhyB2.20230910"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="PhyB2.20231115"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="PhyB2.20231116"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="PhyB2.20231117"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="PhyB2.cas12"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="PhyB2.20231130"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="Psy1.20231130"]<-"Psy1"
AIC_plot$target[AIC_plot$target=="CRTISO.49and50.20231130"]<-"CRTISO.49and50"
AIC_plot$target[AIC_plot$target=="PhyB2.20231130.I72h"]<-"PhyB2"
AIC_plot$target[AIC_plot$target=="Psy1.20231130.I72h"]<-"Psy1"
AIC_plot$target[AIC_plot$target=="CRTISO.49and50.20231130.I72h"]<-"CRTISO.49and50"

#plot relative likelihood (not used)
pAIC2<-ggplot(AIC_plot, aes(x=relL)) + geom_histogram(color="darkblue", fill="lightblue")+facet_wrap("target")+
theme(strip.background = element_rect(fill = "white"), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),
axis.text.x = element_text(angle = 90),aspect.ratio=1) +
ylab("density") + xlab("relative likelihood in support of 4state vs 3state model")+
scale_x_continuous(breaks=c(0.05,1,2,3))
#ggsave(pAIC2,file=outputfile)
#plot delta AIC
pAIC<-ggplot(AIC_plot, aes(x=deltaAIC)) + geom_histogram(color="darkblue", fill="lightblue",breaks=seq(-500,100,25))+facet_wrap("target")+
theme(strip.background = element_rect(fill = "white"), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),aspect.ratio=1) + ylab("density") + xlab("delta AIC (4state - 3state model)")
ggsave(pAIC,file=outputfile)

res1<-AIC_plot %>% filter(target=="PhyB2") %>% filter(relL>0.00001) #4
res2<-AIC_plot %>% filter(target=="Psy1") %>% filter(relL>0.00001) #0
res3<-AIC_plot %>% filter(target=="CRTISO") %>% filter(relL>0.00001) #19
res4<-AIC_plot %>% filter(target=="PhyB2") %>% filter(deltaAIC>0) #0
res5<-AIC_plot %>% filter(target=="Psy1") %>% filter(deltaAIC>0) #0
res6<-AIC_plot %>% filter(target=="CRTISO") %>% filter(deltaAIC>0) #0
print(c(res1,res2,res3,res4,res5,res6))
}


set_filenames<-function(beforejune2023=FALSE)
{
    myfile<<-paste0("/home/labs/alevy/fabrizio/workspace/daniela/resultsv5/",bootstrap_type,"/timecourse_RNP_",mytarget,"/",mymodel,"/results_",mymodel,"_RNP_ind.c",ninduction,"_",mytarget,"/results.tot")
    myfilemean<<-paste0("/home/labs/alevy/fabrizio/workspace/daniela/resultsv5/results_",mymodel,"_RNP_ind.c",ninduction,"_",mytarget,".all_n50000.CI")
    if (beforejune2023)
      	{
      	myfileflow<<-paste0("/home/labs/alevy/fabrizio/workspace/daniela/resultsv5/",bootstrap_type,"/timecourse_RNP_",mytarget,"/",mymodel,"/results_",mymodel,"_RNP_ind.c",ninduction,"_",mytarget,"/results.tot_plot.flow.tab")
      myfileflowmean<<-paste0("/home/labs/alevy/fabrizio/workspace/daniela/resultsv5/plots/plot_",mymodel,"_RNP_ind.c",ninduction,"_",mytarget,"_plot.flow.tab")
      	} else
      	{
      	myfileflow<<-paste0("/home/labs/alevy/fabrizio/workspace/daniela/resultsv5/",bootstrap_type,"/timecourse_RNP_",mytarget,"/",mymodel,"/results_",mymodel,"_RNP_ind.c",ninduction,"_",mytarget,"/results.totall_",timeflow,"h.plot.flow.tab")
      	myfileflowmean<<-paste0("/home/labs/alevy/fabrizio/workspace/daniela/resultsv5/plots/plot_",mymodel,"_RNP_ind.c",ninduction,"_",mytarget,"_",timeflow,"_plot.flow.tab")
      	}
}

write_table_flow<-function(outputfile="table1.tab",orderforplotting=orderforplotting,typetable=1)
{
print(orderforplotting)
xdata.allall<-parse_rate_tables(xdata.allrates.all) %>% left_join(parse_rate_tables(xdata.allflows.all),by=c("model","target","rate","bootstrap"))
#rate == "r11" & model == "4 states" ~ "Precise Repair (from direct DSB)",
#  rate == "r12" & model == "4 states" ~ "Repair-error (from direct DSB)",
#  rate == "r21" ~ "Precise Repair (from processed DSB)",
#  rate == "r22" ~ "Repair-error (from processed DSB)",
if (xdata.allall$model[1]=="4 states")
	{
	adata<-xdata.allall %>% group_by(model,target,rate) %>% filter(rate=="Precise Repair (from direct DSB)") 
	bdata<-xdata.allall %>% group_by(model,target,rate) %>% filter(rate=="Precise Repair (from processed DSB)") 
	adata[,6:ncol(adata)]<-adata[,6:ncol(adata)]+bdata[,6:ncol(bdata)]
	adata$rate<-"Precise Repair (total)"
	xdata.allall<-rbind(xdata.allall,adata)
	adata<-xdata.allall %>% group_by(model,target,rate) %>% filter(rate=="Repair-error (from direct DSB)") 
	bdata<-xdata.allall %>% group_by(model,target,rate) %>% filter(rate=="Repair-error (from processed DSB)") 
	adata[,6:ncol(adata)]<-adata[,6:ncol(adata)]+bdata[,6:ncol(bdata)]
	adata$rate<-"Repair-error (total)"
	xdata.allall<-rbind(xdata.allall,adata)
	adata<-xdata.allall %>% group_by(model,target,rate) %>% filter(rate=="Precise Repair (total)") 
	bdata<-xdata.allall %>% group_by(model,target,rate) %>% filter(rate=="Repair-error (total)") 
	adata[,6:ncol(adata)]<-adata[,6:ncol(adata)]/(adata[,6:ncol(adata)]+bdata[,6:ncol(bdata)])
	adata$rate<-"Repair accuracy (total)"
	xdata.allall<-rbind(xdata.allall,adata)
	orderforplotting<-c(orderforplotting,"Precise Repair (total)","Repair-error (total)","Repair accuracy (total)")
	} else if (xdata.allall$model[1]=="3 states")
	{
	adata<-xdata.allall %>% group_by(model,target,rate) %>% filter(rate=="Precise Repair") 
	bdata<-xdata.allall %>% group_by(model,target,rate) %>% filter(rate=="Repair-error") 
	adata[,6:ncol(adata)]<-adata[,6:ncol(adata)]/(adata[,6:ncol(adata)]+bdata[,6:ncol(bdata)])
	adata$rate<-"Repair accuracy (total)"
	xdata.allall<-rbind(xdata.allall,adata)
	orderforplotting<-c(orderforplotting,"Repair accuracy (total)")
	
	}
xdata.allall.table<-xdata.allall %>% group_by(model,target,rate) %>% summarise(maxx=mean(max),rate01=quantile(rates,0.01),rate05=quantile(rates,0.05),rate5=quantile(rates,0.5),rate95=quantile(rates,0.95),rate99=quantile(rates,0.99), flow_mean=mean(flow_mean),flow01=quantile(flow,0.01,na.rm=T),flow05=quantile(flow,0.05,na.rm=T),flow5=quantile(flow,0.5,na.rm=T),flow95=quantile(flow,0.95,na.rm=T),flow99=quantile(flow,0.99,na.rm=T), pvalue=sum(rates<=0)/sum(rates>=0)) %>% ungroup
xdata.allall.table$rate<-factor(xdata.allall.table$rate,levels=orderforplotting)
xdata.allall.table<-xdata.allall.table %>% #group_by(model) %>% 
arrange(factor(target,levels=c("Psy1","CRTISO","CRTISO(-4bp Processed DSB)","PhyB2")),
#factor(rate,levels=c(orderforplotting)),
rate) # %>% #,
#factor(model,levels=c("3 states","4 states")),
#.by_group = TRUE) %>% #ungroup %>% 
if (typetable==1)
	{
	print("format of table is 1: only 1-99% CI displayed as CIlow-CI-high") #as for main text
	xdata.allall.table <- xdata.allall.table %>% mutate(across(where(is.numeric), ~ ifelse(.<0,0,.)))
	xdata.allall.table <- xdata.allall.table %>% mutate(across(where(is.numeric), ~ round(., digits = 4))) 
	xdata.allall.table<-xdata.allall.table %>% select(model,target,rate,maxx,rate01,rate99,flow_mean,flow01,flow99,pvalue)
	xdata.allall.table$rate01=paste(xdata.allall.table$rate01,xdata.allall.table$rate99,sep="-")
	xdata.allall.table$rate99=NULL
	xdata.allall.table$flow01=paste(xdata.allall.table$flow01,xdata.allall.table$flow99,sep="-")
	xdata.allall.table$flow99=NULL
	} else 
if (typetable==2)
	{
	print("format of table is 2: 1,5,95,99% CI; if 4 states also plot overall repair")
	xdata.allall.table<-xdata.allall.table %>% select(model,target,rate,maxx,rate01,rate05,rate95,rate99,flow_mean,flow01,flow05,flow95,flow99,pvalue)
	xdata.allall.table <- xdata.allall.table %>% mutate(across(where(is.numeric), ~ ifelse(.<0,0,.)))
	#xdata.allall.table[xdata.allall.table<0]<-0
	xdata.allall.table <- xdata.allall.table %>% mutate(across(where(is.numeric), ~ round(., digits = 4)))
	}
write.csv(xdata.allall.table,file=outputfile,row.names=F,quote=FALSE) #sep='\t'
}



xbreaks<-set_breaks()[[1]]
xlabels<-set_breaks()[[2]]
xbreaks_ind<-set_breaks()[[1]]
xlabels_ind<-set_breaks()[[2]]

mypalette=c("firebrick","darkorange","dodgerblue2","dodgerblue4","green4","darkgreen","burlywood4","darksalmon","deeppink3","deeppink4")

set_model<-function(mymodel)
{
if (mymodel=="modelDSBs1i1_3x4") 
    {
    mypalette<<-c("Cutting"="firebrick","Repair-error"="dodgerblue2","Precise Repair"="chartreuse3","Error"="gray65","u (uncut fraction)"="burlywood4","r (induction speed)"="orchid3","d (induction decay)"="purple3")
    orderforplotting<<-c("Cutting","Repair-error","Precise Repair","Error","u (uncut fraction)","r (induction speed)","d (induction decay)")
    induction_params<<-c("Error","u (uncut fraction)","r (induction speed)","d (induction decay)")
    } else
if (mymodel=="modelDSBs1i1_nok12")
    {
    mypalette<<-c("Cutting"="firebrick","Processing"="darkorange1","Repair-error (from direct DSB)"="deepskyblue2","Repair-error (from processed DSB)"="dodgerblue4",
"Precise Repair (from direct DSB)"="chartreuse3","Precise Repair (from processed DSB)"="darkgreen","u (uncut fraction)"="burlywood4",
"r (induction speed)"="orchid3","d (induction decay)"="purple3")
    orderforplotting<<-c("Cutting","Processing","Repair-error (from direct DSB)","Repair-error (from processed DSB)","Precise Repair (from direct DSB)","Precise Repair (from processed DSB)","u (uncut fraction)","r (induction speed)","d (induction decay)")
    induction_params<<-c("u (uncut fraction)","r (induction speed)","d (induction decay)")
    }
}


