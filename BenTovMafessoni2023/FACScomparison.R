#CODE TO GENERATE THE COMPARISON BETWEEN ESTIMATES FOLLOWING INDUCTION CURVES BASED ON FACS AND ESTIMATES FROM THE MODEL

source("plotting_functions.R")
library(tidyverse)
counter<-1
lexperiment<-c(rep("I.72h",3),rep("I.72.5h",3),rep("II.24h",3),rep("II.24.5h",3)) #Loop over the experiments
l_target=c("Psy1","CRTISO.nov2021.49and50","PhyB2.2.nov2021","Psy1.20231005.I72h","CRTISO.49and50.20231005.I72h","PhyB2.20231005.I72h","Psy1.20230910","CRTISO.49and50.20230910","PhyB2.20231115","Psy1.20230909","CRTISO.49and50.20230909","PhyB2.20231117") #Loop over the targets
for (fixedinduction in c(5)){ 
for (xmodel in c("modelDSBs1i1_3x4","modelDSBs1i1_nok12")){ #Loop over the models
for (itarget in 1:length(l_target))
{
###IMPORT ESTIMATES FOR FIXED PARAMS
xtarget<-l_target[itarget]
OUTPUTFOLDER=paste0("~/workspace/daniela/FACS/",fixedinduction,"/results_",xmodel,"_RNP_ind.c3_",xtarget)
OUTPUT=paste0(OUTPUTFOLDER,"/timecourse_rates.all_n50000.CI")
COESTIMATERUN<-paste0("~/workspace/daniela/resultsv5/results_",xmodel,"_RNP_ind.c3_",xtarget,".all_n50000.CI")
xdata_t<-read.table(OUTPUT,header=T,sep='\t')
xdata_t$bootstrap<-0
for ( iboot in 1:100){
OUTPUT=paste0(OUTPUTFOLDER,"/timecourse_rates_bs",iboot,".txt.parsed.tab_n50000.CI")
xdataboot_t<-read.table(OUTPUT,header=T,sep='\t')
xdataboot_t$bootstrap<-iboot
if (iboot==1) {xdataboot<-xdataboot_t} else {xdataboot<-rbind(xdataboot,xdataboot_t)}
}
xdataboot_fi<-xdataboot
xdata_t<-rbind(xdata_t,xdataboot)
xdata_t$model<-xmodel
xdata_t$fixedinduction<-fixedinduction
xdata_t$target<-xtarget
xdata_t$experiment<-lexperiment[itarget]
###IMPORT ESTIMATES FOR COESTIMATED PARAMS
xdata_tco<-read.table(COESTIMATERUN,header=T,sep='\t')
xdata_tco$bootstrap<-0
for ( iboot in 1:100){
OUTPUT=paste0("~/workspace/daniela/resultsv5/stratifiedbootstraps/timecourse_RNP_",xtarget,"/",xmodel,"/results_",xmodel,"_RNP_ind.c3_",xtarget,"/results_",iboot,".txt.parsed.tab_n50000.CI.t.CI") #Loading parameter estimates. Change path if running on other data.
xdataboot_t<-read.table(OUTPUT,header=T,sep='\t')
xdataboot_t$bootstrap<-iboot
if (iboot==1) {xdataboot<-xdataboot_t} else {xdataboot<-rbind(xdataboot,xdataboot_t)}
}
CIs<-xdataboot %>% group_by(rate) %>% summarise(CI5_bootstrap=quantile(max,0.05),CI95_bootstrap=quantile(max,0.95))
xdata_tco<-rbind(xdata_tco,xdataboot)
xdata_t<-cbind(xdata_t,xdata_tco)
names(xdata_t)<-c("estimate_fi","CIlow_fi","CIhigh_fi","rate_fi", "bootstrap_fi","model","fixedinduction","target","experiment","estimate","CIlow","CIhigh","rate","bootstrap")
xdataboot_temp<-xdataboot_fi %>% group_by(rate) %>% summarise(CI5=quantile(max,0.05),CI95=quantile(max,0.95))
xdata_t<-left_join(xdata_t,xdataboot_temp)
if (counter==1) {xdata<-xdata_t} else {xdata<-rbind(xdata,xdata_t)}
counter<-counter+1
}}}
xdata_table<-xdata %>% parse_rate_tables %>% select(experiment,model,target,rate,estimate_fi,CI5,CI95,estimate)
#transform estimates to log scale
minimum_exp<-9
xbreaks<-sapply(minimum_exp:(-1),function(x) 10^(-x))
xlabels<-xbreaks
xbreaks<-log(xbreaks+10^(-minimum_exp),10)
xbreaks[1]<-log(10^(-minimum_exp),10)
xlabels[1]<-0 
xbreaks_ind<-sapply(minimum_exp:(-4),function(x) 10^(-x))
xlabels_ind<-xbreaks_ind
xbreaks_ind<-log(xbreaks_ind+10^(-minimum_exp),10)
xbreaks_ind[1]<-log(10^(-minimum_exp),10)
xlabels_ind[1]<-0 

#Define different palettes and layout for 3 and 4 state models
xdata_table$experiment[xdata_table$experiment=="I.72h"]<-"72h"
xdata_table$experiment[xdata_table$experiment=="I.72.5h"]<-"72h + transformation time"
xdata_table$experiment[xdata_table$experiment=="II.24h"]<-"24h"
xdata_table$experiment[xdata_table$experiment=="II.24.5h"]<-"24h + transformation time"
for (xmodel in c("3 states","4 states")){
for (plotalsoinduction in c(TRUE,FALSE)){
for (myexp in c("I","II")){
if (myexp=="I"){
xdata_table_3states<-xdata_table %>% filter(model==xmodel) %>% filter(experiment=="72h + transformation time" | experiment=="72h")
} else if (myexp=="II"){
xdata_table_3states<-xdata_table %>% filter(model==xmodel) %>% filter(experiment=="24h + transformation time" | experiment=="24h")
}
if (xmodel=="3 states") {
mypalette=c("Cutting"="firebrick","Processing"="darkorange1","Repair-error"="deepskyblue2",
"Precise Repair"="chartreuse3","Error"="gray65","u (uncut fraction)"="burlywood4",
"r (induction speed)"="orchid3","d (induction decay)"="purple3")
orderforplotting=c("Cutting","Processing","Repair-error","Precise Repair","Error","u (uncut fraction)","r (induction speed)","d (induction decay)")
} else if (xmodel=="4 states") 
{
mypalette<-c("Cutting"="firebrick","Processing"="darkorange1","Repair-error (from direct DSB)"="deepskyblue2","Repair-error (from processed DSB)"="dodgerblue4",
"Precise Repair (from direct DSB)"="chartreuse3","Precise Repair (from processed DSB)"="darkgreen","u (uncut fraction)"="burlywood4",
"r (induction speed)"="orchid3","d (induction decay)"="purple3")
orderforplotting<-c("Cutting","Processing","Repair-error (from direct DSB)","Repair-error (from processed DSB)","Precise Repair (from direct DSB)","Precise Repair (from processed DSB)","u (uncut fraction)","r (induction speed)","d (induction decay)")
}
xdata_table_3states_plot<-xdata_table_3states %>% pivot_longer(names_to="induction",values_to="rate_trans",cols=c("estimate","estimate_fi"))
xdata_table_3states_plot$rate_trans[xdata_table_3states_plot$rate_trans<0]<-0
xdata_table_3states_plot$rate_trans<-log(xdata_table_3states_plot$rate_trans+10^(-minimum_exp),10)
xdata_table_3states_plot$rate<-factor(xdata_table_3states_plot$rate,levels=orderforplotting)
xdata_table_3states_plot<-xdata_table_3states_plot %>% arrange(rate)
xdata_table_3states_plot$rate<-fct_inorder(as.factor(as.character(xdata_table_3states_plot$rate)))
mypalette_t<-mypalette[orderforplotting %in% (xdata_table_3states_plot$rate %>% unique)]
induction_params<-c("Error","u (uncut fraction)","r (induction speed)","d (induction decay)")
xdata_table_3states_plot$target<-fct_inorder(as.factor(as.character(xdata_table_3states_plot$target)))
#Choose whether to plot also induction parameters or not
if (plotalsoinduction) {
mypalette_t<-mypalette[orderforplotting %in% (xdata_table_3states_plot$rate %>% unique) ] 
xdata.all.temp<-xdata_table_3states_plot 
} else 
{
mypalette_t<-mypalette[orderforplotting %in% (xdata_table_3states_plot$rate %>% unique) & !orderforplotting %in% induction_params]
xdata.all.temp<-xdata_table_3states_plot %>% filter(!rate %in% induction_params)
}
xdata.all.temp$induction[xdata.all.temp$induction=="estimate_fi"]<-"FACS"
xdata.all.temp$induction[xdata.all.temp$induction=="estimate"]<-"co-estimated"
xcolours<-c("darkseagreen","darkorange")
#Plot split violin plots using github library introdataviz
library(introdataviz)
p<-ggplot(xdata.all.temp, aes(x=rate, y=rate_trans, fill=induction)) +  
introdataviz::geom_split_violin(trim=TRUE, scale="width")+
scale_fill_manual(values = xcolours, name = "Induction curve")+
facet_grid(experiment~target) + 
theme(strip.background = element_rect(fill = "white"),axis.text.x = element_text(angle = 45, hjust=1), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.03),legend.position = "bottom",aspect.ratio=1) + ylab("rate ( event/hour )") +
scale_y_continuous(breaks=xbreaks,labels=xlabels,limits=c(min(xbreaks)-0.01,log(20,10)))
ggsave(p,file=paste0("comparisonFACS_",xmodel,"_",myexp,"_",plotalsoinduction,".pdf"))
}}}



