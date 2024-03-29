#!/usr/bin/env Rscript
#example to run
#./optimize_model_backbone.R  -t CRTISO -i RNP -m modelDSBs1i1_fullimpreciseDSB -e E_errorsfromunbroken -n 10000 -o CRTISO_RNP_optimization.tsv

library(argparser, quietly=TRUE,warn.conflicts=FALSE)
library(tidyverse)
library(optimx)
library(ggplot2)
library(deSolve)

source("likelihood_functions.R")

################################################################################
################## READ ARGUMENTS ##############################################
################################################################################
argv<- arg_parser("Parse arguments")
argv <- add_argument(argv, "-i", help="input_files; a text file specifying in each row an RData file, output from CI.R, for each of the bootstrap estimates")
argv <- add_argument(argv, "-d", help="data_file; time course used to calculate likelihood with optimization.R")
argv <- add_argument(argv, "-o", help="root for output_files")
argv <- add_argument(argv, "-E", help="error matrix")
argv <- add_argument(argv, "-m", help="model")
argv <- add_argument(argv, "-l", help="switch to change likelihood function. To estimate a common error from DSB select 1. Default is no estimate from data, only from controls, to sample 0.2h after induction (0). To set induction curve without delay select 2. To model imprecise DSB as misread precise DSB select 3", default=0)
argv <- add_argument(argv, "-y", help="condition flow to induction curve, i.e. calculate the flow in ideal system in which full induction", default=2)
argv <- add_argument(argv, "-z", help="n parameters in induction curve", default=2)
argv <- add_argument(argv, "-r", help="time resolution; default is 1/100 hours, represented as 0.01", default=0.01)
argv <- add_argument(argv, "-T", help="final time in hours; default is 72h",default=72)
argv <- add_argument(argv, "-w", help="calculate_flow. Default is 0, which indicates FALSE. Select 1 to calculate it.", default=0)
argv <- add_argument(argv, "-n", help="maximum number of curves used for CI, when using likelihood based CI..")
argv <- add_argument(argv, "-j", help="plot processed. choose 0 for no.", default=1)

args <- parse_args(argv)
input_file<-args$i
output_file<-args$o
myerrorE<-args$E
data_file<-args$d
mymodel<-args$m
calculate_flow<-args$w
timeresolution<-as.numeric(args$r)
maxnsims<-as.numeric(args$n)
optimize_errorDSB2indel<-as.numeric(args$l)
normalize.byind<-as.numeric(args$y)
nparamsind<-as.numeric(args$z)
noplotprocessedDSB<-as.numeric(args$j)
finalTime<-as.numeric(args$T)

if (is.na(nparamsind)){ nparamsind<-2 }
if (is.na(calculate_flow)){ calculate_flow<-0 }
if ( calculate_flow!=1){ calculate_flow<-FALSE } else { calculate_flow<-TRUE }
if (is.na(noplotprocessedDSB)){ noplotprocessedDSB<-1 }
if (is.na(optimize_errorDSB2indel)) { optimize_errorDSB2indel<-0 }
if (is.na(normalize.byind)) { normalize.byind<-0 }
if (is.na(timeresolution)) { timeresolution<-0.01 }
if (is.na(finalTime)) { finalTime<-72 }
define_ODE.functions(nparamsind)

#input_file<-"/home/labs/alevy/fabrizio/workspace/daniela/resultsv4/stratifiedbootstraps/timecourse_RNP_Psy1_allb/modelDSBs1i1_mini/results_modelDSBs1i1_mini_RNP_ind.c2_Psy1_allb/listCIfiles.txt"
#data_file<-"/home/labs/alevy/fabrizio/workspace/daniela/input_datasets/timecourse_RNP_Psy1_allb.txt"
#myerrorE<-"/home/labs/alevy/fabrizio/workspace/daniela/error_matrices/error_matrix4_Psy1_errorsfromunbroken.tsv"
#mymodel<-"modelDSBs1i1_mini"

errormatrix<-as.matrix(read.table(myerrorE, header=FALSE))
nameparms<-model2nameparams(mymodel,nparamsind)

################################################################################
################# ADDITIONAL FUNCTIONS FOR CI ##################################
################################################################################

loglik_er_f<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix){
  #ODEfunc: a function of class desolve
  #my_data: a dataframe with "time" course (which should always include 0) as 1st column, and col-types compatible with ODEfunc
  #parms: vector specifying the set of parameters (compatible with ODEfunc)
  #ncells<-apply(mydata[,-1],MARGIN=1,sum)
  yinit<-c(1,rep(0,ncol(my_data)-2))
  out <- ode (times = my_data$time, y = yinit, func = ODEfunc, parms = parms)
  probs<-(out[,2:ncol(out)]) %*% E.matrix
  probs[probs<10^(-16)]<-10^(-16)
  loglik<-sapply(1:nrow(my_data), function(x) dmultinom(x=my_data[x,-1],prob=probs[x,],log=TRUE))
  sum(loglik)
}

################################################################################
################################# PARSE MODEL ##################################
################################################################################
mydata<-read.table(data_file,header=TRUE)
time_courses_begins<-c(which(mydata$time[-1]-mydata$time[-length(mydata$time)]<0)+1)
mydelay<-0
ntypes<-4
xmodel<-0
errormatrix<-as.matrix(read.table(myerrorE, header=FALSE))
define_ODE.functions(nparamsind)
loglik_er_f.pen<-model2likelihoodfunction(mymodel,optimize_errorDSB2indel)
ntypes<-model2ntypes(mymodel)
xmodel<-model2xmodel(mymodel)
ncurves<-model2ncurves(mymodel)
nparms<-length(nameparms) #for how it was built, this (if er1 present) should have 1 less param than nrates - important for 3x4
nheaders<-4
nparms_induction<-nparamsind
nrates<-length(nameparms)
nparamserr<-length(grep("^er",nameparms))
if (optimize_errorDSB2indel==2) { mydelay<-0 }
mypalette_flow=c("Cutting"="firebrick","Processing"="darkorange1","Repair-error (from direct DSB)"="deepskyblue2","Repair-error (from processed DSB)"="dodgerblue4",
"Precise Repair (from direct DSB)"="chartreuse3","Precise Repair (from processed DSB)"="darkgreen","Error"="gray65","u (uncut fraction)"="burlywood4",
"r (induction speed)"="orchid3","d (induction decay)"="purple3")

mypalette_t<-c("intact"="chartreuse3", "preciseDSB"="firebrick","impreciseDSB"="darkorange1","indels"="deepskyblue2")
########################################################################################################
##################### LOAD DATA ########################################################################
########################################################################################################

mydata<-read.table(data_file,header=TRUE)
time_courses_begins<-c(which(mydata$time[-1]-mydata$time[-length(mydata$time)]<0)+1)
print(data_file)
print(mydata)
mydelay<-0

if (ntypes==4)
	{
	mydata.p<-as.data.frame(cbind(mydata$time,t(apply(mydata[,2:5],MARGIN=1,FUN=function(x) x/sum(x)))))
	names(mydata.p)<-c("time","y1","y2","y3","y4")
	mydata0.p<-mydata.p
	names(mydata0.p)<-c("time","intact","indels","preciseDSB","impreciseDSB")
	tidy.mydata0.p<-mydata0.p %>% pivot_longer(names(mydata0.p)[2:5],names_to = "types", values_to = "p")
	if (noplotprocessedDSB!=1)
		{
		mydata0.p<-mydata0.p %>% mutate(preciseDSB=preciseDSB+impreciseDSB,impreciseDSB=NULL)
		tidy.mydata0.p<-mydata0.p %>% pivot_longer(names(mydata0.p)[2:4],names_to = "types", values_to = "p")
		}
	} else if (ntypes==3)
	{
	mydata.p<-as.data.frame(cbind(mydata$time,t(apply(mydata[,2:4],MARGIN=1,FUN=function(x) x/sum(x)))))
	names(mydata.p)<-c("time","y1","y2","y3")
	mydata0.p<-mydata.p
	names(mydata0.p)<-c("time","intact","DSB","indels")
	tidy.mydata0.p<-mydata0.p %>% pivot_longer(names(mydata0.p)[2:4],names_to = "types", values_to = "p")
	}

print(tidy.mydata0.p)
errormatrix<-as.matrix(read.table(myerrorE, header=FALSE))
if ( noplotprocessedDSB == 0 && ntypes==4) { errormatrix[,4]<-0 }

################################################################################
############ PARSE BOOTSTRAP AND EXTRACT BOOTSTRAPPED CURVES ###################
################################################################################

#input_file="/home/labs/alevy/fabrizio/workspace/daniela/resultsv4/stratifiedbootstraps/timecourse_RNP_Psy1_allb/modelDSBs1i1_mini/results_modelDSBs1i1_mini_RNP_ind.c2_Psy1_allb/listCIfiles.txt"
#ls /home/labs/alevy/fabrizio/workspace/daniela/resultsv4/results_modelDSBs1i1_mini_RNP_ind.c2_Psy1_allb.all_n500000.CI > /home/labs/alevy/fabrizio/workspace/daniela/resultsv4/stratifiedbootstraps/timecourse_RNP_Psy1_allb/modelDSBs1i1_mini/results_modelDSBs1i1_mini_RNP_ind.c2_Psy1_allb/listCIfiles.txt
#ls /home/labs/alevy/fabrizio/workspace/daniela/resultsv4/stratifiedbootstraps/timecourse_RNP_Psy1_allb/modelDSBs1i1_mini/results_modelDSBs1i1_mini_RNP_ind.c2_Psy1_allb/results*.tab_n50000.CI >> /home/labs/alevy/fabrizio/workspace/daniela/resultsv4/stratifiedbootstraps/timecourse_RNP_Psy1_allb/modelDSBs1i1_mini/results_modelDSBs1i1_mini_RNP_ind.c2_Psy1_allb/listCIfiles.txt
if (mymodel %in% c("modelDSBs1i1_3x4","modelDSBs1i1_3x4.bytarget")){ xmodel="modelDSBs1i1_3x4" }
if (mymodel %in% c("modelDSBs1i1_3x4nor11")){ xmodel="modelDSBs1i1_3x4nor11" }


#PARSE INPUT FILE AND DECIDE FROM THAT IF LIST OF CI FILES (THEN BOOTSTRAP) OR SINGLE CI.RData FILE
print("PARSE INPUT FILE")
if(length(grep(".RData",input_file)))
{
  print("input_file is a .RData file. Plotting CI from this file.")
  load(input_file)
  simsinCI<-res[[2]][res[[2]]$ok==1,];
  if (sum(res[[2]]$ok==1)<2)
    {
    print("WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("UNRELIABLE CI!!!")
    simsinCI<-res[[2]][998:1000,];
    }
  nsims<-nrow(simsinCI)
  simsinCI<-simsinCI[order(simsinCI$loglik),]#[1:nsims,]
  nsims<-min(maxnsims,nsims)
  xsims<-sample(1:nsims)
  xsims[1]<-2 #commented out 22/11/2023 to avoid crashing when no sims in CI
  #FIT TRAJECTORY FROM ACTUAL DATA#
  restemp<-res[[1]]
  #parsing to avoid duplicate parameters
  while( sum(duplicated(restemp$rate))>0)
	{
  	dupl.rates<-unique(restemp$rate[duplicated(restemp$rate)])
        for ( i.dup in dupl.rates)
	    {
	    n.dup<-sum(restemp$rate==i.dup)
	    dup.ending<-c("",paste0(".",1:(n.dup-1)))
	    restemp$rate[restemp$rate==i.dup]<-paste0(restemp$rate[restemp$rate==i.dup],dup.ending)
	    }
	}
  bestmodels_t<- restemp %>% select(rate,max) %>% pivot_wider(names_from=rate,values_from=max) %>% ungroup
  print(bestmodels_t)
  } else
  {
  print("input_file is not a .RData file. Using it as a list of CI files to extract curves.")
  listfiles<-unlist(read.table(input_file)[,1])
  for (iit in 1:length(listfiles)){
    print(iit)
    res<-read.table(listfiles[iit],sep="\t",header=TRUE)
    restemp<-res  %>% select(rate,max) %>% pivot_wider(names_from=rate,values_from=max)
    if (iit==1) { allmaxs<-restemp } else { allmaxs<-rbind(allmaxs,restemp) }
    }
    nsims<-nrow(allmaxs)
    simsinCI<-allmaxs
    xsims<-2:nsims
    #FIT TRAJECTORY FROM ACTUAL DATA#
    i<-1; bestmodels_t<-allmaxs[i,]
  }


bestmodels.fitted.0er<-predict_models(bestmodels_t,ntypes=ntypes,nparms=nrates,nheaders=0,errormatrix=diag(ntypes),mymodel=xmodel,timest=seq(0,finalTime,timeresolution),noplotprocessedDSB=noplotprocessedDSB)
########FIT TRAJECTORIES AND INDUCTION FROM BOOTSTRAPS #####################
print("FIT TRAJECTORIES AND INDUCTION FROM BOOTSTRAPS")
names(simsinCI)[1:length(bestmodels_t)]<-names(bestmodels_t)
simsinCIall<-simsinCI
x.dupls<-unique(sub(".*\\.",".",names(bestmodels_t)))
x.dupls<-x.dupls[grep("\\.",x.dupls)]

for (x.dupl in c(".0",x.dupls))
  {
  if (x.dupl!=".0")
  	{
	newcols<-names(simsinCIall)[grep(paste0("\\",x.dupl),names(simsinCIall))]
	oldcols<-newcols %>% sapply(function(x) sub(paste0("\\",x.dupl),"",x)) %>% unname
	simsinCI<-simsinCIall[,-which(names(simsinCIall) %in% oldcols)]
	names(simsinCI)<-names(simsinCI) %>% sapply(function(x) sub(paste0("\\",x.dupl),"",x)) %>% unname
	} else { simsinCI<-simsinCIall }
    for ( isim in xsims)
        {
        if ( (isim %% 10)==0){print(isim)};
	bestmodels_t.cfitted.sims<-predict_induction(simsinCI[isim,],nparms=nparms-nparms_induction-nparamserr,nheaders=0,nparms_induction=nparms_induction,induction_function=induction_curve_vectorized_default)
        induction_for_normalized_k11<-bestmodels_t.cfitted.sims %>% filter(time==6) %>% select(curve_fitted) %>% as.numeric
        if (normalize.byind==1) { k11<-simsinCI[isim,1]/induction_for_normalized_k11 } else { k11<-as.numeric(simsinCI[isim,1]) }
	if ("k12" %in% nameparms)
        	{
                if (normalize.byind==1){ k12<-simsinCI[isim,2]/induction_for_normalized_k11 } else { k12<-as.numeric(simsinCI[isim,2]) }
		}
	if (normalize.byind==1) 
		{
	        if ("k12" %in% nameparms)
        	        {
	                bestmodels_t.cfitted.sims$curve_fitted<-bestmodels_t.cfitted.sims$curve_fitted*(k11+k12)
        	        } else
                	{
	                bestmodels_t.cfitted.sims$curve_fitted<-bestmodels_t.cfitted.sims$curve_fitted*k11
        	        }
		}
        if (isim==2) { bestmodels.cfitted.CI<-bestmodels_t.cfitted.sims } else { bestmodels.cfitted.CI<-rbind(bestmodels.cfitted.CI,bestmodels_t.cfitted.sims) } #changed on 22/11/2023 from 2
        simsinCI[isim,1]<-k11
        if ("k12" %in% nameparms) {simsinCI[isim,2]<-k12}
        bestmodels_t.fitted.sims<-predict_models(simsinCI[isim,],ntypes=ntypes,nparms=nrates,nheaders=0,errormatrix=errormatrix,mymodel=xmodel,noplotprocessedDSB=noplotprocessedDSB) #"E_errorsfromunbroken"
        if (isim==2) { bestmodels.fitted.CI<-bestmodels_t.fitted.sims } else { bestmodels.fitted.CI<-rbind(bestmodels.fitted.CI,bestmodels_t.fitted.sims)} #changed on 22/11/2023 from 2
        }
bestmodels.fitted.CI<-bestmodels.fitted.CI %>% group_by(time,types) %>% summarise(lowCI=min(p),highCI=max(p)) %>% ungroup %>% mutate(replicate=x.dupl)
bestmodels.cfitted.CI<-bestmodels.cfitted.CI %>% group_by(time) %>% summarise(lowCI=min(curve_fitted),highCI=max(curve_fitted)) %>% ungroup %>% mutate(replicate=x.dupl)
  if (x.dupl==".0")
	{
	bestmodels.fitted.CI.tot<-bestmodels.fitted.CI 
	bestmodels.cfitted.CI.tot<-bestmodels.cfitted.CI
	} else
	{
	bestmodels.fitted.CI.tot<-rbind(bestmodels.fitted.CI.tot,bestmodels.fitted.CI)
        bestmodels.cfitted.CI.tot<-rbind(bestmodels.cfitted.CI.tot,bestmodels.cfitted.CI)
	}
  }
bestmodels.fitted.CI<-bestmodels.fitted.CI.tot
bestmodels.cfitted.CI<-bestmodels.cfitted.CI.tot

#####################3 CALCULATE MEAN INDUCTION ###############
print("CALCULATE MEAN INDUCTION")
bestmodels_t0<-bestmodels_t
for (x.dupl in c(".0",x.dupls))
  {
  if (x.dupl!=".0")
        {
        newcols<-names(bestmodels_t0)[grep(paste0("\\",x.dupl),names(bestmodels_t0))]
        oldcols<-newcols %>% sapply(function(x) sub(paste0("\\",x.dupl),"",x)) %>% unname
        bestmodels_t<-bestmodels_t0[,-which(names(simsinCIall) %in% oldcols)]
        names(bestmodels_t)<-names(bestmodels_t) %>% sapply(function(x) sub(paste0("\\",x.dupl),"",x)) %>% unname
        }
  bestmodels.cfitted<-predict_induction(bestmodels_t,nparms=nparms-nparms_induction-nparamserr,nheaders=0,nparms_induction=nparms_induction,induction_function=induction_curve_vectorized_default)
  #extract induction at time 6h to calculate k11-no-induction
  induction_for_normalized_k11<-bestmodels.cfitted %>% filter(time==6) %>% select(curve_fitted) %>% as.numeric
  if (normalize.byind==1)
        {
        k11<-bestmodels_t$k11/induction_for_normalized_k11
        bestmodels_t$k11<-k11
        if ("k12" %in% nameparms)
                {
                k12<-bestmodels_t$k12/induction_for_normalized_k11
                bestmodels_t$k12<-k12
                }
     	if ("k12" %in% nameparms)
        	{
	        #rescale mean induction in terms of induction*k11 (cutting flow)
        	bestmodels.cfitted$curve_fitted<-bestmodels.cfitted$curve_fitted*(k11+k12)
	        } else
        	{
	        bestmodels.cfitted$curve_fitted<-bestmodels.cfitted$curve_fitted*(k11)
        	}
	}
  bestmodels.cfitted <- bestmodels.cfitted %>% mutate(replicate=x.dupl)
  #calculate mean trajectory with k11-no-induction
  bestmodels.fitted<-predict_models(bestmodels_t,ntypes=ntypes,nparms=nrates,nheaders=0,errormatrix=errormatrix,mymodel=xmodel,noplotprocessedDSB=noplotprocessedDSB)
  bestmodels.fitted <- bestmodels.fitted %>% mutate(replicate=x.dupl)
  #calculate mean trajectory with no errors
  bestmodels_t00<-bestmodels_t
#---SET THE VALUES AS BELOW IF FLOW WITH FULL INDUCTION IS DESIRED:---
#  bestmodels_t00[which(names(bestmodels_t0)=="r0")]<-999999999
#  bestmodels_t00[which(names(bestmodels_t0)=="r2")]<-0
#  bestmodels_t00[which(names(bestmodels_t0)=="x0")]<-0.0000001
#---------------------------------------------------------------------
  if (mymodel %in% c("modelDSBs1i1_3x4","modelDSBs1i1_3x4nor11","modelDSBs1i1_3x4.bytarget")){ bestmodels_t00[names(bestmodels_t00)=="er1"]<-0 }
  print(bestmodels_t00)
  print(xmodel)
  bestmodels.fitted.0er<-predict_models(bestmodels_t00,ntypes=ntypes,nparms=nrates,nheaders=0,errormatrix=diag(ntypes),mymodel=xmodel,timest=seq(0,finalTime,timeresolution),noplotprocessedDSB=noplotprocessedDSB)
  bestmodels.fitted.0er %>% mutate(replicate=x.dupl)
  if (x.dupl==".0")
	{
	bestmodels.fitted.tot<-bestmodels.fitted 
	bestmodels.fitted.0er.tot<-bestmodels.fitted.0er
	bestmodels.cfitted.tot<-bestmodels.cfitted
	} else
	{
	bestmodels.fitted.tot<-rbind(bestmodels.fitted.tot,bestmodels.fitted)
        bestmodels.fitted.0er.tot<-rbind(bestmodels.fitted.0er.tot,bestmodels.fitted.0er)
        bestmodels.cfitted.tot<-rbind(bestmodels.cfitted.tot,bestmodels.cfitted)
	}
  }
  bestmodels.cfitted<-left_join(bestmodels.cfitted.tot,bestmodels.cfitted.CI)
  bestmodels.fitted <- bestmodels.fitted.tot %>% left_join(bestmodels.fitted.CI)

#output_file="./"
print("PLOT")
print(tidy.mydata0.p)
if (finalTime==24){ xbreaks=c(0,2,4,6,12,24) } else if (finalTime==72) { xbreaks=c(0,6,12,24,36,48,72) }
#====PLOT TRAJECTORIES====
if (length(unique(bestmodels.fitted$replicate))>1)
{
myplot <- bestmodels.fitted %>% ggplot(aes(x=time,y=p,colour=types)) +
xlim(0,finalTime+0.2) + scale_y_continuous("p", breaks=c(0,0.01,0.1,0.25,0.5,0.75,1), trans='sqrt') + scale_x_continuous("time (hours)", breaks=xbreaks, limits=c(0,finalTime+0.5))+
scale_color_manual(values=mypalette_t)+scale_fill_manual(values=mypalette_t)+
geom_ribbon(aes(ymin = lowCI, ymax = highCI,fill=types), alpha = 0.2,outline.type="both",linetype="blank")+ #,linetype="1F")+
geom_line()+geom_point(data=tidy.mydata0.p)+
theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),strip.background = element_rect(fill = "white"), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.01),legend.position = "bottom",aspect.ratio=1,text = element_text(size = 16))+
#theme_bw()+theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank())
facet_wrap(~replicate)
} else
{
print("running this")
myplot <- bestmodels.fitted %>% ggplot(aes(x=time,y=p,colour=types)) +
xlim(0,finalTime+0.2) + 
scale_y_continuous("p", breaks=c(0,0.01,0.1,0.25,0.5,0.75,1), trans='sqrt') + 
scale_x_continuous("time (hours)", breaks=xbreaks, limits=c(0,finalTime+0.5))+
scale_color_manual(values=mypalette_t)+scale_fill_manual(values=mypalette_t)+
geom_ribbon(aes(ymin = lowCI, ymax = highCI,fill=types), alpha = 0.2,outline.type="both",linetype="blank")+
geom_line()+geom_point(data=tidy.mydata0.p)+
theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),strip.background = element_rect(fill = "white"), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.01),legend.position = "bottom",aspect.ratio=1,text = element_text(size = 24))
#theme_bw()+theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank())
}

ggsave(myplot,filename=paste0(output_file,"_plot.trajectories.pdf"))

#====PLOT INDUCTION CURVES====
#datatmpnor11<-bestmodels.cfitted %>% filter(target==mytarget,induction==myinduction,model=="modelDSBs1i1_realimpnor11")
myplot <- bestmodels.cfitted %>% ggplot(aes(x=time,y=curve_fitted)) +
xlim(0,finalTime+0.2) + 
scale_y_continuous("precise cuts/intact DNA per hour", breaks=c(0,0.01,0.1,0.25,0.5,0.75,1), trans='sqrt') + 
scale_x_continuous("time (hours)", breaks=xbreaks, limits=c(0,finalTime+0.5)) +
geom_ribbon(aes(ymin = lowCI, ymax = highCI), alpha=0.1,linetype="dotted")+
geom_line()+
theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),strip.background = element_rect(fill = "white"), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.01),legend.position = "bottom",aspect.ratio=1,text = element_text(size = 24))
#theme(strip.background = element_rect(fill = "white"), panel.background = element_rect(fill = 'white',color='black'),panel.grid.major = element_line(color = 'grey', linetype = 'longdash',size=0.01),legend.position = "bottom",aspect.ratio=1,text = element_text(size = 24))
ggsave(myplot,filename=paste0(output_file,"_plot.induction.pdf"))

#====CALCULATE FLOW====
if (calculate_flow && noplotprocessedDSB==1){
print("calculate flow")
print(bestmodels_t)
y.fitted<-bestmodels.fitted.0er %>% pivot_wider(names_from="types",values_from="p",values_fill=0)
if (ntypes==4){
names(y.fitted)<-c("time","y1","y2","y3","y4")
y.fitted <- y.fitted %>% group_by(time) %>% summarise(y1=sum(y1),y2=sum(y2),y3=sum(y3),y4=sum(y4)) #%>% select(y1,y2,y3,y4)
#print(y.fitted)
} else if (ntypes==3){
names(y.fitted)<-c("time","y1","y2","y3")
y.fitted <- y.fitted %>% group_by(time) %>% summarise(y1=sum(y1),y2=sum(y2),y3=sum(y3))
}
bestmodels_t0<-bestmodels_t00
for (ipar in 1:(nparms-nparms_induction-nparamserr)){bestmodels_t0[ipar]<-0}

flow_l<-list()
#substitute param to obtain a full induction so that 1.
bestmodels_t0[which(names(bestmodels_t0)=="r0")]<-999999999999
bestmodels_t0[which(names(bestmodels_t0)=="r2")]<-0
bestmodels_t0[which(names(bestmodels_t0)=="K")]<-1
if (nparamsind==3){ bestmodels_t0[which(names(bestmodels_t0)=="K")]<-0 }
bestmodels.cfitted<-predict_induction(bestmodels_t,nparms=nparms-nparms_induction-nparamserr,nheaders=0,nparms_induction=nparms_induction,induction_function=induction_curve_vectorized,times=seq(0,finalTime,timeresolution),yno1=1-y.fitted$y1)
#rescale k11 to have induction inside
for (ipar in 1:(nparms-nparms_induction-nparamserr)){
bestmodels_tt<-bestmodels_t0
bestmodels_tt[ipar]<-bestmodels_t[ipar]
bestmodels_tt.df<-sapply(1:length(seq(0,finalTime,timeresolution)), function(x) unlist(bestmodels_tt)) %>% t %>% as.data.frame
if (names(bestmodels_t0)[ipar]=="k11" || names(bestmodels_t0)[ipar]=="k12" )
{
#bestmodels_tt.df[,ipar]<-bestmodels_tt.df[,ipar]*as.numeric(bestmodels.cfitted$curve_fitted)
flow_l[[nameparms[ipar]]]<-bestmodels_t[ipar]*sum(as.numeric(bestmodels.cfitted$curve_fitted)*timeresolution)
} else 
{
changess<-lapply(1:nrow(y.fitted), function(x) predict_models(bestmodels_tt.df[x,],ntypes=ntypes,nparms=nrates,nheaders=0,errormatrix=diag(ntypes),mymodel=xmodel,yinit=unlist(y.fitted[x,-1]),timest=c(0,timeresolution)))

#flow_l[[nameparms[ipar]]]<-sum(sapply(changes,function(x) sum(abs(x$p[x$time==0.1]-x$p[x$time==0]))/2))

flows<-sapply(changess,function(x) max(abs(x$p[x$time==timeresolution]-x$p[x$time==0])))
flow_l[[nameparms[ipar]]]<-sum(flows)
}

#abundance_l<-apply(y.fitted[,-1], MARGIN=2,FUN=function(x) mean(x))
#changes_l[[nameparms[ipar]]]<-sum(sapply(changes,function(x) x$p[x$time==0.0002]-x$p[x$time==0]))
#print("changes")
#print(head(changes))
}
#print("abundance")
#print(unlist(abundance_l))
#print("flow")
#print(unlist(flow_l))
#write.table(unlist(flow_l),file=paste0(output_file,"_plot.flow.tab"))
AIC<-0
myerror<-try( AIC<-2*(ncol(res[[2]])-4)-2*res[[2]]$maxll[1] )
print(myerror)
AICtable<-data.frame(AICt="AIC",AICv=AIC)
write.table(AICtable,file=paste0(output_file,"_plot.flow.tab"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(unlist(flow_l),file=paste0(output_file,"_plot.flow.tab"),quote=FALSE,sep="\t",append=TRUE,col.names=FALSE)

#calculate flow for extreme values:
#for ( i in 1:length(extremes_sims)){
#restemp<-simsinCI[extremes_sims[i],]
#worsemodel_t<- restemp %>% select(rate,max) %>% pivot_wider(names_from=rate,values_from=max) %>% ungroup
#bestmodels.fitted.0er<-predict_models(worsemodel_t,ntypes=ntypes,nparms=nrates,nheaders=0,errormatrix=diag(ntypes),mymodel=mymodel,timest=seq(0,72,timeresolution))
#}
}

