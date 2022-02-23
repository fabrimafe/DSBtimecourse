#!/usr/bin/env Rscript
#example to run
#./optimize_model_backbone.R  -t CRTISO -i RNP -m modelDSBs1i1_fullimpreciseDSB -e E_errorsfromunbroken -n 10000 -o CRTISO_RNP_optimization.tsv

library(argparser, quietly=TRUE,warn.conflicts=FALSE)
library(tidyverse)
library(optimx)
library(ggplot2)
library(deSolve)

################################################################
# DRAW THE MODEL (an ODE function) AND GENERATE SOME RANDOM DATA
################################################################

source("likelihood_functions.R")

################################################################################
################## READ ARGUMENTS ##############################################
################################################################################


argv<- arg_parser("Parse arguments")

argv <- add_argument(argv, "-i", help="input_file; RData file, output from CI.R")
argv <- add_argument(argv, "-d", help="data_file; time course used to calculate likelihood with optimization.R")
argv <- add_argument(argv, "-o", help="roof for output_files")
argv <- add_argument(argv, "-E", help="error matrix")
argv <- add_argument(argv, "-m", help="model")
argv <- add_argument(argv, "-n", help="max accuracy for CI. Recommended >=100. Beyond 1000 calculations can be very slow.")
argv <- add_argument(argv, "-l", help="switch to change likelihood function. To estimate a common error from DSB select 1. Default is no estimate from data, only from controls, to sample 0.2h after induction (0). To set induction curve without delay select 2. To model imprecise DSB as misread precise DSB select 3", default=0)
argv <- add_argument(argv, "-y", help="condition flow to induction curve, i.e. calculate the flow in ideal system in which full induction", default=2)
argv <- add_argument(argv, "-z", help="n parameters in induction curve", default=2)
argv <- add_argument(argv, "-r", help="time resolution; default is 0.01")


args <- parse_args(argv)

input_file<-args$i
output_file<-args$o
myerrorE<-args$E
data_file<-args$d
mymodel<-args$m
timeresolution<-as.numeric(args$r)
maxnsims<-as.numeric(args$n)
optimize_errorDSB2indel<-as.numeric(args$l)
normalize.byind<-as.numeric(args$y)
nparamsind<-as.numeric(args$z)
if (is.na(nparamsind)){ nparamsind<-2 }
define_ODE.functions(nparamsind)


load(input_file)
if (is.na(optimize_errorDSB2indel)) { optimize_errorDSB2indel<-0 }
if (is.na(normalize.byind)) { normalize.byind<-0 }
if (is.na(timeresolution)) { timeresolution<-0.01 }



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
#mydata<-read.table(input.file, header=TRUE)
errormatrix<-as.matrix(read.table(myerrorE, header=FALSE))
nameparms<-model2nameparams(mymodel,nparamsind)

if (mymodel=="model5i1" || mymodel=="model5i1_nor11")
                {
                ntypes<-3
                };
if ( mymodel=="modelDSBs1i1_realimprecise.inductionx3")
                {
                loglik_er_f.pen<-loglik_er_f.pen_modelinductionx3
                xmodel<-get("modelDSBs1i1_realimprecise")
                };
nparms<-length(nameparms)
nheaders<-4
nparms_induction<-nparamsind

if ( optimize_errorDSB2indel==1 )
        {
        loglik_er_f.pen<-loglik_er_f.pen_errorDBS2indel_4states_m1
        nameparms<-c(nameparms,"er1");
        print("error model is 1")
        } else
if ( optimize_errorDSB2indel==2 )
        {
        loglik_er_f.pen<-loglik_er_f.nopen
        print("error model is 2")
        }  else
if ( optimize_errorDSB2indel==3 )
        {
        loglik_er_f.pen<-loglik_er_f.pen_errorimpDSB2pDSB_4states_m1
	nameparms<-c(nameparms,"er1")
        };
nrates<-length(nameparms)

if (optimize_errorDSB2indel==2) { mydelay<-0 }

if (!"function" %in% is(xmodel)) { xmodel<-get(mymodel) }
########################################################################################################
##################### LOAD DATA ########################################################################
########################################################################################################

mydata<-read.table(data_file,header=TRUE)
time_courses_begins<-c(which(mydata$time[-1]-mydata$time[-length(mydata$time)]<0)+1)
mydelay<-0

if (ntypes==4){
mydata.p<-as.data.frame(cbind(mydata$time,t(apply(mydata[,2:5],MARGIN=1,FUN=function(x) x/sum(x)))))
names(mydata.p)<-c("time","y1","y2","y3","y4")
mydata0.p<-mydata.p
names(mydata0.p)<-c("time","intact","indels","preciseDSB","impreciseDSB")
tidy.mydata0.p<-mydata0.p %>% pivot_longer(names(mydata0.p)[2:5],names_to = "types", values_to = "p")
} else if (ntypes==3){
mydata.p<-as.data.frame(cbind(mydata$time,t(apply(mydata[,2:4],MARGIN=1,FUN=function(x) x/sum(x)))))
names(mydata.p)<-c("time","y1","y2","y3")
mydata0.p<-mydata.p
names(mydata0.p)<-c("time","intact","DSB","indels")
tidy.mydata0.p<-mydata0.p %>% pivot_longer(names(mydata0.p)[2:4],names_to = "types", values_to = "p")
}        

#mydata<-read.table(input.file, header=TRUE)
errormatrix<-as.matrix(read.table(myerrorE, header=FALSE))
#if (mymodel=="model5i1" || mymodel=="model5i1_nor11")
#        {
#        ntypes<-3 ; load("RData/error_matrices3_l.RData");
#        error_mat_l<-error_matrices3_l
#        }
#
if (optimize_errorDSB2indel==2) { mydelay<-0 } 

#################################################################################
##################### EXTRACT CURVES FOR CI ########################################
#################################################################################
#nsims<-100
simsinCI<-res[[2]][res[[2]]$ok==1,]; 
nsims<-nrow(simsinCI)
simsinCI<-simsinCI[order(simsinCI$loglik),]#[1:nsims,] 
nsims<-min(maxnsims,nsims)
xsims<-sample(1:nsims)
for ( isim in 1:nsims)
        {
	if ( (isim %% 100)==0){print(isim)};
        bestmodels_t.cfitted.sims<-predict_induction(simsinCI[xsims[isim],],nparms=nparms-nparms_induction,nheaders=0,nparms_induction=nparms_induction,induction_function=induction_curve_vectorized)
	induction_for_normalized_k11<-bestmodels_t.cfitted.sims %>% filter(time==6) %>% select(curve_fitted) %>% as.numeric
	if (normalize.byind==1){ k11<-simsinCI[xsims[isim],1]/induction_for_normalized_k11 } else { k11<-simsinCI[xsims[isim],1] }
	if ("k12" %in% nameparms)
		{
		if (normalize.byind==1){ k12<-simsinCI[xsims[isim],2]/induction_for_normalized_k11 } else { k12<-simsinCI[xsims[isim],2] }
		bestmodels_t.cfitted.sims$curve_fitted<-bestmodels_t.cfitted.sims$curve_fitted*(k11+k12)
		} else 
		{
		bestmodels_t.cfitted.sims$curve_fitted<-bestmodels_t.cfitted.sims$curve_fitted*(k11)
		}
	if (isim==1) { bestmodels.cfitted.CI<-bestmodels_t.cfitted.sims } else { bestmodels.cfitted.CI<-rbind(bestmodels.cfitted.CI,bestmodels_t.cfitted.sims) }
        simsinCI[xsims[isim],1]<-k11
	if ("k12" %in% nameparms) {simsinCI[xsims[isim],2]<-k12}
	bestmodels_t.fitted.sims<-predict_models(simsinCI[xsims[isim],],ntypes=ntypes,nparms=nrates,nheaders=0,errormatrix=errormatrix,mymodel=mymodel) #"E_errorsfromunbroken"
	if (isim==1) { bestmodels.fitted.CI<-bestmodels_t.fitted.sims } else { bestmodels.fitted.CI<-rbind(bestmodels.fitted.CI,bestmodels_t.fitted.sims) }
	#print(simsinCI[xsims[isim],])
	#print(bestmodels_t.fitted.sims %>% filter(time==0))
	#print(bestmodels_t.fitted.sims %>% filter(time==72))
        }
bestmodels.fitted.CI<-bestmodels.fitted.CI %>% group_by(time,types) %>% summarise(lowCI=min(p),highCI=max(p)) %>% ungroup
bestmodels.cfitted.CI<-bestmodels.cfitted.CI %>% group_by(time) %>% summarise(lowCI=min(curve_fitted),highCI=max(curve_fitted)) %>% ungroup
#print(bestmodels.fitted.CI)    
#print(tail(bestmodels.fitted.CI))


restemp<-res[[1]]
bestmodels_t<- restemp %>% select(rate,max) %>% pivot_wider(names_from=rate,values_from=max) %>% ungroup
print(bestmodels_t)


#calculate mean induction
bestmodels.cfitted<-predict_induction(bestmodels_t,nparms=nparms-nparms_induction,nheaders=0,nparms_induction=nparms_induction,induction_function=induction_curve_vectorized)
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
	}
if ("k12" %in% nameparms)	
	{ 
	#rescale mean induction in terms of induction*k11 (cutting flow)
	bestmodels.cfitted$curve_fitted<-bestmodels.cfitted$curve_fitted*(k11+k12)
	} else 
	{
	bestmodels.cfitted$curve_fitted<-bestmodels.cfitted$curve_fitted*(k11)
	}
bestmodels.cfitted<-left_join(bestmodels.cfitted,bestmodels.cfitted.CI)
#calculate mean trajectory with k11-no-induction
bestmodels.fitted<-predict_models(bestmodels_t,ntypes=ntypes,nparms=nrates,nheaders=0,errormatrix=errormatrix,mymodel=mymodel)
bestmodels.fitted <- bestmodels.fitted %>% left_join(bestmodels.fitted.CI)
#calculate mean trajectory with no errors
bestmodels.fitted.0er<-predict_models(bestmodels_t,ntypes=ntypes,nparms=nrates,nheaders=0,errormatrix=diag(ntypes),mymodel=mymodel,timest=seq(0,72,timeresolution))


#====PLOT TRAJECTORIES====
myplot <- bestmodels.fitted %>% ggplot(aes(x=time,y=p,colour=types)) + 
xlim(0,72.2) + scale_y_continuous("p", breaks=c(0,0.01,0.1,0.25,0.5,0.75,1), trans='sqrt') + scale_x_continuous("time (hours)", breaks=c(0,6,12,24,36,48,72), limits=c(0,72.5))+
geom_ribbon(aes(ymin = lowCI, ymax = highCI), alpha = 0.1 ,linetype="dotted")+ 
geom_line()+geom_point(data=tidy.mydata0.p)+
theme_bw()  

ggsave(myplot,filename=paste0(output_file,"_plot.trajectories.pdf"))


#+geom_point(data=datatmp) + geom_line(data=datatmpnor11,linetype="dashed") +scale_y_continuous(trans = 'sqrt') + theme_bw() + 
 #       geom_ribbon(aes(ymin = lowCI, ymax = highCI), alpha = 0.1 ,linetype="dotted") + scale_y_continuous("p", breaks=c(0,0.01,0.1,0.25,0.5,0.75,1), trans='sqrt') + scale_x_continuous("time (hours)", breaks=c(0,6,12,24,36,48,72), limits=c(0,72.5)) + ggtitle(paste(mytarget,myinduction))


#====PLOT INDUCTION CURVES====
#datatmpnor11<-bestmodels.cfitted %>% filter(target==mytarget,induction==myinduction,model=="modelDSBs1i1_realimpnor11")
myplot <- bestmodels.cfitted %>% ggplot(aes(x=time,y=curve_fitted)) + 
xlim(0,72.2) +
scale_y_continuous("p", trans='sqrt') + scale_x_continuous("time (hours)", breaks=c(0,6,12,24,36,48,72), limits=c(0,72.5)) + 
geom_line()+xlim(0,72.2) +scale_y_continuous("precise cuts/intact DNA per hour",trans = 'sqrt') + theme_bw() +
geom_ribbon(aes(ymin = lowCI, ymax = highCI), alpha=0.1,linetype="dotted")  
ggsave(myplot,filename=paste0(output_file,"_plot.induction.pdf"))


#====CALCULATE FLOWS====
print("calculate flow")
print(bestmodels_t)
#print(bestmodels.fitted)
#y.fitted<-bestmodels.fitted %>% pivot_wider(names_from="types",values_from="p")
y.fitted<-bestmodels.fitted.0er %>% pivot_wider(names_from="types",values_from="p",values_fill=0)
#print(y.fitted)
if (ntypes==4){
names(y.fitted)<-c("time","y1","y2","y3","y4")
y.fitted <- y.fitted %>% group_by(time) %>% summarise(y1=sum(y1),y2=sum(y2),y3=sum(y3),y4=sum(y4)) #%>% select(y1,y2,y3,y4)
#print(y.fitted)
} else if (ntypes==3){
names(y.fitted)<-c("time","y1","y2","y3")
y.fitted <- y.fitted %>% group_by(time) %>% summarise(y1=sum(y1),y2=sum(y2),y3=sum(y3))
}
bestmodels_t0<-bestmodels_t
for (ipar in 1:(nparms-nparms_induction)){bestmodels_t0[ipar]<-0}

flow_l<-list()
#substitute param to obtain a full induction so that 1.
bestmodels_t0[which(names(bestmodels_t0)=="r0")]<-999999999
bestmodels_t0[which(names(bestmodels_t0)=="r2")]<-0
bestmodels_t0[which(names(bestmodels_t0)=="x0")]<-0.0000001
bestmodels.cfitted<-predict_induction(bestmodels_t,nparms=nparms-nparms_induction,nheaders=0,nparms_induction=nparms_induction,induction_function=induction_curve_vectorized,times=seq(0,72,timeresolution))
#rescale k11 to have induction inside
for (ipar in 1:(nparms-nparms_induction)){
bestmodels_tt<-bestmodels_t0
bestmodels_tt[ipar]<-bestmodels_t[ipar]
bestmodels_tt.df<-sapply(1:length(seq(0,72,timeresolution)), function(x) unlist(bestmodels_tt)) %>% t %>% as.data.frame
if (names(bestmodels_t0)[ipar]=="k11" || names(bestmodels_t0)[ipar]=="k12" )
{
bestmodels_tt.df[,ipar]<-bestmodels_tt.df[,ipar]*as.numeric(bestmodels.cfitted$curve_fitted)
}
changess<-lapply(1:nrow(y.fitted), function(x) predict_models(bestmodels_tt.df[x,],ntypes=ntypes,nparms=nrates,nheaders=0,errormatrix=diag(ntypes),mymodel=mymodel,yinit=unlist(y.fitted[x,-1]),timest=c(0,timeresolution)))
#flow_l[[nameparms[ipar]]]<-sum(sapply(changes,function(x) sum(abs(x$p[x$time==0.1]-x$p[x$time==0]))/2))
flow_l[[nameparms[ipar]]]<-sum(sapply(changess,function(x) max(abs(x$p[x$time==timeresolution]-x$p[x$time==0]))))
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
AIC<-2*(ncol(res[[2]])-4)-2*res[[2]]$maxll[1]
AICtable<-data.frame(AICt="AIC",AICv=AIC)
write.table(AICtable,file=paste0(output_file,"_plot.flow.tab"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(unlist(flow_l),file=paste0(output_file,"_plot.flow.tab"),quote=FALSE,sep="\t",append=TRUE,col.names=FALSE)


