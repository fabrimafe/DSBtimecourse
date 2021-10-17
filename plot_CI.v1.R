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

#-->Logistic decay -> discard and change to make decay constant as Daniela suggests and Brinkam does
#parms_induction<-c(K=3.5,x0=1,x2=0.7,r0=0,r2=0.02)
#induction_curve<-function(x,K,x0,x2,r0,r2) K/(1+exp(-r0*(x-x0)))/(1+exp(+r2*(x-x2)))
#induction_curve_vectorized<-function(x,params)
#  {
#  K<-params[1];x0<-params[2];x2<-params[3];r0<-params[4];r2<-params[5];
#  return(K/(1+exp(-r0*(x-x0)))/(1+exp(+r2*(x-x2))))
#  }
#x<-seq(0,2,0.01)
#plot(x,induction_curve(x,parms_induction[1],parms_induction[2],parms_induction[3],parms_induction[4],parms_induction[5]),type="l")
#-->Half-life decay as in Brinkman. r2 is 1/half life
#load("RData/notebook_2errormatrices.RData")
parms_induction<-c(K=1,x0=1,r0=10,r2=1)
#induction_curve_logistic<-function(x,K,x0,x2,r0,r2) K/(1+exp(-r0*(x-x0)))/(1+exp(+r2*(x-x2)))
induction_curve<-function(x,K,x0,r0,r2) K*2^(-x*r2)/(1+exp(-r0*(x-x0)))
induction_curve_vectorized<-function(x,params) 
  {
  return(induction_curve(x,params[1],params[2],params[3],params[4]))
  }

x<-seq(0,72,0.01)
plot(x,induction_curve(x,parms_induction[1],parms_induction[2],parms_induction[3],parms_induction[4]),type="l")

model5i1 <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -(k11*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[2]
        dy2 <- -(r11+r12)*y[2]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy3 <- r12*y[2]
        list(c(dy1, dy2, dy3))
      })
}
model5i1_nor11 <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -(k11*induction_c(t,K,x0,r0,r2) )*y[1]
        dy2 <- -r12*y[2]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy3 <- r12*y[2]
        list(c(dy1, dy2, dy3))
      })
}

modelDSBs1i1_impfromcut <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -((k11+k12)*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]+r21*y[4]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r11+r12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy4 <- -(r21+r22)*y[4]+(k12*induction_c(t,K,x0,r0,r2))*y[1]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_impfrompDSB <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -(k11*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]+r21*y[4]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r11+r12+k12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy4 <- -(r21+r22)*y[4]+k12*y[3]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_imp2DSB <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -(k11*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]+r21*y[4]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r11+r12)*y[3]+k12*y[4]
        dy4 <- -(r21+r22+k12)*y[4]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_realimprecise <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -((k11+k12)*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]+r21*y[4]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r11+r12+rr12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy4 <- -(r21+r22)*y[4]+(k12*induction_c(t,K,x0,r0,r2))*y[1]+rr12*y[3]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_fullimpreciseDSB <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -((k11+k12)*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]+r21*y[4]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r11+r12+rr12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]+rr21*y[4]
        dy4 <- -(r21+r22+rr21)*y[4]+(k12*induction_c(t,K,x0,r0,r2))*y[1]+rr12*y[3]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_fullimpreciseDSB_nor11 <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -((k11+k12)*induction_c(t,K,x0,r0,r2) )*y[1]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r12+rr12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]+rr21*y[4]
        dy4 <- -(r22+rr21)*y[4]+(k12*induction_c(t,K,x0,r0,r2))*y[1]+rr12*y[3]
        list(c(dy1, dy2, dy3, dy4))
      })
}

modelDSBs1i1_realimpnor11 <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -(k11+k12)*induction_c(t,K,x0,r0,r2)*y[1]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r12+rr12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy4 <- -r22*y[4]+(k12*induction_c(t,K,x0,r0,r2))*y[1]+rr12*y[3]
        list(c(dy1, dy2, dy3, dy4))
      })
}

###################################################################################
################### DEFINE FUNCTIONS ##############################################
###################################################################################

predict_models <- function(df,nparms=6,ntypes=6,nheaders=3,errormatrix=error_matrices3_l,mymodel=0,yinit=0,timest=0) {
dfl<-unlist(df)
res_parms<-as.numeric(dfl[(nheaders+1):(nheaders+nparms)]);
names(res_parms)<-names(dfl[(nheaders+1):(nheaders+nparms)])
if (length(timest)>1){times<-timest} else {times<-seq(0,72,0.1)}
yini<-c(y1 = 1, y2 = 0, y3 = 0,y4 = 0, y5 = 0, y6 = 0)
if (ntypes==3) { yini<-c(y1 = 1, y2 = 0, y3 = 0) }
if (ntypes==4) { yini<-c(y1 = 1, y2 = 0, y3 = 0, y4 = 0) }
if (length(yinit)>1){yini<-yinit}
if (is(xmodel)[1]!="function"){
mydata.fitted <- ode (times = times, y = yini, func = get(as.character(df[["model"]])[1]), parms = res_parms)} else {
mydata.fitted <- ode (times = times, y = yini, func = mymodel, parms = res_parms)
}
mydata.fitted.df<-as.data.frame(mydata.fitted)
mydata.fitted.df[,2:ncol(mydata.fitted.df)]<-t(sapply(1:nrow(mydata.fitted.df), function(x) as.matrix(mydata.fitted.df)[x,2:ncol(mydata.fitted.df)] %*% errormatrix))
if (ntypes==6) { names(mydata.fitted.df)<-c("time","x0","x+","x-","y0","y+","y-") }
if (ntypes==3) { names(mydata.fitted.df)<-c("time","intact","DSB","indels") } 
if (ntypes==4) { names(mydata.fitted.df)<-c("time","intact","indels","preciseDSB","impreciseDSB") } 
tidy.mydata.fitted.df<-mydata.fitted.df %>% pivot_longer(names(mydata.fitted.df)[2:(ntypes+1)],names_to = "types", values_to = "p")
return(tidy.mydata.fitted.df)
}


predict_induction <- function(df,nparms=6,nparms_induction=5,nheaders=3,induction_function=induction_curve_vectorized) {
#nparms=3;nparms_induction=5;ntypes=6;induction_function=induction_curve_vectorized
#df<-bestmodels_l[1,]
dfl<-unlist(df)
res_parms<-as.numeric(dfl[(nheaders+nparms+1):(nheaders+nparms+nparms_induction)]);
names(res_parms)<-names(dfl[(nheaders+nparms+1):(nheaders+nparms+nparms_induction)])
times<-seq(0,72,0.1)
mydata.fitted <-data.frame(time=times,curve_fitted=induction_function(times,res_parms))
return(mydata.fitted)
}

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

loglik_er_f.pen<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix,induction_curve=induction_curve_vectorized){
  penalty<-induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  loglik_er_f(parms,my_data,ODEfunc,E.matrix)-penalty
  }

loglik_er_f.nopen<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix,induction_curve=induction_curve_vectorized){
  #penalty<-induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  #if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  loglik_er_f(parms,my_data,ODEfunc,E.matrix) #-penalty
  }

loglik_er_f.pen_errorDBS2indel_4states_m1<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix,induction_curve=induction_curve_vectorized){
  erDSB2indel<-parms[length(parms)]
  penalty<-0
  if (erDSB2indel>0.5) { penalty<-100000*(erDSB2indel-0.7)^2 }
  E.matrix_t<-E.matrix
  E.matrix_t[3,1]<-erDSB2indel/2
  E.matrix_t[4,1]<-erDSB2indel/2
  E.matrix_t[3,2]<-erDSB2indel/2
  E.matrix_t[3,2]<-erDSB2indel/2
  E.matrix_t[3,3]<-1-erDSB2indel
  E.matrix_t[4,4]<-1-erDSB2indel
  penalty<-penalty+induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  loglik_er_f(parms,my_data,ODEfunc,E.matrix)-penalty
  }

#function for unconstrained optimization
loglik_er_f.pen.unconstrainedopt<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix,induction_curve=induction_curve_vectorized){
  penalty<-induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  penalty<-penalty+sum(apply(cbind(parms,rep(0,length(parms))),MARGIN=1,min)^2)*10^16
  loglik_er_f(parms,my_data,ODEfunc,E.matrix)-penalty
  }


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



args <- parse_args(argv)

input_file<-args$i
output_file<-args$o
myerrorE<-args$E
data_file<-args$d
mymodel<-args$m
maxnsims<-as.numeric(args$n)


load(input_file)
#mytarget_i<-inputf$target[1]
#myinduction<-inputf$induction[1]
#myerror<-inputf$error[1]
#mymodel<-inputf$model[1]
optimize_errorDSB2indel<-0
if ( length(grep("ester1",input_file))==1){ optimize_errorDSB2indel<-1 }
if ( length(grep("ester2",input_file))==1){ optimize_errorDSB2indel<-2 }
if ( length(grep("ester3",input_file))==1){ optimize_errorDSB2indel<-3 }

########################################################################################################
##################### LOAD DATA ########################################################################
########################################################################################################

mydata<-read.table(data_file,header=TRUE)
mydata.p<-as.data.frame(cbind(mydata$time,t(apply(mydata[,2:5],MARGIN=1,FUN=function(x) x/sum(x)))))
names(mydata.p)<-c("time","y1","y2","y3","y4")
mydata0.p<-mydata.p
names(mydata0.p)<-c("time","intact","indels","preciseDSB","impreciseDSB")
tidy.mydata0.p<-mydata0.p %>% pivot_longer(names(mydata0.p)[2:5],names_to = "types", values_to = "p")
        
mydelay<-0
ntypes<-4

#mydata<-read.table(input.file, header=TRUE)
errormatrix<-as.matrix(read.table(myerrorE, header=FALSE))
if (mymodel=="model5i1" || mymodel=="model5i1_nor11")
        {
        ntypes<-3 ; load("RData/error_matrices3_l.RData");
        error_mat_l<-error_matrices3_l
        }

if (optimize_errorDSB2indel==2) { mydelay<-0 } 


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


if (mymodel=="modelDSBs1i1_fullimpreciseDSB")
    {
    nameparms <-c("k11","k12","rr12","rr21","r11","r12","r21","r22","K","x0","r0","r2")
  	} else 
if ( mymodel=="modelDSBs1i1_impfromcut" || mymodel=="modelDSBs1i1_impfrompDSB" || mymodel=="modelDSBs1i1_imp2DSB")
    {
	nameparms <-c("k11","k12","r11","r12","r21","r22","K","x0","r0","r2")
	} else
if ( mymodel=="modelDSBs1i1_realimprecise")
	{
	nameparms <-c("k11","k12","rr12","r11","r12","r21","r22","K","x0","r0","r2")
    } else
if ( mymodel=="modelDSBs1i1_realimpnor11")
    {
	nameparms <-c("k11","k12","rr12","r12","r22","K","x0","r0","r2")
	} else
if ( mymodel=="modelDSBs1i1_fullimpreciseDSB_nor11")
	{
    nameparms <-c("k11","k12","rr12","rr21","r12","r22","K","x0","r0","r2")
	} else
if ( mymodel=="model5i1")
	{
    nameparms <-c("k11","r11","r12","K","x0","r0","r2")
	} else
if ( mymodel=="model5i1_nor11")
	{
    nameparms <-c("k11","r12","K","x0","r0","r2")
	};
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
	};
xmodel=get(mymodel);
nrates<-length(nameparms)

    if (mymodel=="modelDSBs1i1_fullimpreciseDSB"){ 
        nameparms <-c("k11","k12","rr12","rr21","r11","r12","r21","r22","K","x0","r0","r2");
        ntypes<-4;nparms<-12;nheaders<-4;nparms_induction<-4;
        } else 
    if (mymodel=="modelDSBs1i1_realimprecise"){
        nameparms <-c("k11","k12","rr12","r11","r12","r21","r22","K","x0","r0","r2")
        ntypes<-4;nparms<-11;nheaders<-4;nparms_induction<-4;
        } else 
    if (mymodel=="modelDSBs1i1_realimpnor11"){
        nameparms <-c("k11","k12","rr12","r12","r22","K","x0","r0","r2")
        ntypes<-4;nparms<-9;nheaders<-4;nparms_induction<-4;
        } else 
    if (mymodel=="modelDSBs1i1_fullimpreciseDSB_nor11"){
        nameparms <-c("k11","k12","rr12","rr21","r12","r22","K","x0","r0","r2")
        ntypes<-4;nparms<-10;nheaders<-4;nparms_induction<-4;
        }

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
	k11<-simsinCI[xsims[isim],1]/induction_for_normalized_k11
	k12<-simsinCI[xsims[isim],2]/induction_for_normalized_k11
	bestmodels_t.cfitted.sims$curve_fitted<-bestmodels_t.cfitted.sims$curve_fitted*(k11+k12)
	if (isim==1) { bestmodels.cfitted.CI<-bestmodels_t.cfitted.sims } else { bestmodels.cfitted.CI<-rbind(bestmodels.cfitted.CI,bestmodels_t.cfitted.sims) }
        simsinCI[xsims[isim],1]<-k11
        simsinCI[xsims[isim],2]<-k12
	bestmodels_t.fitted.sims<-predict_models(simsinCI[xsims[isim],],ntypes=ntypes,nparms=nparms,nheaders=0,errormatrix=errormatrix,mymodel=get(mymodel)) #"E_errorsfromunbroken"
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
print(head(bestmodels_t))


#calculate mean induction
bestmodels.cfitted<-predict_induction(bestmodels_t,nparms=nparms-nparms_induction,nheaders=0,nparms_induction=nparms_induction,induction_function=induction_curve_vectorized)
#extract induction at time 6h to calculate k11-no-induction
induction_for_normalized_k11<-bestmodels.cfitted %>% filter(time==6) %>% select(curve_fitted) %>% as.numeric
k11<-bestmodels_t$k11/induction_for_normalized_k11
k12<-bestmodels_t$k12/induction_for_normalized_k11
bestmodels_t$k11<-k11
bestmodels_t$k12<-k12
#rescale mean induction in terms of induction*k11 (cutting flow)
bestmodels.cfitted$curve_fitted<-bestmodels.cfitted$curve_fitted*(k11+k12)
bestmodels.cfitted<-left_join(bestmodels.cfitted,bestmodels.cfitted.CI)
#calculate mean trajectory with k11-no-induction
bestmodels.fitted<-predict_models(bestmodels_t,ntypes=ntypes,nparms=nparms,nheaders=0,errormatrix=errormatrix,mymodel=get(mymodel))
bestmodels.fitted <- bestmodels.fitted %>% left_join(bestmodels.fitted.CI)
#calculate mean trajectory with no errors
bestmodels.fitted.0er<-predict_models(bestmodels_t,ntypes=ntypes,nparms=nparms,nheaders=0,errormatrix=diag(4),mymodel=get(mymodel),timest=seq(0,2,0.0002))


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
names(y.fitted)<-c("time","y1","y2","y3","y4")
y.fitted <- y.fitted %>% group_by(time) %>% summarise(y1=sum(y1),y2=sum(y2),y3=sum(y3),y4=sum(y4)) #%>% select(y1,y2,y3,y4)
#print(y.fitted)
bestmodels_t0<-bestmodels_t
for (ipar in 1:(nparms-nparms_induction)){bestmodels_t0[ipar]<-0}

flow_l<-list()
#abundance_l<-list()
#change_l<-list()
for (ipar in 1:(nparms-nparms_induction)){
bestmodels_tt<-bestmodels_t0
bestmodels_tt[ipar]<-bestmodels_t[ipar]
changes<-apply(y.fitted[,-1], MARGIN=1,FUN=function(x) predict_models(bestmodels_tt,ntypes=ntypes,nparms=nparms,nheaders=0,errormatrix=diag(4),mymodel=get(mymodel),yinit=c(x),timest=c(0,0.0002)))
#flow_l[[nameparms[ipar]]]<-sum(sapply(changes,function(x) sum(abs(x$p[x$time==0.1]-x$p[x$time==0]))/2))
flow_l[[nameparms[ipar]]]<-sum(sapply(changes,function(x) sum(max(x$p[x$time==0.0002]-x$p[x$time==0]))))
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
write.table(unlist(flow_l),file=paste0(output_file,"_plot.flow.tab"))


