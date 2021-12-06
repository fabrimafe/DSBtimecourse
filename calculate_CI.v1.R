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

###################################################################################
################### DEFINE FUNCTIONS ##############################################
###################################################################################

predict_models <- function(df,nparms=6,ntypes=6,nheaders=3,errormatrix=error_matrices3_l,mymodel=0) {
dfl<-unlist(df)
res_parms<-as.numeric(dfl[(nheaders+1):(nheaders+nparms)]);
names(res_parms)<-names(dfl[(nheaders+1):(nheaders+nparms)])
times<-seq(0,72,0.1)
yini<-c(y1 = 1, y2 = 0, y3 = 0,y4 = 0, y5 = 0, y6 = 0)
if (ntypes==3) { yini<-c(y1 = 1, y2 = 0, y3 = 0) }
if (ntypes==4) { yini<-c(y1 = 1, y2 = 0, y3 = 0, y4 = 0) }
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

loglik_er_f.pen_errorimpDSB2pDSB_4states_m1<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix,induction_curve=induction_curve_vectorized){
  erDSB2indel<-parms[length(parms)]
  penalty<-0
  if (erDSB2indel>0.5) { penalty<-100000*(erDSB2indel-0.7)^2 }
  E.matrix_t<-E.matrix
  E.matrix_t[3,3]<-1-erDSB2indel
  E.matrix_t[3,4]<-erDSB2indel
  penalty<-penalty+induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  loglik_er_f(parms,my_data,ODEfunc,E.matrix_t)-penalty
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

argv <- add_argument(argv, "-i", help="input_file; optimization file, output from optimization.R")
argv <- add_argument(argv, "-d", help="data_file; time course used to calculate likelihood with optimization.R")
argv <- add_argument(argv, "-o", help="output_file")
argv <- add_argument(argv, "-n", help="number of sampled points")
argv <- add_argument(argv, "-E", help="error matrix")
argv <- add_argument(argv, "-m", help="model - mandatory, with no default")
argv <- add_argument(argv, "-l", help="switch to change likelihood function. To estimate a common error from DSB select 1. Default is no estimate from data, only from controls, to sample 0.2h after induction (0). To set induction curve without delay select 2. To model imprecise DSB as misread precise DSB select 3", default=0)


args <- parse_args(argv)

n.max<-as.numeric(args$n)
input_file<-args$i
output_file<-args$o
myerrorE<-args$E
data_file<-args$d
mymodel<-args$m
optimize_errorDSB2indel<-as.numeric(args$l)


inputf<-read.table(input_file,header=TRUE)
if (is.na(optimize_errorDSB2indel)) { optimize_errorDSB2indel<-0 }






########################################################################################################
##################### LOAD DATA ########################################################################
########################################################################################################
xmodel<-0
mydata<-read.table(data_file,header=TRUE)
time_courses_begins<-c(which(mydata$time[-1]-mydata$time[-length(mydata$time)]<0)+1)
mydelay<-0
ntypes<-4

#mydata<-read.table(input.file, header=TRUE)
errormatrix<-as.matrix(read.table(myerrorE, header=FALSE))
nameparms<-model2nameparams(mymodel)

if (mymodel=="model5i1" || mymodel=="model5i1_nor11")
                {
                ntypes<-3
                };
if ( mymodel=="modelDSBs1i1_realimprecise.inductionx3")
                {
                loglik_er_f.pen<-loglik_er_f.pen_modelinductionx3
                xmodel<-get("modelDSBs1i1_realimprecise")
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
	}  else 
if ( optimize_errorDSB2indel==3 )
	{
	loglik_er_f.pen<-loglik_er_f.pen_errorimpDSB2pDSB_4states_m1
	};

if (optimize_errorDSB2indel==2) { mydelay<-0 } 

if (!"function" %in% is(xmodel)) { xmodel<-get(mymodel) }
nrates<-length(nameparms)

#################################################################################
##################### START OPTIMIZATION ########################################
#################################################################################

maxl<-max(inputf$value)
print(maxl)
#bestmodel<-inputf[tail(order(inputf$value)),]
bestmodel<-inputf[maxl==inputf$value,][1,]
#bestmodels_l_rates_t<-bestmodels_l_rates %>%  filter(target==mytarget_i,induction==myinduction)

print("start calculations");
res<-generate_CI(bestmodel,inputasrates=FALSE,likfunction=function(zz) loglik_er_f.pen(zz,my_data=mydata,ODEfunc=xmodel,E.matrix=errormatrix),npermutations=n.max,returnonlyCI=TRUE,addpreviouslysampled="inputf")

#print(paste0("finished calculations, now saving to file ",output_file))
#print(res[[1]])
write.table(res[[1]],file=output_file,quote=FALSE,row.names=FALSE,sep="\t")
#write.table(res[[2]],file=output_file,quote=FALSE,row.names=FALSE,sep="\t")
save(res,file=paste0(output_file,".RData"))


