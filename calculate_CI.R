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

argv <- add_argument(argv, "-i", help="input_file; optimization file, output from optimization.R")
argv <- add_argument(argv, "-d", help="data_file; time course used to calculate likelihood with optimization.R")
argv <- add_argument(argv, "-o", help="output_file")
argv <- add_argument(argv, "-n", help="number of sampled points", default=10000)
argv <- add_argument(argv, "-z", help="n parameters in induction curve", default=2)
argv <- add_argument(argv, "-E", help="error matrix")
argv <- add_argument(argv, "-m", help="model - mandatory, with no default")
argv <- add_argument(argv, "-l", help="switch to change likelihood function. To estimate a common error from DSB select 1. Default is no estimate from data, only from controls, to sample 0.2h after induction (0). To set induction curve without delay select 2. To model imprecise DSB as misread precise DSB select 3", default=0)
argv <- add_argument(argv, "-k", help="maximum value of cutting rate k11", default=10)
argv <- add_argument(argv, "-u", help="file with pre-defined induction curve")


args <- parse_args(argv)

n.max<-as.numeric(args$n)
k.max<-as.numeric(args$k)
input_file<-args$i
output_file<-args$o
myerrorE<-args$E
data_file<-args$d
mymodel<-args$m
optimize_errorDSB2indel<-as.numeric(args$l)
nparamsind<-as.numeric(args$z)
fixedinduction_file<-args$u


if (is.na(k.max)){ k.max<-10 }
inputf<-read.table(input_file,header=TRUE)
if (is.na(optimize_errorDSB2indel)) { optimize_errorDSB2indel<-0 }
if (is.na(nparamsind)){ nparamsind<-2 }
if (is.na(fixedinduction_file)){ fixedinduction<-FALSE } else
        {
        fixedinduction<-TRUE
        fixedinduction_table<-read.table(fixedinduction_file,header=TRUE)
        #preset_induction<-unname(unlist(fixedinduction_table$max))
}
define_ODE.functions(nparamsind)

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
nameparms<-model2nameparams(mymodel,nparamsind)
loglik_er_f.pen<-model2likelihoodfunction(mymodel,optimize_errorDSB2indel)
ntypes<-model2ntypes(mymodel)
xmodel<-model2xmodel(mymodel)

normalize_k11.t<-FALSE
if (nparamsind==4){ normalize_k11.t=TRUE }
if (optimize_errorDSB2indel==2) { mydelay<-0 } 

nrates<-length(nameparms)


if (fixedinduction)
		{
                print("## case 1: likelihood function using predefined induction parameters. Constrain by preventing induction parameters to change.")
                x_loglik<-function(zz)
                        {
                        for (ifi in 1:nrow(fixedinduction_table)){
                        zz[which(nameparms==fixedinduction_table$rate[ifi])]<-fixedinduction_table$max[ifi]
                        }
                        loglik_er_f.pen(zz,my_data=mydata,ODEfunc=xmodel,E.matrix=errormatrix,nind=nparamsind)
                        }
		} else
		{
                print("## case 2: likelihood functions for unconstrained optimization")
		x_loglik<-function(zz) 
			{
			loglik_er_f.pen(zz,my_data=mydata,ODEfunc=xmodel,E.matrix=errormatrix,nind=nparamsind)
			}
		}

#################################################################################
##################### START OPTIMIZATION ########################################
#################################################################################

inputf<-inputf[inputf$r0>0 & inputf$k11>0,]
inputf<-inputf[!is.na(inputf$k11),]
maxl<-max(inputf$value,na.rm=T)
print(maxl)
#bestmodel<-inputf[tail(order(inputf$value)),]
bestmodel<-inputf[maxl==inputf$value,][1,]
#bestmodels_l_rates_t<-bestmodels_l_rates %>%  filter(target==mytarget_i,induction==myinduction)
print("start calculations");
#res<-generate_CI(bestmodel,inputasrates=FALSE,likfunction=function(zz) loglik_er_f.pen(zz,my_data=mydata,ODEfunc=xmodel,E.matrix=errormatrix,nind=nparamsind),npermutations=n.max,returnonlyCI=TRUE,addpreviouslysampled="inputf",normalize_k11=normalize_k11.t)
res<-generate_CI(bestmodel,inputasrates=FALSE,likfunction=function(zz) x_loglik(zz),npermutations=n.max,returnonlyCI=TRUE,addpreviouslysampled="inputf",normalize_k11=normalize_k11.t)

#print(paste0("finished calculations, now saving to file ",output_file))
#print(res[[1]])
write.table(res[[1]],file=output_file,quote=FALSE,row.names=FALSE,sep="\t")
#write.table(res[[2]],file=output_file,quote=FALSE,row.names=FALSE,sep="\t")
save(res,file=paste0(output_file,".RData"))


