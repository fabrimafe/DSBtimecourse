#!/usr/bin/env Rscript
#example to run 
#module load R/4.1.0
#./timeseriesboootstraps.R -i timecourse_RNP_PhyB2.txt -o bs/ -n 100 -m 1
#note that in current implementation the ME bootstrap cannot run when for 1 time point I have only 1 sample. add a nice warning

#JACKNIFE CODE (OBSOLETE)
#calculate jacknife CI following http://www.comparingpartitions.info/?link=Tut9
#files_v<-dir()[grep("*.a",dir())]
#library(tidyverse)
#library(nnet)
#res_jk<-sapply(files_v, function(x) read.table(x,header=TRUE) %>% filter(value==max(value))) %>% t
#nparms<-which(colnames(res_jk)=="value")-1
#mat<-res_jk[,1:nparms]
#rownames(mat)<-NULL
#mat_t<-matrix(unlist(mat),nrow=nrow(mat));colnames(mat_t)<-colnames(mat)
#pseudovalues_mat<-apply(mat_t,MARGIN=2,FUN=function(x) mean(x)*length(x)-x*(length(x)-1))
#apply(pseudovalues_mat,MARGIN=2,FUN=function(x) mean(x))
#apply(pseudovalues_mat,MARGIN=2,FUN=function(x) 1.96*sqrt((sum((x-mean(x))^2/(length(x)-1)))/length(x)))
#apply(mat,MARGIN=2,FUN=function(x) mean(x))

####################################################
##########DEFINE INPUT FLAGS########################
####################################################

library(argparser, quietly=TRUE,warn.conflicts=FALSE)
argv<- arg_parser("Parse arguments")

argv <- add_argument(argv, "-i", help="input timecourse file")
argv <- add_argument(argv, "-o", help="output folder")
argv <- add_argument(argv, "-n", help="number of permutations")
argv <- add_argument(argv, "-m", help="type of bootstrap method: default is time stratified bootstrap (2); to select stationary bootstrap choose 0; 1 selects maximum entropy (ME) bootstrap; 3 selects a stratified bootstrap for which not only time but also potential batch-effects are included, and these must be specified in an additional column", default=2)
argv <- add_argument(argv, "-r", help="number of replicated runs. needed when running ME bootstrap")


args <- parse_args(argv)

npermutations<-as.numeric(args$n)
nreplicatesperpoint<-as.numeric(args$r)
inputf<-args$i
destination_path<-args$o
type_of_resampling<-as.numeric(args$m) #default is stationary bootstrap
df<-read.table(inputf,header=TRUE)

stationary_bootstrap<-function(x,p=0.1,excess_n=100,minimum.ntimes=1){
####MAIN FUNCTION FOR IMPLEMENTING STATIONARY BOOTSTRAP#####
#x is a sorted vector with times
#p is the probability of the geometric distribution controlling the average lenght of block as 1/p
#excess_n is a hyperparameter that controls the excess of lenghts to be generated. Decrease from default to speed up. Increase if bugs because too low, potentially when very short series, very low ps and huge n of resamplings.
#return new order of samples doing a stationary block bootstrap
final.ntimes<-0
while(final.ntimes<minimum.ntimes){
    nsims<-round(excess_n*length(x)*p)+3
    sims<-rgeom(nsims,p)+1
    while(sims[1]>=length(x)){
    sims<-rgeom(nsims,p)+1
    }
    nblocks<-min(which(cumsum(sims)>=length(x)))
    sims<-sims[1:nblocks]
    if ( sum(sims)==length(x)) {length_last_block<-0} else {length_last_block<-sims[nblocks]-(cumsum(sims)[nblocks]-length(x))}
    inits<-sample(1:length(x),nblocks,replace=TRUE)
    xorder<-1:length(x)
    new.samples<-lapply(1:length(inits),function(x) inits[x]:(inits+sims-1)[x])
    if (length_last_block!=0){ new.samples[[length(new.samples)]]<-new.samples[[length(new.samples)]][1:length_last_block] }
    new.samples<-sort(unlist(new.samples))
    new.samples[new.samples>max(xorder)]<-new.samples[new.samples>max(xorder)]-max(xorder)
    #new.samples is index of samples that are used
    #now, I want to add step in which samples at same times are taken at random at each iteration, to avoid unwanted spurious correlation between having tx_1 with earlier points and tx_2 with later, at least across new bootstrap events
    new.samples<-sort(sapply(new.samples,function(z) sample(which(x==x[new.samples[z]]),1)))
    final.ntimes<-length(unique(x[new.samples]))
    }
    new.samples[1:length(x)] 
}

if ( type_of_resampling == 0 ){
print("method of resampling: stationary bootstrap")

stationary.bootstraps<-sapply(1:npermutations,function(z) stationary_bootstrap(df$time,p=0.3,excess_n=10,minimum.ntimes=3))

for (ip in 1:npermutations){
df.new<-df[stationary.bootstraps[,ip],]
write.table(df.new,file=paste0(destination_path,"/timecourse.bs",ip,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
#,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
}

} else if ( type_of_resampling == 1 )
{
print("method of resampling: maximum entropy bootstrap")

library(meboot)
#install.packages("meboot")
library("meboot")
library("plm")

for (ip in 1:npermutations){
df.t<-cbind(ind1=1:nrow(df),ind2=df$time,df)
boot_l<-list()
for (ib in 1:(ncol(df)-1)){
x<-df.t[,c(1,2,3,3+ib)]
x<-pdata.frame(x)
boot_l[[ib]]<-round(meboot(x = x, reps = nreplicatesperpoint, colsubj = 2,coltimes=3, coldata = 4,trim=list(xmin=0),force.clt=FALSE,expand.sd=FALSE,reachbnd=FALSE)[,1]) ############ use this scheme. If raw count, I fear that small counts could be paired with high 
}
mat.t<-matrix(unlist(boot_l),ncol=ncol(df)-1)
mat.t[mat.t<0]<-0
df.new<-as.data.frame(cbind(df$time,mat.t))
names(df.new)<-names(df)
write.table(df.new,file=paste0(destination_path,"/timecourse.bs",ip,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
}
} else if ( type_of_resampling == 2 )
{
print("method of resampling: time-stratified bootstrap")

stratified.bootstrap<-function(df){
#CODE IMPLEMENTING TIME-STRATIFIED BOOTSTRAP
#( i.e. simple bootstrapping scheme in which data points are resampled within each time point )
times<-unique(sort(df$time))
y<-c(unlist(sapply(times, function(x) { y<-which(df$time==x); y<-sample(y,length(y),replace=TRUE); return(y)})))
return(df[y,])
}

for (ip in 1:npermutations){
df.new<-stratified.bootstrap(df)
write.table(df.new,file=paste0(destination_path,"/timecourse.bs",ip,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
}
} else if ( type_of_resampling == 3 )
{
print("method of resampling: time and batch-stratified bootstrap")

#main function for resampling at each time point
stratified.bootstrap<-function(df){
times<-unique(sort(df$time))
y<-c(unlist(sapply(times, function(x) { y<-which(df$time==x); y<-sample(y,length(y),replace=TRUE); return(y)})))
return(df[y,])
}

print(df)
time_courses_begins<-c(which(df$time[-1]-df$time[-length(df$time)]<0)+1)
print(time_courses_begins)
mydata.1<-df[1:(time_courses_begins[1]-1),]
mydata.2<-df[time_courses_begins[1]:nrow(df),]

#loop and do again for each permutation
for (ip in 1:npermutations){
df.new.1<-stratified.bootstrap(mydata.1)
df.new.2<-stratified.bootstrap(mydata.2)
write.table(rbind(df.new.1,df.new.2),file=paste0(destination_path,"/timecourse.bs",ip,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
}

}

warnings()
