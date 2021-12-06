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

modelDSBs1i1_realnor21 <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -((k11+k12)*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]
        dy2 <- r12*y[3]+r22*y[4]
        dy3 <- -(r11+r12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy4 <- -r22*y[4]+(k12*induction_c(t,K,x0,r0,r2))*y[1]
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


modelDSBs1i1_3x4 <- function(t, y, parms,induction_c=induction_curve) {
  with(as.list(c(y, parms)), 
       {
        dy1 <- -(k11*induction_c(t,K,x0,r0,r2) )*y[1]+r11*y[3]
        dy2 <- r12*y[3]
        dy3 <- -(r11+r12)*y[3]+(k11*induction_c(t,K,x0,r0,r2))*y[1]
        dy4 <- +0
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
if (is(xmodel)[1]!="function"){mymodel_t<-get(as.character(df[["model"]])[1])} else {
mymodel_t<-mymodel}
mydata.fitted <- ode (times = times, y = yini, func = mymodel_t, parms = res_parms)
mydata.fitted.df<-as.data.frame(mydata.fitted)
errormatrix_t<-errormatrix
if (mymodel_t=="modelDSBs1i1_3x4")
{
errormatrix_t[3,3]<-1-res_parms[nparms]
errormatrix_t[3,4]<-res_parms[nparms]
}
mydata.fitted.df[,2:ncol(mydata.fitted.df)]<-t(sapply(1:nrow(mydata.fitted.df), function(x) as.matrix(mydata.fitted.df)[x,2:ncol(mydata.fitted.df)] %*% errormatrix_t))
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
  E.matrix_t[3,2]<-erDSB2indel
  E.matrix_t[4,2]<-erDSB2indel
  E.matrix_t[3,3]<-1-erDSB2indel
  E.matrix_t[4,4]<-1-erDSB2indel
  penalty<-penalty+induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  loglik_er_f(parms,my_data,ODEfunc,E.matrix_t)-penalty
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

loglik_er_f.pen_modelinductionx3<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix,induction_curve=induction_curve_vectorized){
  mydata.1<-my_data[1:(time_courses_begins[1]-1),]		
  mydata.2<-my_data[time_courses_begins[1]:(time_courses_begins[2]-1),]		
  mydata.3<-my_data[time_courses_begins[2]:nrow(mydata),]		
  E.matrix.1<-E.matrix[1:4,]
  E.matrix.2<-E.matrix[5:8,]
  E.matrix.3<-E.matrix[9:12,]
  #xmodel<-get("modelDSBs1i1_realimprecise")
  parms.1<-parms[c(seq(1,length(parms)-6,3),(length(parms)-3):length(parms))]
  parms.2<-parms[c(seq(2,length(parms)-5,3),(length(parms)-3):length(parms))]
  parms.3<-parms[c(seq(3,length(parms)-4,3),(length(parms)-3):length(parms))]
  penalty<-induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  loglik_er_f(parms.1,mydata.1,ODEfunc,E.matrix.1)+loglik_er_f(parms.2,mydata.2,ODEfunc,E.matrix.2)+loglik_er_f(parms.3,mydata.3,ODEfunc,E.matrix.3)-penalty
  }

#function for unconstrained optimization
loglik_er_f.pen.unconstrainedopt<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix,induction_curve=induction_curve_vectorized){
  penalty<-induction_curve_vectorized(0,parms[(length(parms)-3):length(parms)])
  if (penalty>0.00001){ penalty<-(10^7)*(penalty-0.00001)^2 } else {penalty<-0}
  penalty<-penalty+sum(apply(cbind(parms,rep(0,length(parms))),MARGIN=1,min)^2)*10^16
  loglik_er_f(parms,my_data,ODEfunc,E.matrix)-penalty
  }



model2nameparams<-function(mymodel){
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
	if ( mymodel=="modelDSBs1i1_realnor21")
		{
	        nameparms <-c("k11","k12","r11","r12","r22","K","x0","r0","r2")
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
		} else
	if ( mymodel=="modelDSBs1i1_3x4")
		{
    		nameparms <-c("k11","r11","r12","K","x0","r0","r2")
		} else
	if ( mymodel=="modelDSBs1i1_realimprecise.inductionx3")
		{
	        nameparms <-c("k11","k11","k11","k12","k12","k12","rr12","rr12","rr12","r11","r11","r11","r12","r12","r12","r21","r21","r21","r22","r22","r22","K","x0","r0","r2")
		loglik_er_f.pen<-loglik_er_f.pen_modelinductionx3
		xmodel<-get("modelDSBs1i1_realimprecise")
		};
	print(mymodel)
	print(nameparms)
	return(nameparms)
}


loglik_er_f<-function(parms,my_data=mydata,ODEfunc=model1,E.matrix=error_matrix){
  #ODEfunc: a function of class desolve
  #my_data: a dataframe with "time" course (which should always include 0) as 1st column, and col-types compatible with ODEfunc
  #parms: vector specifying the set of parameters (compatible with ODEfunc)
  #ncells<-apply(mydata[,-1],MARGIN=1,sum)
  yinit<-c(1,rep(0,ncol(my_data)-2))
  out <- ode (times = my_data$time, y = yinit, func = ODEfunc, parms = parms)
  probs<-(out[,2:ncol(out)]) %*% E.matrix
  if( sum(is.nan(probs))>0){ return(-999999999)} else {
  probs[probs<10^(-9)]<-10^(-9)
  loglik<-sapply(1:nrow(my_data), function(x) dmultinom(x=my_data[x,-1],prob=probs[x,],log=TRUE))
  return(sum(loglik))
  }
}


generate_CI<-function(bestmodels_l_rates_t,inputasrates=TRUE,likfunction,npermutations=1000,nameparams=nameparms, exploration.radius=1,returnonlyCI=FALSE,addpreviouslysampled="0",normalize_k11=TRUE)
  {
  print("calculate CI")
  res<-c();res_temp<-c()
  counterpar<-1
  xxparams<-rep(0,length(nameparams))
  names(xxparams)<-nameparams
  if (inputasrates)
    {
    for (xrate in 1:length(nameparams))
            {
            xxparams[counterpar]<- bestmodels_l_rates_t %>% .[xrate,] %>% select(rates) %>% as.numeric()
            maxl<- bestmodels_l_rates_t %>% select(value) %>% head(1) %>% as.numeric()
            counterpar<-counterpar+1
          }
    } else
    {
    for (xrate in nameparams)
        {
        xxparams[counterpar]<-bestmodels_l_rates_t[[xrate]]
        counterpar<-counterpar+1
        }
    };
        print("Generate parameters to explore")
        counter00<-1
        maxl2<-likfunction(xxparams)
        xxprobs0<-xxparams^(-3)/sum(xxparams^(-3))
        mysd<-sapply(xxparams,function(z) (z+0.001)/exploration.radius/1000)
        newpars.m1<-t(sapply(1:round(npermutations/5), function(x) rnorm(n=length(xxparams),mean=xxparams,sd=mysd)))
        mysd<-sapply(xxparams,function(z) (z+0.001)/exploration.radius/10)
        newpars.m2<-t(sapply(1:round(npermutations/5), function(x) rnorm(n=length(xxparams),mean=xxparams,sd=mysd)))
        mysd<-sapply(xxparams,function(z) (z+0.001)/exploration.radius)
        newpars.m3<-t(sapply(1:round(npermutations/5), function(x) rnorm(n=length(xxparams),mean=xxparams,sd=mysd)))
        mysd<-sapply(xxparams,function(z) (z+0.001)/exploration.radius*10)
        newpars.m4<-t(sapply(1:round(npermutations/5), function(x) rnorm(n=length(xxparams),mean=xxparams,sd=mysd)))
        mysd<-sapply(xxparams,function(z) (z+0.001)/exploration.radius*100)
        newpars.m5<-t(sapply(1:round(npermutations/5), function(x) rnorm(n=length(xxparams),mean=xxparams,sd=mysd)))
        newpars.m5[newpars.m5>1000]<-1000
        newpars<-rbind(newpars.m5,newpars.m4,newpars.m3,newpars.m2,newpars.m1)
        newpars[newpars<0]<-0
        colnames(newpars)<-nameparams
        print("Compute likelihood for new parameters")
        for (iit in 1:npermutations)
            {
            newpars_t<-newpars[iit,]
            newpars_t<-sapply(newpars_t,function(x) max(x,0))
            names(newpars_t)<-nameparams
            if (newpars_t[nameparams=="r0"]==0) { newpars_t[nameparams=="r0"]<-0.0001 }
            if (sum(newpars_t>50)>0) { newpars_t[newpars_t>50]<-50 }
            #print(c(iit,newpars_t))
            #reslik<-tryCatch(likfunction(newpars_t)); #not working on subfunctions?
            reslik<-likfunction(newpars_t)
            #loglik_er_f(newpars,my_data=mydata,ODEfunc=xmodel,E.matrix=error_matrices3_l[[mytarget]][[myerror]])
            if (!"error" %in% class(reslik)){
            res_temp<-c(newpars_t,reslik,maxl,maxl2);
            if (counter00==1) {res<-res_temp } else { res<-rbind(res,res_temp) };
            counter00<-counter00+1;}
            }
        res<-as.data.frame(res);row.names(res)<-NULL;
        names(res)<-c(nameparams,"loglik","maxll","maxll.here");
        if (addpreviouslysampled!="0")
            {
            dft<-get(addpreviouslysampled) %>% select(all_of(c(nameparams,"value")));
            names(dft)<-c(nameparams,"loglik");dft$maxll<-res$maxll[1];dft$maxll.here<-res$maxll.here[1];
            res<-rbind(res,dft)
            }
        if (normalize_k11)
                {
                print("normalize k11 in terms of cutting flow")
                maxk11_norm<-as.numeric(xxparams[["k11"]]*induction_curve_vectorized(6,xxparams[(length(nameparams)-3):(length(nameparams))])[1])
                if ("k12" %in% nameparams) { maxk12_norm<-as.numeric(xxparams[["k12"]]*induction_curve_vectorized(6,xxparams[(length(nameparams)-3):(length(nameparams))])[1])}
                for (ikk in 1:nrow(res))
                        {
                        res$k11[ikk]<-as.numeric(res$k11[ikk]*induction_curve_vectorized(6,res[ikk,(length(nameparams)-3):(length(nameparams))])[1])
                        if ("k12" %in% nameparams) { res$k12[ikk]<-as.numeric(res$k12[ikk]*induction_curve_vectorized(6,res[ikk,(length(nameparams)-3):(length(nameparams))])[1]) }
                        }
                print("normalization of k11 done")
                }
        res$ok<-0
        res$ok[abs(res$maxll-res$loglik)<1.92]<-1
        if (!returnonlyCI){return(res[res$ok==1,])} else {
        #Monte Carlo exploration returns points retained in CI, so biased towards underestimation.
        #To correct run mid-point interpolation.
        maxCI.biased<-apply(res[res$ok==1,1:length(nameparams),],MARGIN=2,FUN=max)
        minCI.biased<-apply(res[res$ok==1,1:length(nameparams),],MARGIN=2,FUN=min)
        meanCI.biased<-apply(res[res$ok==1,1:length(nameparams),],MARGIN=2,FUN=median)
        maxCI.unbiased.l<-c()
        minCI.unbiased.l<-c()
        for (ipar in 1:length(nameparams))
          {
          print(nameparms[ipar])
          maxpoint<-res[res$ok==1 & res[,ipar]==maxCI.biased[ipar],1:length(nameparams)] %>% head(1) %>% unlist
          df<-as.data.frame(res[res$ok==0 & res[,ipar]>maxCI.biased[ipar],1:length(nameparams)])
          if (nrow(df)>0)
            {
            maxCI.unbiased<-apply(df,MARGIN=1,FUN=function(x) sum((x-maxpoint)^2/meanCI.biased))
	    tmin<-which.min(maxCI.unbiased)[1]
	    if (length(tmin)>0){
	    maxCI.unbiased.l[ipar]<-(df[tmin,ipar]+maxCI.biased[ipar])/2} else
	    { maxCI.unbiased.l[ipar] <-NA }
            } else { maxCI.unbiased.l[ipar]<-maxCI.biased[ipar] }
          minpoint<-res[res$ok==1 & res[,ipar]==minCI.biased[ipar],1:length(nameparams)] %>% head(1) %>% unlist
          if (minCI.biased[ipar]==0){minCI.unbiased.l[ipar]<-0} else {
            if (nrow(df)>0)
                {
                df<-as.data.frame(res[res$ok==0 & res[,ipar]<minCI.biased[ipar],1:length(nameparams)])
                minCI.unbiased<-apply(df,MARGIN=1,FUN=function(x) sum((x-minpoint)^2/meanCI.biased))
                minCI.unbiased.l[ipar]<-(df[which.min(minCI.unbiased)[1],ipar]+minCI.biased[ipar])/2
                }
                else { minCI.unbiased.l[ipar]<-minCI.biased[ipar] }
          }
        }
        xxparams[which(nameparams=="k11")]<-maxk11_norm
        if ("k12" %in% nameparams) { xxparams[which(nameparams=="k12")]<-maxk12_norm }
        return(list(data.frame(max=xxparams,CIlow=minCI.unbiased.l,CIhigh=maxCI.unbiased.l,rate=nameparams),res))
        }
}


