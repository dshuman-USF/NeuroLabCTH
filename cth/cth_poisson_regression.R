#!/usr/bin/Rscript

options.warn <- 1

argv <- commandArgs (trailingOnly=T)
if (length(argv) < 2)
{
   cat("Usage: cth_poison_regression.R infile outfile\n")
   q()
}
if (argv[1]==argv[2])
{
   cat("in and out names the same, will clobber in file\n")
   q()
}
q <- 10
i <- grep ('<-', argv)
if (length (i) == 1) {
    eval (parse (text=argv[i]))
    argv <- argv[-i]
}


library(mgcv)
library(splines)
library(AER)

source(argv[1])           # read histm matrix
histm<-round(histm)
rows<-nrow(histm)
cols<-ncol(histm)
outcols<-cols+1
allrows<-seq(1:rows)
q<-20                  # 20 seems reasonable upper limit
q_steps <- seq(4,q,1)  # what is min num of segs?
append_to <- F
maxt <- length(histm[1,])
mint <- 0
t <- mint:maxt

fit_func <- function (ikn) {
    knts <- sort(c(mint,unique(ikn),maxt))
    if (knts[1] < mint) { return (NA) }
    if (knts[length(knts)] > maxt) { return (NA) }
    X <- cSplineDes (t, knts)
    m<-gam (hist ~ X-1, family=poisson)
#    cat ("logLik: ", logLik(m), "\n")
    -logLik (m)
}


for(cth in allrows)
{
   cat("cth: ",cth,'\n')
   hist <- histm[cth,]
   hist[maxt+1] <- hist[1]
#x11()
#plot(hist,t='s')
   done <- F
   while (done == F)
   {
      prev <- 0
      for (guess in q_steps)
      {
         cat("trying", guess, " ")
         val0 <- quantile(rep(mint:maxt,hist),seq(1/guess,(guess-1)/guess,l=guess-1))
         optfit <- optim(val0,fit_func,control=list(trace=0,maxit=5000))
cat("converge: ",optfit$convergence,"\n");
if (optfit$convergence != 0)
   cat("converge FAIL\n");
         ikn <- sort(unique(optfit$par))
         knts <- sort(c(mint,ikn,maxt))
print(knts)
         X <- cSplineDes (t, knts)
         m<-gam (hist ~ X-1, family=poisson)
         d_res <- dispersiontest(m)
         est <- d_res$estimate
cat("dispersion:", est,'\n')
         if (prev == 0)
         {
            prev <- est
            prev_m <- m
         }
             # 1->less than 1 crossing 
         if (est <= 1.0 && prev > 1.0)
         {
            cat("best est is ",est," best iter ",guess,'\n')
            currpt <- format(fitted(m),digits=17)
            done <- T
            break
         }
         else if (est < 1.0 && est < prev) # 1st guess may be < 1
         {
            cat("best est is ",prev," best iter ",guess-1,'\n')
            currpt <- format(fitted(prev_m),digits=17)
            done <- T
            break
         }
         prev <- est
         prev_m <- m
      }
      if (done == F)  # no est <= 1, pick last in hopes it is smallest disp
      {
         cat("no est < 1, best est is ",est," best iter ",guess,'\n')
         currpt <- format(fitted(m),digits=17)
         done <- T
      }
   }
   if (append_to == F)
   {
      write(currpt,file=argv[2],ncolumns=outcols)
      append_to <- T
   }
   else
      write(currpt,file=argv[2],ncolumns=outcols,append=TRUE)
}

