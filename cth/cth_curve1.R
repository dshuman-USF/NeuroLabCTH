#!/usr/bin/Rscript

#options(warn=1)

# stop that plot window focus stealing (or at least put it back where it was)
plot <- function(...) 
{
   graphics::plot(...)
   system("wmctrl -a :ACTIVE:")
}

go <- function()
{
   graphics.off()
}

# From R wavethresh package, name changed from guyrot, which made no sense
# by G P Nason
# rotate vector v by n steps.  n positive rotates right, n negatve rotates left
rotvec <- function(v,n)
{
   l <- length(v)
   n <- n%%l
   if (n == 0) 
      return(v)
   tmp <- v[(l - n + 1):l]
   v[(n + 1):l] <- v[1:(l - n)]
   v[1:n] <- tmp
   v
}

jupp <-function (v,a,b) {
    if (length(v) < 2)
       return(v);
    d <- c(v,b) - c(a,v)
    n <- length (v)
    log(d[2:(n+1)] / d[1:n])
}

unjupp <-function (v,a,b) {
    if (length(v) < 2)
       return(v);
    h <- exp(v)
    n <- length (v)
    z <- 1
    for (i in 1:n) {
        z <- 1 + h[n+1-i] * z
    }
    d <- rep (0, n)
    d[1] <- (b - a) / z
    for (i in 1:(n-1)) {
        d[i+1] <- d[i] * h[i];
    }
    a + cumsum (d)
}

# Modified, original version from AER package
loc_dispersiontest <- function(refit,rebin,alternative=c("greater","two.sided","less")) 
{
    alternative <- match.arg(alternative)
    y <- rebin
    yhat <- refit
    aux <- ((y - yhat)^2 - y)/yhat

    STAT <- sqrt(length(aux)) * mean(aux)/sd(aux)
    NVAL <- c(dispersion = 1)
    EST <- c(dispersion = mean(aux) + 1)

    rval <- list(statistic = c(z = STAT), p.value = switch(alternative, 
        greater = pnorm(STAT, lower.tail = FALSE), two.sided = pnorm(abs(STAT), 
            lower.tail = FALSE) * 2, less = pnorm(STAT)), estimate = EST, 
        null.value = NVAL, alternative = alternative, method = switch(alternative, 
            greater = "Overdispersion test", two.sided = "Dispersion test", 
            less = "Underdispersion test"), data.name = deparse(substitute(refit)))
    class(rval) <- "htest"
    return(rval)
}


# function to be minimized
fit_func <- function (ikn,cth,mint,maxt,t)
{
   knts <- sort(unique(c(mint,unjupp(ikn,mint,maxt),maxt)))
   if (knts[1] < mint) {return (NA) }
   if (knts[length(knts)] > maxt) {return (NA) }
   X <- .Call('csplinedes',t,knts)
   m <- gam (cth ~ X-1, family=poisson)
   -logLik (m)
}

# This gets called for each CTH.  It iterates for a fixed number of knots, then
# picks the one with the largest p.value.
# If not in debugging mode, multiple instances of this run on separate cores as
# separate processes.  Nothing we do here will affect globals, etc., in the
# parent process.  Returns the fitted curve as a set of Y values in text to the
# parent or an error if it occured.
cth_poisson <- function(cth)
{
   t_time <- proc.time()
   set.seed(33620)         # usf zip code
   cthnum <- cth[1]
   currwin <- 1
   assign("cthnum",cthnum,envir=.GlobalEnv)
   cat("CTH: ",cthnum, "Starting CTH","\n");

   if (Debug)
      browser()
   cth <- cth[2:length(cth)]       # first item is cth num, remove it
   i_bins <- which(cth == -1) - 1  # -1 is marker for end of i/start of e binsk
   cth <- cth[cth != -1]           # get rid of it
   tot_bins <- length(cth)
   subs <- 6    # 6 x 6 plot window
   # If all nz bins have same value, or there are just a few bins with a value
   # (e.g. 200 bins with 1, 3 with 2), we wind up with a line.  Handle this
   # case by re-binning for dispersion test by combining bins.
   browser()
   cmax <- max(cth)
   cth_hist <- vector('integer',cmax)
   for (bar in 1:cmax)
     cth_hist[bar] <- length(which(cth==bar))
   nz_bins <- length(which(cth > 0))
   thresh <- nz_bins * 0.05   # arbitrary 5% , maybe a constant, like 10?
   sparsebins <- length(which(cth_hist < thresh))
   if (sparsebins > 0)
   {
      loc_disp = TRUE
      cat("CTH: ",cthnum," Using local dispersion test\n")
      rb <- floor(sqrt(nz_bins))
      rebin <- vector(mode='integer',len=rb+1)
      idxes <- vector(mode='integer',len=rb+1)

      for (bin in 1:tot_bins)
      {
         idx <- floor((bin/tot_bins)*rb + 1)
         if (idxes[idx] == 0)
            idxes[idx] <- bin
         rebin[idx] <- rebin[idx] + cth[bin] 
      }
      rebin[rb] <- rebin[rb] + rebin[rb+1]   # grab last tick
      rebin <- rebin[1:rb]
   }
   else
      loc_disp = FALSE

   maxt <- length(cth)
   mint <- 0
   t <- .5:(length (cth) - .5)

   if (Draw==TRUE)
   {
      x11(width=14,height=11)
      currwin <- dev.cur()
      par(mfrow=c(subs,subs), mar=c(2,2,2,2))
      plot(cth,t='s',ylim=c(0,max(cth)+max(cth)*.2),col='gray38')
      hlims=par('usr')
   }

   best_pval <- -Inf 
   best_intrn <- 0
   best_knts <- 0
   best_est <- 0

   for (intrn_knt in q_steps)
   {
      stime <- proc.time()
      cat("CTH: ",cthnum,"Trying",intrn_knt+1,"knots\n")

        # find best set of random starting knots
      best_g <- Inf;
      best_val0 <- 0;
      for (g in 1:100) # 100 arbitrary, trying 500 makes almost no difference
      {
         val0 <- sort(runif(n=intrn_knt,min=1/intrn_knt,max=maxt-1))
         val0 <- jupp(val0,mint,maxt)
         l_lik <- fit_func(val0,cth,mint,maxt,t)
         if (l_lik < best_g)
         {
            best_g <- l_lik
            best_val0 <- val0
         }
         if (intrn_knt == 0)  # this never changes between iterations, so 1 is enough
            break;
      }
      val0 <- best_val0

#         if (Debug)
         cat("CTH: ",cthnum,"Initial knots [",unjupp(val0,mint,maxt),"]\n")

      optfit <- optim(val0,fit_func,cth=cth,mint=mint,maxt=maxt,t=t,control=list(trace=0,maxit=5000))
      val0 <- unjupp(val0,mint,maxt)  # for plots
      tlen <- capture.output(print(proc.time()-stime))
      cat("CTH: ",cthnum," iteration time:\nCTH: ",cthnum,tlen[1],"\nCTH: ",cthnum,tlen[2],"\n")
      if (!is.null(optfit))
      {
         if (length(optfit$convergence) > 0 && optfit$convergence != 0)
            cat("CTH: ",cthnum,"converge FAIL\n")
            cat("CTH: ",cthnum," Required ",optfit$counts[1]," iterations\n")

         ikn <- unjupp(optfit$par,mint,maxt);
         ikn <- sort(unique(ikn))
         knts <- sort(c(mint,ikn,maxt))
         X <- .Call('csplinedes',t,knts)
         m <- gam (cth ~ X-1, family=poisson)
         fit <- exp(X %*% m$coefficients)
            # the gam fitting function appears to do something like this
            # when exp returns numbers that underflow to zero.  This later causes
            # div by zero error
         fit[which(fit < .Machine$double.eps)] <- .Machine$double.eps
         scale <- max(cth)/max(fit)

          # which dispersion test to use?
         if (loc_disp)
         {
            refit <- vector(mode="double",len=length(rebin))
            for (idx in 1:(length(idxes)-1))
               refit[idx] <- sum(fit[idxes[idx]:(idxes[idx+1]-1)])
            d_res <- loc_dispersiontest(refit,rebin)
         }
         else
            d_res <- dispersiontest(m)
         est <- d_res$estimate
         pval <- d_res$p.value
         cat("CTH: ",cthnum," q: ",intrn_knt+1," CTH:",cthnum,"disp: ",est," p value: ",pval,"\n")
         subt <- sprintf("CTH: %d  k: %d  disp: %g   p: %g",cthnum,intrn_knt+1,est,pval);
#         if (PDF)
         {
            name <- sprintf("plots/v11_plot_%d_knts_%d.pdf",cthnum,intrn_knt+1)
            cat("Saving",name,"\n")
            pdf(name)
            plot(cth,t='s',main=subt,col='gray38')
            lines(fit*scale,t='l',col="red")
            dev.off()
            dev.set(which=currwin)  # dev.off switches to another win if there is one
                                    # when we save the 1st subplot, so set it back
         }
         if (Draw==TRUE)
         {
            rbow1=rbow
            pcolor1=rbow
            for (crvs in 1:ncol(X))
            {
               if (crvs == 1)
                  plot(X[,crvs],t='l',col=rbow1[1],ylim=c(min(X),max(X)))
               else
                  lines(X[,crvs],t='l',col=rbow1[1])
               rbow1=rotvec(rbow1,-1);
            }

            newlim <- c(hlims[1],hlims[2],min(hlims[3],min(fit)),max(hlims[4],max(fit)))
            if (length(which(is.infinite(newlim))) == 0)
            {
               par(usr=newlim)
               plot(fit*scale,t='l',col=pcolor1[1],ylim=newlim[3:4],yaxs='i',main=subt)
               yval=matrix(0,ncol=length(knts))
               lines(knts,yval,col="darkgoldenrod1",t='p')
               if (length(val0) > 0)
               {
                  yval=matrix(0,ncol=length(val0))
                  lines(val0,yval,col="green",t='p')
               }
               curmfg=par('mfg')
               par(mfg=(c(1,1,subs,subs)))   # id of 1st subplot
               lines(fit*scale,t='l',col=pcolor1[1])
               par(mfg=curmfg,new=F)         # new subplot
               pcolor1=rotvec(pcolor1,-1)
            }
            else
               cat("Splines contain Inf values, can't plot\n")
         }

         if (pval > best_pval)
         {
            best_pval <- pval
            best_intrn <- intrn_knt
            best_knts <- knts
            best_est <- est
         }
         if (pval >= 0.1)
         {
            break;
         }
      }
      else
      {
         print("optfit is NULL")
      }
   }

   pval <- best_pval   # pick best pval 
   knts <- best_knts
   intrn_knt <- best_intrn
   est <- best_est
   X <- .Call('csplinedes',t,knts)
   m <- gam (cth ~ X-1, family=poisson)
   winhead=sprintf("CTH: %d knts: %d  dispersion: %g   p.value: %g",cthnum,intrn_knt+1,est,pval);

  if (Debug)
     browser()
    # and the winnah is. . .
   if (length(which(is.na(knts) == 0)))
   {
      tot_pts <- 1000   # we don't need 40,000 pts for showing plots later
      wid <- tot_pts/2  # i & e widths
      pts <- seq(.5,length(cth)-.5,length=tot_pts)  #
      X1 <- .Call('csplinedes',pts,knts)
      v1 <- X1 %*% m$coefficients
      fit <- exp(v1)
      d_res <- dispersiontest(m)
      est <- d_res$estimate
      pval <- d_res$p.value

       # no matter what start and end values we use for pts above
       # the resulting Y values are 1 based. To have the curve line up
       # with the histogram in the plot, make a set of X values starting at 1
      plot_x <- seq(1,maxt,length=tot_pts)

        # i and e are usually very different durations, which results in the
        # the e width appearing to be about twice the width of i.  The
        # curves are displayed with other CTH plots in other programs that use
        # the same bin size for i and e, so they don't look "the same."
        # The next bit refits the knots so the i and e widths are the same in the
        # plots.  As a secondary goal, we also are reducing the 40,000 or so
        # bins/points to a much smaller number.  Experiments suggest it is hard
        # to tell the difference bween 40,000 points and 1000 points.  This keeps
        # the output file reasonably small.
        # make jump from end of i to start of e same step size as e step size
        #  e_step <- (last_e - (last_i + e_step)) / (wid-1), solve for e_step

      pt_i <- seq(.5,i_bins,length=wid)
      e_step <- ((tot_bins-i_bins)/(wid-1)) / (1+(1/(wid-1)))
      pt_e <- seq(pt_i[length(pt_i)]+e_step,tot_bins,len=wid)
      ptt <- c(pt_i,pt_e)
      Xt1 <- .Call('csplinedes',ptt,knts)
      vt1 <- Xt1 %*% m$coefficients
      fitt <- exp(vt1)

if (Debug)
   browser()

#      if (PDF)
      {
         subt <- sprintf("CTH: %d  k: %d  disp: %g   p: %g",cthnum,intrn_knt+1,est,pval);
         name <- sprintf("plots/v11_final_plot_%d_knts_%d.pdf",cthnum,intrn_knt+1)
         cat("Saving",name,"\n")
         pdf(name)
         par(mfrow=c(1,2), mar=c(2,2,2,2))
         plot(cth,t='s',main=subt,col='gray38')
         lines(plot_x,fit*scale,t='l',col="red")
         plot(fitt,t='s',col='blue')
         dev.off()
      }
      if (length(which(is.infinite(fit))) > 0)
      {
         infmsg <- sprintf("CTH: %d  Best fit has Inf values, not included in output file",cthnum)
         print(infmsg)
         fit <- NA
      }
   }
   else
      fit <- NA

   if (Draw==TRUE && length(which(is.na(fit))) == 0)
   {
      x11(width=7,height=5,xpos=1020,ypos=60,title=winhead)
      newlim <- c(hlims[1],hlims[2],min(hlims[3],min(fit)),max(hlims[4],max(fit)))
      plot(cth,t='s',ylim=newlim[3:4],yaxs='i',col='gray38')
      lines(plot_x,fit*scale,t='l',col="red")
      yval=matrix(0,ncol=length(knts))
      lines(knts,yval,col="darkgoldenrod1",t='p')
      if (length(val0) > 0)
      {
         yval=matrix(0,ncol=length(val0))
         lines(val0,yval,col="green",t='p')
      }

      x11(width=7,height=5,xpos=1020,ypos=600,title=winhead)
      plot(pts,fit*scale,col="red",t='l',ylim=newlim[3:4],yaxs='i')
      yval=matrix(0,ncol=length(knts))
      lines(knts,yval,col="darkgoldenrod1",t='p')
      if (length(val0) > 0)
      {
         yval=matrix(0,ncol=length(val0))
         lines(val0,yval,col="green",t='p')
      }

      x11(width=7,height=5,xpos=1020,ypos=700,title=winhead)
      rbow1=rbow
      for (crvs in 1:ncol(X))
      {
         if (crvs == 1)
            plot(X[,crvs],t='l',col=rbow1[1],ylim=c(min(X),max(X)))
         else
            lines(X[,crvs],t='l',col=rbow1[1])
         rbow1=rotvec(rbow1,-1);
      }
      x11(width=7,height=5,xpos=1025,ypos=610,title="Plot of 1000 pts, e phase adjusted")
      plot(fitt*scale,col="blue",t='l',ylim=newlim[3:4],yaxs='i')
   }

   tlen <- capture.output(print(proc.time()-t_time))
   cat("CTH: ",cthnum," total iteration time:\nCTH: ",cthnum,tlen[1],"\nCTH: ",cthnum,tlen[2],"\n")

   # return value
   # stick some info on knts on end of return array
   c(format(fitt,digits=8),c(pval,est,intrn_knt+1))
}

# End of functions.
# IT ALL STARTS HERE

# can over-ride these on command line
q <- 10  # rather arbitrary

Draw <- F    
PDF <- T    
Debug <- F

# for debugging in rstudio

# lines that should be curves
#argv <- c("cth641_to_r.txt", "cth641_from_r.txt","Draw<-T")
#argv <- c("cth389_to_r.txt", "cth389_from_r.txt","Draw<-T")

#argv <- c("cth538_to_r.txt", "cth538_from_r.txt","Draw<-T")

# cth with lots of knots on left side
#argv <- c("cth131_to_r.txt", "cth131_from_r.txt","Draw<-T")

# near flat with lots of bumps
#argv <- c("cth503_to_r.txt", "cth503_from_r.txt","Draw<-T")


#argv <- c("rtest1_to_r.txt", "rtest1_from_r.txt","Draw<-T")
#argv <- c("c1_to_r.txt", "c1_from_r.txt","Draw<-T")
#argv <- c("allv8_14_to_r.txt","allv8_14_from_r.txt","Draw<-T")
#argv <- c("allv8_15_to_r.txt","allv8_15_from_r.txt","Draw<-T")
#argv <- c("allv8_25_to_r.txt","allv8_25_from_r.txt","Draw<-T")
#argv <- c("c75_to_r.txt", "c75_from_r.txt","Draw<-T")
#argv <- c("c18_to_r.txt", "c18_from_r.txt","Draw<-T")
#argv <- c("c821_20_to_r.txt", "c821_20_from_r.txt","Draw<-T")
#argv <- c("c821_10_20_to_r.txt", "c821_10_20_from_r.txt","Draw<-T")
#argv <- c("c856_to_r.txt", "c856_from_r.txt","Draw<-T"")
#argv <- c("allexp_v10_674_to_r.txt", "allexp_v10_674_from_r.txt","Draw<-T")

#argv <- c("shortlist_to_r.txt", "shortlist_from_r.txt")
#argv <- c("shortlist_to_r.txt", "shortlist_from_r.txt","Draw<-T")
#argv <- c("short988_to_r.txt", "short988_from_r.txt")
#argv <- c("short988_to_r.txt", "short988_from_r.txt","Draw<-T")
#argv <- c("c988_10_20_to_r.txt", "c988_10_20_from_r.txt","Draw<-T")
#argv <- c("c988_10_10_to_r.txt", "c988_10_10_from_r.txt","Draw<-T")

#argv <- c("c501_10_20_to_r.txt", "c501_10_20_from_r.txt","Draw<-T")
#argv <- c("1exp20_to_r.txt", "1exp20_from_r.txt")
#argv <- c("cth_20_to_r.txt", "cth_20_from_r.txt","Draw<-T")
#argv <- c("flattie.txt", "flattie.out.txt","Draw<-T")
#argv <- c("R501_to_r.txt", "R501_from_r.txt","Draw<-T")
#argv <- c("del_me_to_r.txt", "del_me_from_r.txt","Draw<-T")
#argv <- c("cthsynth1.txt", "ctcsynth1.txt","Draw<-T")
#argv <- c("cths1.txt", "ctcs1.txt","Draw<-T")
#argv <- c("cths1.txt", "ctcs1.txt")
#argv <- c("cths2.txt", "ctcs2.txt","Draw<-T")
#argv <- c("cthbad.txt", "ctcbad.txt","Draw<-T")
#argv <- c("noisy1.txt", "noisy1_c.txt", "Draw<-T")
#argv <- c("cth788_100.txt", "ctc788_100.txt", "Draw<-T")
#argv <- c("cth390_100.txt", "ctc390_100.txt", "Draw<-T")
#argv <- c("cth726_100.txt", "ctc726_100.txt", "Draw<-T")
#argv <- c("testxxx_to_r.txt", "testxxx_from_r", "Draw<-T")
#argv <- c("20clean_curve_to_r.txt", "20clean_curve_from_r.txt", "Draw<-T")

if (!exists("argv"))
   argv <- commandArgs (trailingOnly=T)

if (length(argv) < 2)
{
   cat("Usage: cth_curve.R infile outfile\n")
   q()
}
if (argv[1]==argv[2])
{
   cat("in and out names the same, will clobber _to_r file\n")
   q()
}
i <- grep ('<-', argv)  # cmd line can over-rides, starting q
if (length (i) > 0)
{
   eval (parse (text=argv[i]))
   argv <- argv[-i]
}

print(argv)

if (Draw == T)
{
   Debug <- T   # won't draw unless serialize cth generation
#   PDF <- T
}
suppressMessages(library(SparseM))
suppressMessages(library(AER))
suppressMessages(library(mgcv))
suppressMessages(library(splines))
suppressMessages(library(parallel))
library(Rcpp)
dyn.load('csplinedes.so')

# Read cth info.  Format is:
# Identifier for data
# data
# Break these up into a list of lines, convert to ints,
# and put them all in a list of matrices.
# Assumes NO blank lines in file.
cth_data <- readLines(argv[1])
#cth_rows <- strsplit(cth_data,'\n')
#cth_list <- vector('list',length(cth_rows))
#for (r in 1:length(cth_rows))
#{
#   tmp <- strsplit(cth_rows[[r]],' ')
#   cth_list[[r]] <- matrix(unlist(lapply(tmp,as.integer)),nrow=1)
#}
#

cth_list <- vector('list')
for (r in seq(1,length(cth_data),2))
{
   tmp <- strsplit(cth_data[r+1],' ')
   cth_list[[ cth_data[r] ]] <- matrix(unlist(lapply(tmp,as.integer)),nrow=1)
}

if (Draw)
{
   rbow=c("red","palegreen4","salmon","purple","limegreen","orange","indianred4","green4","darkslateblue","darkolivegreen3","blue","chartreuse","burlywood4","aquamarine3","darkmagenta","darkorchid2","green","lightgoldenrod3","mediumvioletred","yellow4")
}

rows <- length(cth_list)
cols <- length(cth_list[[1]])-2  # 1 num is cth num, 1 num is boundary between i & e
                                 # this is # of cols of data
if (q > cols)
   q <- cols

q_steps <- seq(0,q,1)  # internal knots

# for debugging need to use 1 core and make preschedule=T
# so we don't create a child process (which won't let us do plots)
n <- names(cth_list)
if (Debug || Draw) {
   cat("Using 1 core to serialize cth processing\n")
   res <- mclapply(cth_list,FUN=cth_poisson,mc.preschedule=T,mc.cores=1)
} else {
   cat("Using ", detectCores()-2," cores\n")
   res <- mclapply(cth_list,FUN=cth_poisson,mc.preschedule=F,mc.cores=detectCores()-2)
}
print("Saving. . .")
blank <- vector("double",4)
f <- file(argv[2],"w")
for (save in 1:rows)
{
   tmp <- as.numeric(res[[save]])   # some fits fail, 
   if (is.numeric(tmp) && !is.na(tmp))
   {
      cat(file=f, sep="", names(res[save]),"\n")
      cat(file=f, sep=" ", res[[save]],"\n")
   }
   else
   {
      cat("cth ",save, " not numeric, using zero vector","\n","error text is",res[[save]])
      cat(file=f, sep="", names(res[save]),"\n")
      cat(file=f, sep=" ", blank,"\n")
   }
}
close(f)
cat("Saved ",save," curves\n")
print("DONE")

if (Draw==TRUE)
{
      # rstudio does not handle stdin very well
   if (Sys.getenv("RSTUDIO") == "1")
      browser()
   else
   {
      print("Enter q to quit")
      key <- ' '
      while (key != 'q')
      {
         key <- scan("stdin", what=character(), n=1,quiet=T)
      }
   }
}

