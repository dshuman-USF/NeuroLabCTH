#!/usr/bin/Rscript --vanilla

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
       return(v)
    d <- c(v,b) - c(a,v)
    n <- length (v)
    log(d[2:(n+1)] / d[1:n])
}

unjupp <-function (v,a,b) {
    if (length(v) < 2)
       return(v)
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
# refit is the model
# rebin is the original data
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
# The offset var is a global of just zeros.  It's rather large
# and the same for every call.
fit_func <- function (ikn,cth,mint,maxt,t)
{
   knts <- sort(unique(c(mint,unjupp(ikn,mint,maxt),maxt)))
   if (knts[1] < mint) {return (NA) }
   if (knts[length(knts)] > maxt) {return (NA) }
   X <- csplinedesx(t,knts)
   mx <- glmfitpx(X,cth,offset)
   lmx_lik = -(mx$rank - (-2 * sum(dpois(cth,mx$fitted.values,log=TRUE)) + 2 * mx$rank)/2)
   lmx_lik
}

# This gets called for each CTH.  It iterates up to a max number of knots.  It
# stops if the p.value reaches a target value or if the max number of
# iterations is reached.  In the latter case, it picks the iteration with the
# largest p value.  If not in debugging mode, multiple instances of this run on
# separate cores as separate processes.  The mclapply function will wait for
# one to finish before starting another, so we will not have thousands of
# processes for large data sets.  Nothing we do here will affect globals,
# etc., in the parent process.
# Returns the fitted curve as a set of Y values in text to the parent or an
# error if it occured.

cth_poisson <- function(cthname,cth_list)
{
   # the openGL based plotting functions can not handle anything larger or
   # smaller than the max/min single precision float values because most GPUs
   # can only handle single precision.  If we have numbers so huge or so tiny
   # (and sometimes we do), constrain them.
   # There does not seem to be a way to get these IEEE values in R.
   # Though the max case is handled here, in practice only the min case has occured.
   max_single <- 3.4028e+38
   min_single <- 1.1755e-38
   clip_max <- max_single/2
   clip_min <- min_single/2

   if (Debug)
      browser()
     # if we have already done this one, assume a restart and we're done with this one
   save_name <- sprintf("plots/%s.crv",cthname)
   if (file.exists(save_name))
   {
      cat("Using result from previous run for",save_name,"\n");
      return(cthname)
   }
   else
   {
      cat("We must have failed on ",save_name,"\n");
   }

   cth=cth_list[[cthname]]
   cthnum <- cth[1]              # this CTH's number
   cth <- cth[2:length(cth)]     # remove it
   i_bins <- which(cth == -1) - 1  # -1 is marker for end of i/start of e bins
   cth <- cth[cth != -1]           # get rid of it
   if (sum(cth) < 5)  # need at least a few points. . .
   {
      crv_f <- file(save_name,"w")
      cat("CTH is all zeros, saving dummy final fit ",save_name,"\n")
      cat(file=crv_f, sep="", cthname,"\n")
      blank <- vector("double",4)
      cat(file=crv_f, sep=" ", blank,"\n")
      close(crv_f)
      return(cthname)
   }

   t_time <- proc.time()
   set.seed(33620)               # usf zip code
   currwin <- 1
   assign("cthnum",cthnum,envir=.GlobalEnv)
   cat("CTH: ",cthnum, "Starting CTH","\n")
   tot_bins <- length(cth)
   e_bins = tot_bins-i_bins 
   subs <- 6    # 6 x 6 plot window
     # vector of zeros for fit_func
   offset <- rep(0,tot_bins)
   assign("offset",offset,envir=.GlobalEnv)

   # If all nz bins have same value, or if there are just a few bins with a
   # value (e.g. 400 bins with 1, 3 with 2), we wind up with a line because the
   # pval is larger than our threshold.  What seems to work is that if at least
   # 2 bins have 2% of the total spikes, we can proceed, otherwise we do a
   # custom dispersion test by combining bins and later summing them to create
   # the model.
   if (Debug)
      browser()
   cmax <- max(cth)
   cth_hist <- vector('integer',cmax)
   for (bar in 1:cmax)
     cth_hist[bar] <- length(which(cth==bar))
   nz_bins <- length(which(cth > 0))
   thresh <- nz_bins * 0.02   # empirically determined 2% catches all current CTH cases 
   thresh_bins <- length(which(cth_hist > thresh))
   if (thresh_bins < 2)
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

stime <- proc.time()
   for (intrn_knt in q_steps)
   {
      cat("CTH: ",cthnum,"Trying",intrn_knt+1,"knots\n")

        # find best set of random starting knots
      best_g <- Inf;
      best_val0 <- 0;
      t0 <- proc.time()
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
      t1 <- proc.time()
      cat("CTH: ",cthnum," knot guess time: ",t1-t0,"\n")

#         if (Debug)
         cat("CTH: ",cthnum,"Initial knots [",unjupp(val0,mint,maxt),"]\n")

      t0 <- proc.time()
      optfit <- optim(val0,fit_func,cth=cth,mint=mint,maxt=maxt,t=t,control=list(trace=0,maxit=5000))
      t1 <- proc.time()
      val0 <- unjupp(val0,mint,maxt)  # for plots
      tlen <- capture.output(print(t1-t0))
      cat("CTH: ",cthnum," iteration time:\nCTH: ",cthnum,tlen[1],"\nCTH: ",cthnum,tlen[2],"\n")
      if (!is.null(optfit))
      {
         if (length(optfit$convergence) > 0 && optfit$convergence != 0)
            cat("CTH: ",cthnum,"converge FAIL\n")
            cat("CTH: ",cthnum," Required ",optfit$counts[1]," iterations\n")

         ikn <- unjupp(optfit$par,mint,maxt)
         ikn <- sort(unique(ikn))
         knts <- sort(c(mint,ikn,maxt))
         X <- csplinedesx(t,knts)
         m <- gam (cth ~ X-1, family=poisson)
         fit <- exp(X %*% m$coefficients)
            # the gam fitting function appears to do something like this
            # when exp returns numbers that underflow to zero.  This eventually
            # causes div by zero error
         fit[which(fit < .Machine$double.eps)] <- .Machine$double.eps
         scale <- max(cth) / max(fit)

          # which dispersion test to use?
         if (loc_disp)
         {
            refit <- vector(mode="double",len=length(rebin))
            for (idx in 1:(length(idxes)-1))
               refit[idx] <- sum(fit[idxes[idx]:(idxes[idx+1]-1)])
            d_res <- loc_dispersiontest(refit,rebin)
         }
         else
         {
            d_res <- dispersiontest(m)
         }
         est <- d_res$estimate
         pval <- d_res$p.value
         cat("CTH: ",cthnum," q: ",intrn_knt+1," CTH:",cthnum,"disp: ",est," p value: ",pval,"\n")
         subt <- sprintf("CTH: %d  k: %d  disp: %g   p: %g",cthnum,intrn_knt+1,est,pval)

#         if (PDF)
#         {
             # the v17 version is arbitrary, it helps me keep track of which files
             # were generated using which params
            name <- sprintf("plots/v17_plot_%d_knts_%d.pdf",cthnum,intrn_knt+1)
            cat("Saving",name,"\n")
            pdf(name)
            plot(cth,t='s',main=subt,col='gray38')
            lines(fit*scale,t='l',col="red")
            yval=matrix(0.5,ncol=length(knts))
            lines(knts,yval,pch=21,bg="darkgoldenrod1",col="darkgoldenrod1",t='p')
            if (length(val0) > 0)
            {
               yval=matrix(0,ncol=length(val0))
               lines(val0,yval,pch=21,col="green",bg='green',t='p')
            }
            dev.off()
            if (Draw==TRUE)
               dev.set(which=currwin)  # dev.off switches to another win if there is one
                                       # when we save the 1st subplot, so set it back
#         }
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
               rbow1=rotvec(rbow1,-1)
            }

            newlim <- c(hlims[1],hlims[2],min(hlims[3],min(fit*scale)),max(hlims[4],max(fit*scale))) 
            if (length(which(is.infinite(newlim))) == 0)
            {
               par(usr=newlim)
               plot(fit*scale,t='l',col=pcolor1[1],ylim=newlim[3:4],yaxs='i',main=subt)
               yval=matrix(0.5,ncol=length(knts))
               lines(knts,yval,pch=21,bg="darkgoldenrod1",col="darkgoldenrod1",t='p')
               if (length(val0) > 0)
               {
                  yval=matrix(0,ncol=length(val0))
                  lines(val0,yval,pch=21,bg='green',col="green",t='p')
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
         if (Debug)
            browser()
         if (is.na(pval))  # this can happen with very sparse CTHs
         {
            cat("CTH: ", cthnum," pval is NA, calling it quits.\n")
            best_pval <- 0;
            best_intrn <- intrn_knt
            best_knts <- knts
            best_est <- est
            break;
         }
         if (pval > best_pval)
         {
            best_pval <- pval
            best_intrn <- intrn_knt
            best_knts <- knts
            best_est <- est
         }
           # 0.9 is the v17 pval 
         if (pval >= 0.9)
         {
            break;
         }
      }
      else
      {
         print("optfit is NULL")
      }
   }
   if (Debug)
      browser()

   pval <- best_pval   # pick best pval 
   knts <- best_knts
   intrn_knt <- best_intrn
   est <- best_est
   X <- csplinedesx(t,knts)
   m <- gam (cth ~ X-1, family=poisson)
   winhead=sprintf("CTH: %d knts: %d  dispersion: %g   p.value: %g",cthnum,intrn_knt+1,est,pval)

  if (Debug)
     browser()
    # and the winnah is. . .
   if (length(which(is.na(knts) == 0)))
   {
      tot_pts <- 4000   # we don't need 40,000 pts for showing plots later
      wid <- tot_pts/2  # i & e widths
      pts <- seq(.5,length(cth)-.5,length=tot_pts)  #
      X1 <- csplinedesx(pts,knts)
      v1 <- X1 %*% m$coefficients
      fit <- exp(v1)

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
        # to tell the difference bween 40,000 points and 4000 points.  This keeps
        # the output file smaller and makes plotting a lot faster.
        # make jump from end of i to start of e same step size as e step size
        #  e_step <- (last_e - (last_i + e_step)) / (wid-1), solve for e_step
      pt_i <- seq(.5,i_bins,length=wid)
      e_step <- ((e_bins)/(wid-1)) / (1+(1/(wid-1)))
      pt_e <- seq(pt_i[length(pt_i)]+e_step,tot_bins,len=wid)
      ptt <- c(pt_i,pt_e)
      Xt1 <- csplinedesx(ptt,knts)
      vt1 <- Xt1 %*% m$coefficients
      final_fit <- exp(vt1)    # this is what we are here for

      if (Debug)
         browser()

#      if (PDF)
#      {
         subt <- sprintf("CTH: %d  k: %d  disp: %g   p: %g",cthnum,intrn_knt+1,est,pval)
         name <- sprintf("plots/v17_final_plot_%d_knts_%d.pdf",cthnum,intrn_knt+1)
         cat("Saving",name,"  ",subt,"\n")
         pdf(name)
         par(mfrow=c(1,2), mar=c(2,2,2,2))
         plot(cth,t='s',main=subt,col='gray38')
         lines(plot_x,fit*scale,t='l',col="red")
         plot(final_fit,t='s',col='blue')
         if (length(val0) > 0)
         {
            yval=matrix(0,ncol=length(val0))
            lines(val0,yval,pch=21,bg='green',col="green",t='p')
         }
         dev.off()
#      }
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
      yval=matrix(0.5,ncol=length(knts))
      lines(knts,yval,pch=21,bg="darkgoldenrod1",col="darkgoldenrod1",t='p')
      if (length(val0) > 0)
      {
         yval=matrix(0,ncol=length(val0))
         lines(val0,yval,pch=21,bg='green',col="green",t='p')
      }

      x11(width=7,height=5,xpos=1020,ypos=600,title=winhead)
      plot(pts,fit*scale,col="red",t='l',ylim=newlim[3:4],yaxs='i')
      yval=matrix(0.5,ncol=length(knts))
      lines(knts,yval,pch=21,bg="darkgoldenrod1",col="darkgoldenrod1",t='p')
      if (length(val0) > 0)
      {
         yval=matrix(0,ncol=length(val0))
         lines(val0,yval,pch=21,bg='green',col="green",t='p')
      }

      x11(width=7,height=5,xpos=1020,ypos=700,title=winhead)
      rbow1=rbow
      for (crvs in 1:ncol(X))
      {
         if (crvs == 1)
            plot(X[,crvs],t='l',col=rbow1[1],ylim=c(min(X),max(X)))
         else
            lines(X[,crvs],t='l',col=rbow1[1])
         rbow1=rotvec(rbow1,-1)
      }
      if (i_bins > e_bins)
         x11(width=7,height=5,xpos=1025,ypos=610,title="I bins adjusted to be same width as E bins")
      else
         x11(width=7,height=5,xpos=1025,ypos=610,title="E bins adjusted to be same width as I bins")
      plot(final_fit*scale,col="blue",t='l',ylim=newlim[3:4],yaxs='i')
   }

   tlen <- capture.output(print(proc.time()-t_time))
   cat("CTH: ",cthnum," total time for all iterations :\nCTH: ",cthnum,tlen[1],"\nCTH: ",cthnum,tlen[2],"\n")

   # stick some info on knts on end of return array
   save_val <- c(format(final_fit,digits=8),c(pval,est,intrn_knt+1))

   crv_f <- file(save_name,"w")
   cat("Saving final fit ",save_name,"\n")
   cat(file=crv_f, sep="", cthname,"\n")
   tmp <- as.numeric(save_val)  # some fits fail, 
   if (is.numeric(tmp) && !is.na(tmp))
   {
      tmp[which(tmp >= max_single)] = clip_max
      tmp[which(tmp <= min_single)] = clip_min
      cat(file=crv_f, sep=" ", save_val,"\n")
   }
   else
   {
      blank <- vector("double",4)
      cat(file=crv_f, sep=" ", blank,"\n")
   }
   close(crv_f)
   cthname
}

# End of functions.
# IT ALL STARTS HERE

# can over-ride these on command line

# q is the internal knot count
# 14 is the v17 iteration limit  Note this means
# there will be up to 15 curves internally 0 to 14,
# but 1 to 15 in file names
q <- 14  # rather arbitrary 

Draw <- F 
PDF <- T    
Debug <- F

# for debugging in rstudio

#argv <- c("2004-537-4_to_r.txt", "2004-537-4_from_r.txt","Draw<-F","Debug<-T")
# argv <- c("k04_915_to_r.txt","k04_915_from_r.txt","Draw<-F","Debug<-T")
#argv <- c("allexp-100_ctl_to_r.txt","allexp-100_ctl_from_r.txt","Draw<-F","Debug<-T")
#argv <- c("justforfun_to_r.txt", "justforfun_from_r.txt","Draw<-T","Debug<-T")
#argv <- c("cth64_to_r.txt", "cth64_from_r.txt","Draw<-F","Debug<-F")
#    argv <- c("cth64_to_r.txt", "cth64_from_r.txt","Draw<-F","Debug<-F")
#argv <- c("cth64_to_r.txt", "cth64_from_r.txt")

#argv <- c("cth524_to_r.txt", "cth524_from_r.txt","Draw<-T","Debug<-T")

#argv <- c("cth8_to_r.txt", "cth8_from_r.txt","Draw<-T","Debug<-T")
# this a a fairly simple E phase - why so hard to fit?
#argv <- c("cth466_to_r.txt", "cth466_from_r.txt","Draw<-T","Debug<-T")

# lines that should be curves
#argv <- c("cth352_to_r.txt", "cth352_from_r.txt","Draw<-T")
#argv <- c("cth641_to_r.txt", "cth641_from_r.txt","Draw<-T")
#argv <- c("cth389_to_r.txt", "cth389_from_r.txt","Draw<-T","Debug<-T")
#argv <- c("cth806_to_r.txt", "cth806_from_r.txt","Draw<-T")
#argv <- c("cth538_to_r.txt", "cth538_from_r.txt","Draw<-T")

# cth with lots of knots on left side
#argv <- c("cth131_to_r.txt", "cth131_from_r.txt","Draw<-T")

# near flat with lots of bumps
#argv <- c("cth503_to_r.txt", "cth503_from_r.txt","Draw<-T")

#argv <- c("cth75_to_r.txt", "cth75_from_r.txt","Draw<-T","Debug<-T")

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

cat("infile:  ",argv[1], "\noutfile: ",argv[2],"\n")

if (Draw == T)
{
#   PDF <- T
}
suppressMessages(library(SparseM))
suppressMessages(library(AER))
suppressMessages(library(mgcv))
suppressMessages(library(splines))
suppressMessages(library(parallel))
# our local C++ functions
#suppressMessages(library(csplinedesx))
#suppressMessages(library(glmfitpx))
library(csplinedesx)
library(glmfitpx)

# Read cth info.  Format is 2 lines:
# Identifier for data
# data
# Break these up into a list of lines, convert data to ints,
# and put them all in a list of matrices.
# Assumes NO blank lines in file.
cth_data <- readLines(argv[1])
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
                                 # this is number of of cols of real data
if (q > cols)
   q <- cols

q_steps <- seq(0,q,1)  # internal knots

# for debugging need to use 1 core and make preschedule=T
# so we don't create a child process (which won't let us do plots)
nlist <- names(cth_list)

if (Debug || Draw) {
   cat("Using 1 core to serialize cth processing\n")
   res <- mclapply(nlist,cth_poisson,cth_list,mc.preschedule=T,mc.cores=1)
} else {
   cat("Using ", detectCores()-2," cores\n")
   res <- mclapply(nlist,cth_poisson,cth_list,mc.preschedule=F,mc.cores=detectCores()-2)
}

cat("Saving. . .\n")
out_file <- file(argv[2],"w")
if (!isOpen(out_file))
{
   cat("Error opening output file ",argv2[]," nothing saved, exiting. . .\n")
   quit(save="no",status=1)
}

if (Debug)
   browser();

count <- 0;
for (save in 1:rows)
{
   res_name <- sprintf("plots/%s.crv",res[[save]])
   if (file.exists(res_name))
   {
      res_file <- file(res_name,"r")
      fres <- readLines(res_name)
      if (length(fres) != 2)
      {
         cat("The file ", res_name, " seems to be invalid, skipping\n")
         next;
      }
      else
      {
         cat(file=out_file, sep="", fres[1],"\n")
         cat(file=out_file, sep=" ", fres[2],"\n")
         count=count+1
      }
      close(res_file)
   }
   else
      cat("File ",res_name," does not exist, skipping. . .\n")
}
close(out_file)
cat("Saved ",count," curves\n")
cat("DONE\n")
quit(save="no",status=0)

if (Draw==TRUE)
{
      # rstudio does not handle stdin very well
   if (Sys.getenv("RSTUDIO") == "1")
      browser()
   else
   {
      print("Enter q to quit")  # hang around so we can look at plots
      key <- ' '
      while (key != 'q')
      {
         key <- scan("stdin", what=character(), n=1,quiet=T)
      }
   }
}

