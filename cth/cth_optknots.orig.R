#!/usr/bin/Rscript

args <- commandArgs (trailingOnly=T)
q <- 10
i <- grep ('<-', args)
if (length (i) == 1) {
    eval (parse (text=args[i]))
    args <- args[-i]
}

jupp <-function (v,a,b) {
    d <- c(v,b) - c(a,v)
    n <- length (v)
    log(d[2:(n+1)] / d[1:n])
}

unjupp <-function (v,a,b) {
    h <- exp(v)
    n <- length (v)
    z <- 1
    for (i in 1:n) {
        z <- 1 + h[n+1-i] * z
    }
    d <- rep (0, n)
    d[1] = (b - a) / z
    for (i in 1:(n-1)) {
        d[i+1] <- d[i] * h[i];
    }
    a + cumsum (d)
    
}

library(mgcv)
library(splines)
library(AER)
hist <- scan (args[1])
# hist <- scan ("pnppcth.out")
# q <- 4
hist <- hist[hist != -1];
maxt <- length (hist)
t <- .5:(length (hist) - .5)
val0<-quantile(rep(t,hist),seq(1/q,(q-1)/q,l=q-1))

fit_func <- function (ikn) {
    knts <- c(0,unjupp (ikn,0,maxt), maxt)
    X <- cSplineDes (t, knts)
    m <- gam (hist ~ X-1, family=poisson)
    -logLik (m)
}
optfit<-optim(jupp(val0,0,maxt),fit_func,control=list(trace=1,maxit=5000))

knts <- c(0,unjupp (optfit$par,0,maxt), maxt)
X <- cSplineDes (t, knts)
m <- gam (hist ~ X-1, family=poisson)
dispersiontest(m)
# z = -0.3124, p-value = 0.6226
# alternative hypothesis: true dispersion is greater than 1
# sample estimates:
# dispersion 
#  0.9978829 
plot (hist,type='s')
lines (fitted (m),col='red')
format (knts,digit=17)
write (format(knts,digits=17), file="cth_optknots.out", n=1)
# 829 2003-12-11
# q <- 6 1.012011 0.04333
# q <- 7 1.013531 0.02696
# q <- 8 1.01185  0.04546
# q <- 9 1.011485 0.05066
# 821 2005-02-17
# q <- 4 1.0412   0.000241
# q <- 5 1.03037  0.00753
# q <- 6 1.029742 0.008581
# q <- 7 1.029735 0.008623
# q <- 10 1.029612 0.008784
# write (format(knts,digits=17), file="2005-02-17_821_10", n=1)
# q <- 11 0.9996191 0.4856
# write (format(knts,digits=17), file="2005-02-17_821_11", n=1)
# pdf ("2005-02-17_821_11.pdf")
# plot (hist,type='s')
# lines (fitted (m),col='red')
# dev.off()
# q <- 12 1.001415 0.4473
# write (format(knts,digits=17), file="2005-02-17_821_12", n=1)
# q <- 20 0.9987038 0.5493
# write (format(knts,digits=17), file="2005-02-17_821_20", n=1)
