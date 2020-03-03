
suppressMessages(library(SparseM))
suppressMessages(library(AER))
suppressMessages(library(mgcv))
suppressMessages(library(splines))
suppressMessages(library(parallel))
library(Rcpp)

dyn.load('glmfit.so');

x=seq(1,1000);
y=seq(2000,3000);
offset=seq(0,0,length=1000);
browser();

m <- .Call('glmfitp',x,y,offset);

# to do, need to build object for logLik 

 ret <- -logLik(mobj);

browser();



