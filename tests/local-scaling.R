library(splcr)
baseline <- getFromNamespace('.csh_baseline','splcr')
calibrate <- getFromNamespace('.csh_lc_baseline','splcr')
hazard <- getFromNamespace('.csh_lc_hazard','splcr')
cumulative <- getFromNamespace('.csh_lc_cumulative','splcr')
breaks <- getFromNamespace('.csh_breaks','splcr')
u <- c(.8,1.7,2.4,3.8,4.6)
d <- c(.1,.13,.14,.12,.11)
for(kernel in c('uniform','triangular','epanechnikov','quartic')) {
 a <- calibrate(baseline(u,d,5,2,kernel,'constant'))
 b <- calibrate(baseline(c(.4,1.4,2.8,4.2),c(.07,.09,.1,.14),5,2,kernel,'constant'))
 stopifnot(max(abs(diff(cumulative(a,a$edges))-d))<1e-10)
 points <- sort(unique(c(breaks(a$base,5),breaks(b$base,5),a$edges,b$edges)))
 integrate_parts <- function(fun,left,right) {
  pts <- sort(unique(c(left,points[points>left & points<right],right)))
  sum(vapply(seq_len(length(pts)-1L),function(i) stats::integrate(fun,pts[i],pts[i+1L],rel.tol=1e-9)$value,numeric(1)))
 }
 areas <- vapply(seq_along(d),function(i)integrate_parts(function(t)hazard(a,t),a$edges[i],a$edges[i+1L]),numeric(1))
 stopifnot(max(abs(areas-d))<1e-8)
 density <- function(t) (hazard(a,t)+hazard(b,t))*exp(-cumulative(a,t)-cumulative(b,t))
 mass <- integrate_parts(density,0,5)+exp(-cumulative(a,5)-cumulative(b,5))
 stopifnot(abs(mass-1)<1e-8)
 stopifnot(all(hazard(a,seq(0,5,length.out=101))>0),abs(cumulative(a,6)-cumulative(a,5)-hazard(a,5))<1e-10)
 cat(kernel,'local integral/probability checks passed\n')
}
set.seed(77)
x<-rnorm(120,2)
t1<-rexp(120,.15*exp(.3*x)); t2<-rexp(120,.12);ct<-rexp(120,.08)
dat<-data.frame(time=pmin(t1,t2,ct),status=ifelse(ct<pmin(t1,t2),0,ifelse(t1<t2,1,2)),x=x)
fold<-getFromNamespace('.assign_folds','splcr')(dat$status,5,9,TRUE)
for(kernel in c('uniform','triangular','epanechnikov','quartic')) {
 res<-spl_csh(~x,dat,kernel=kernel,fold_id=fold,return_scores=TRUE)
 stopifnot(res$algorithm=='local_scaled',res$kernel==kernel,all(is.finite(res$scores)),abs(res$Lhat-mean(res$scores))<1e-12)
 converted<-dat; converted$time<-365*converted$time
 alt<-spl_csh(~x,converted,kernel=kernel,fold_id=fold,return_scores=TRUE)
 ev<-dat$status>0
 for(f in unique(fold)) ev[fold==f & dat$time>max(dat$time[fold!=f])]<-FALSE
 stopifnot(max(abs(alt$scores-res$scores+ev*log(365)))<1e-8)
 cat(kernel,'CV/time-unit checks passed\n')
}
res<-spl_compare(list(M1=~x,M2=~1),dat,fold_id=fold)
stopifnot(nrow(res)==8,all(is.finite(res$Lhat)))
singleton<-getFromNamespace('.assign_folds','splcr')(c(0,0,1,1,2),2,1,TRUE)
stopifnot(length(singleton)==5,all(singleton>0))
expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
expect_error(spl_csh(~x,dat,causes=1))
expect_error(spl_csh(~x,dat,fold_id=rep(1,nrow(dat))))
expect_error(spl_csh(~x,dat,fold_id=rep(1.5,nrow(dat))))
sparse<-dat;sparse$status<-0L;sparse$status[1]<-1L
expect_error(spl_csh(~x,sparse,fold_id=fold))
