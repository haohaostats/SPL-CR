.csh_kernel <- function(z, kernel) {
  v <- switch(kernel, uniform=0*z+.5, triangular=1-abs(z),
              epanechnikov=.75*(1-z*z), quartic=15/16*(1-z*z)^2)
  if (is.null(v)) stop('Unknown kernel')
  v[abs(z)>1+1e-14] <- 0
  pmax(v,0)
}
.csh_primitive <- function(z,kernel) {
  z <- pmin(pmax(z,-1),1)
  switch(kernel,uniform=z/2,triangular=z-z*abs(z)/2,
         epanechnikov=.75*(z-z^3/3),quartic=15/16*(z-2*z^3/3+z^5/5))
}
.csh_baseline <- function(times,jumps,tau,k=5,kernel='epanechnikov',tail='error') {
  if (!kernel %in% c('uniform','triangular','epanechnikov','quartic')) stop('Unknown kernel')
  if (!tail %in% c('error','constant')) stop('Unknown tail policy')
  if (length(k)!=1 || !is.finite(k) || k<1 || k!=as.integer(k)) stop('Invalid k')
  if (length(times)!=length(jumps) || any(!is.finite(c(times,jumps,tau))) || any(jumps<0)) stop('Invalid baseline')
  keep <- jumps>0; times<-times[keep]; jumps<-jumps[keep]
  if (!length(times) || times[1]<=0 || any(diff(times)<=0) || tau<max(times)) stop('No valid positive jumps')
  structure(list(u=times,d=jumps,tau=tau,k=k,kernel=kernel,tail=tail),class='.csh_baseline')
}
.csh_params <- function(object,t) {
  u<-object$u;k<-object$k

  nl<-findInterval(t,u,left.open=TRUE)
  nr<-findInterval(t,u)
  li<-nl-k+1L;ri<-nr+k
  L<-ifelse(li>=1,u[pmax(li,1)],0)
  R<-ifelse(ri<=length(u),u[pmin(ri,length(u))],object$tau)
  list(L=L,R=R,b=pmax(t-L,R-t))
}
.csh_hazard <- function(object,t) {
  if (any(!is.finite(t)) || any(t<0)) stop('Invalid prediction time')
  if (object$tail=='error' && any(t>object$tau)) stop('Prediction exceeds training follow-up')
  t<-pmin(t,object$tau); ans<-numeric(length(t))
  if (!length(t)) return(ans)
  for (start in seq.int(1,length(t),by=1024)) {
    ids<-start:min(start+1023,length(t)); x<-t[ids]; b<-.csh_params(object,x)$b
    if (any(b<=0)) stop('Degenerate bandwidth')
    Z<-.csh_primitive((object$tau-x)/b,object$kernel)-.csh_primitive(-x/b,object$kernel)
    if (any(Z<=0)) stop('Zero boundary mass')
    z<-outer(x,object$u,'-')/b
    ans[ids]<-as.vector(.csh_kernel(z,object$kernel)%*%object$d)/(b*Z)
  }
  ans
}
.csh_breaks <- function(object,end) {

  end<-min(end,object$tau)
  edges<-sort(unique(c(0,object$u[object$u<end],end)))
  out<-edges; anchors<-c(0,object$u,object$tau)
  if (length(edges)<2) return(edges)
  for (i in seq_len(length(edges)-1)) {
    a<-edges[i];b<-edges[i+1];p<-.csh_params(object,(a+b)/2); mid<-(p$L+p$R)/2
    if (mid>a && mid<b) out<-c(out,mid)
    hi<-min(b,mid)
    if (hi>a) {cross<-(anchors+p$R)/2;out<-c(out,cross[cross>a & cross<hi])}
    lo<-max(a,mid)
    if (b>lo) {cross<-(anchors+p$L)/2;out<-c(out,cross[cross>lo & cross<b])}
  }
  sort(unique(out))
}
.csh_gauss_cache <- new.env(parent=emptyenv())
.csh_gauss <- function(order) {
  key<-as.character(order)
  if (exists(key,.csh_gauss_cache,inherits=FALSE)) return(get(key,.csh_gauss_cache))
  if (order<2 || order!=as.integer(order)) stop('Invalid quadrature order')
  ii<-seq_len(order-1); off<-ii/sqrt(4*ii^2-1)
  J<-matrix(0,order,order);J[cbind(ii,ii+1)]<-off;J[cbind(ii+1,ii)]<-off
  e<-eigen(J,symmetric=TRUE);o<-order(e$values)
  res<-list(x=e$values[o],w=2*e$vectors[1,o]^2);assign(key,res,.csh_gauss_cache);res
}
.csh_cumulative <- function(object,t,order=16) {
  if (!length(t)) return(numeric())
  if (any(!is.finite(t)) || any(t<0)) stop('Invalid prediction time')
  if (object$tail=='error' && any(t>object$tau)) stop('Prediction exceeds training follow-up')
  x<-pmin(t,object$tau)
  if (!max(x)) return(numeric(length(x)))
  edges<-sort(unique(c(.csh_breaks(object,max(x)),x)))
  half<-diff(edges)/2;mid<-(head(edges,-1)+tail(edges,-1))/2
  q<-.csh_gauss(order);points<-outer(mid,q$x,function(m,z)m)+outer(half,q$x,'*')
  hv<-matrix(.csh_hazard(object,as.vector(points)),nrow=length(mid))
  areas<-half*as.vector(hv%*%q$w)
  ans<-c(0,cumsum(areas))[match(x,edges)]
  if (object$tail=='constant') ans<-ans+(t-x)*.csh_hazard(object,object$tau)
  ans
}

.csh_lc_baseline <- function(object,order=32) {
  u<-object$u
  edges<-c(0,(head(u,-1)+tail(u,-1))/2,object$tau)
  if(any(diff(edges)<=0)) stop('Degenerate calibration cells')
  G<-.csh_cumulative(object,edges,order=order)
  area<-diff(G)
  if(any(!is.finite(area)) || any(area<=0)) stop('Zero or invalid hazard mass in a calibration cell')
  structure(list(base=object,edges=edges,G_edges=G,d=object$d,scale=object$d/area,
                 tau=object$tau,tail=object$tail,order=order),class='csh_local')
}
.csh_lc_cell <- function(object,t) pmin(findInterval(t,object$edges,rightmost.closed=TRUE),length(object$d))
.csh_lc_hazard <- function(object,t) {
  if(any(!is.finite(t))||any(t<0))stop('Invalid time')
  if(object$tail=='error'&&any(t>object$tau))stop('Prediction beyond training range')
  x<-pmin(t,object$tau);cell<-.csh_lc_cell(object,x)
  .csh_hazard(object$base,x)*object$scale[cell]
}
.csh_lc_cumulative <- function(object,t) {
  if(any(!is.finite(t))||any(t<0))stop('Invalid time')
  if(object$tail=='error'&&any(t>object$tau))stop('Prediction beyond training range')
  x<-pmin(t,object$tau);cell<-.csh_lc_cell(object,x)
  local_G<-.csh_cumulative(object$base,x,order=object$order)-object$G_edges[cell]
  before<-c(0,cumsum(object$d))[cell]
  H<-before+object$scale[cell]*local_G
  if(object$tail=='constant')H<-H+(t-x)*.csh_lc_hazard(object,object$tau)
  H
}
