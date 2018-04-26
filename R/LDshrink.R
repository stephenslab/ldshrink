

LDshrink <- function(haplo_panel,map_data,m=85,Ne=11490.672741,cutoff=1e-3,isGeno=NA,cov_2_cor=T,na.rm=T){
  if(is.na(isGeno)){
    isGeno <- max(haplo_panel,na.rm = na.rm)>1
  }
  stopifnot(!is.na(isGeno))
  Genomult <- ifelse(isGeno,0.5,1)
  haplo_panel <- scale(haplo_panel,center=T,scale=F)
  if(cov_2_cor){
  return(cov2cor(shrinkCov(cov(haplo_panel,use=ifelse(na.rm,"complete.obs","all.obs"))*Genomult,map_data,m,Ne,cutoff)))
  }else{
    return(shrinkCov(cov(haplo_panel,use=ifelse(na.rm,"complete.obs","all.obs"))*Genomult,map_data,m,Ne,cutoff))
  }
}