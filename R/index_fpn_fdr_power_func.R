

index_fpn_fdr_power_func <- function(sat_index,true_set){
  if(length(true_set)==0){
    fpn <- ifelse(length(sat_index)>0,1,0)
    return(fpn)
  }else{
    fdr <-length(which(sat_index %!in% true_set))/max(length(sat_index),1)
    power <- length(which(sat_index %in% true_set))/length(true_set)
    return(c(fdr,power))
  }
}
