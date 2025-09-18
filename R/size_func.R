size_func <- function(sat_index,true_set,sim.num){
  if(length(true_set)==0){
    fpn <- length(sat_index)
    return(fpn)
  }else{
    size <- length(which(sat_index %!in% true_set))/sim.num
    power <- length(which(sat_index %in% true_set))/length(true_set)
    return(c(size,power))
  }
}
