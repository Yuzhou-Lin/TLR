

generate_pvalues_indp_from_norm <- function(sim.num,ws,m,delta){
  pras <- c(m/sqrt(1+delta^2),m/sqrt(1+delta^2)*delta)
  combinations <- as.matrix(expand.grid(rep(list(c(0, 1)),2)))[c(2,3,1,4),]
  # Calculate the number of units to be assigned to each set based on probabilities
  sim.num_per_set <- round(ws * sim.num)
  sets <- vector("list", length = length(ws))
  current_index <- 1
  for (i in 1:length(ws)) {
    if (sim.num_per_set[i] > 0) {
      sets[[i]] <- current_index:(current_index + sim.num_per_set[i] - 1)
      current_index <- current_index + sim.num_per_set[i]
    } else {
      sets[[i]] <- integer(0)  # Empty set
    }
  }
  mu_M <- sample(c(1,-1),length(c(unlist(sets[1]),unlist(sets[4]))),replace = T)*rnorm(length(c(unlist(sets[1]),unlist(sets[4]))),pras[1],1)
  mu_Y <- sample(c(1,-1),length(c(unlist(sets[2]),unlist(sets[4]))),replace = T)*rnorm(length(c(unlist(sets[2]),unlist(sets[4]))),pras[2],1)
  Z.M <- rnorm(sim.num,0,1); Z.M[c(unlist(sets[1]),unlist(sets[4]))] <- rnorm(length(c(unlist(sets[1]),unlist(sets[4]))),mu_M,1);
  Z.Y <- rnorm(sim.num,0,1); Z.Y[c(unlist(sets[2]),unlist(sets[4]))] <- rnorm(length(c(unlist(sets[2]),unlist(sets[4]))),mu_Y,1);
  p.M <- 2* (1 - pnorm(abs(Z.M))); p.Y <- 2* (1 - pnorm(abs(Z.Y)))
  p.M[p.M==0] <- 1e-17 ; p.Y[p.Y==0] <- 1e-17
  return(list(pvalues = cbind(p.M,p.Y),zvalues=cbind(Z.M,Z.Y)))
}
