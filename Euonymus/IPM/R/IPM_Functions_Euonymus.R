## four parameter survival curve
predict.surv <- function(x, model) {
  
  linear_p <- predict(model, data.frame(Ht14 = x))[ ,1]
  
  return(linear_p)
}

# make nice colour palettes
make.ramp <- function(tmpcol, levels = 10, transp = NA) {
  col.out <- colorRampPalette(tmpcol)(levels)
  if(!is.na(transp)) {
    col.out <- paste(col.out, transp, sep = "")
  }
  return(col.out)
}

# probability of reproduction
Reproduction.prob<-function(x,params) {
  u=exp(params$prob.repro.int+params$prob.repro.slope*x)
  return(u/(1+u))
}
# Fecundity for SB1 and SB2
SB.Go<-function(x, params, SB.Year) {
  if(SB.Year==2){
    Reproduction.prob(x,params)* 
    exp(params$seed.int+params$seed.slope*x)*
    params$E2
  }else if(SB.Year==3){
    Reproduction.prob(x,params)* 
    exp(params$seed.int+params$seed.slope*x)*
    params$E3
  }
}

# 5. Seeds exiting
SB.Emerge<-function(x,params) {
    dnorm(x,params$recruit.size.mean,params$recruit.size.sd)
}

# Makes the vec-permutation matrix
make.K <- function(G, S) {
  n.t <- matrix(seq(G * S), ncol = 1) # ex. pop vector
  n.tk <- matrix(n.t, ncol = G, byrow = FALSE)
  n.tk.new <- matrix(t(n.tk), ncol = 1, byrow = FALSE)  # this is how you want it to look in the end.
  
  comp <- cbind(n.t, n.tk.new)  # this creates the index to move from one to the other.
  K.tmp <- matrix(0, G*S, G*S)   # empty transition matrix
  K <- replace(K.tmp, comp, 1) 
  return(K)
}


# cumulative probability of moving from one size to another
# grow.prob <- function(size, size.next, params, dist = "gamma") {
#   u <- pgamma((size.next - size), params[1], params[2])
#   return(u)
# }
# Growth using truncated density distributions to correct for 
# eviction.
grow.prob.trunc<-function(size2,size1,model,lower.size.bound,upper.size.bound){
  sse <- sum(resid(model)^2)
  sdhat <- sqrt(sse/df.residual(model))
  Gdata <- data.frame(Ht14 = size1)
  size2.pred <- predict(model, newdata = Gdata, type ='response')
  return(dtrunc(size2,
                spec = 'norm',
                a = lower.size.bound,
                b = upper.size.bound,
                mean = size2.pred,
                sd = sdhat))
}

# Correct for eviction in the IPM script rather than the function itself
grow.prob<-function(size2,size1,model){
  sse <- sum(resid(model)^2)
  sdhat <- sqrt(sse/df.residual(model))
  Gdata <- data.frame(Ht14 = size1)
  size2.pred <- predict(model, newdata = Gdata, type ='response')
  return(dnorm(size2,
               mean = size2.pred,
               sd = sdhat))
}



# find the BA of an individual
find.BA <- function(sizes){
  r <- sizes/2
  BA <- (pi * (r ^ 2)) / 100000 
  return(BA)
}

# make the F matrices
make.F.matrices <- function(G, Y, offspring.sizes, fec.thresh, f.params){
  BA <- find.BA(Y)
  f.thresh <- which(abs(Y - fec.thresh) == min(abs(Y - fec.thresh)))  # inedx of Y that is closest to fec.thresh
  BA[1:f.thresh] <- 0  # those trees with DBH < threshold get 0 fecundity
  
  Fmats <- vector('list', G)
  for(i in 1:G){
    repro <- BA * f.params[i]  
    # make a matrix - prob of each adult size producing each offspring size
    Fmats[[i]] <- t(repro %o% offspring.sizes)
  }
  return(Fmats)
}

# longevity statistics
get.longevity <- function(S, G, Utilde){
  I <- diag(S*G)  # Identity matrix dimensions of Utilde
  Utilde <- as.matrix(Utilde)
  Ntilde <- ginv(I - Utilde)  # Fundamental matrix of Utilde
  # moments
  m1 <- t(rep(1, S*G)) %*% Ntilde
  m2 <- m1 %*% (2*Ntilde - diag(S*G))
  m3 <- m1 %*% (6*(Ntilde%*%Ntilde) - 6*Ntilde + diag(S*G))
  # variance
  var <- colSums(2 * Ntilde %*% Ntilde - Ntilde) - colSums(Ntilde) * colSums(Ntilde)
  
  return(list(mu = m1, var = var, m3 = m3))
}


### passage times to sets of states
get.passage.time <- function(size.range, Y, G.int, Utilde, S, G){
  
  # Define the set of states - move from actual sizes to corresponding size bins
  min.c <- which(abs(Y- size.range[1]) == min(abs(Y- size.range[1])))
  max.c <- which(abs(Y - size.range[length(size.range)]) == 
                   min(abs(Y - size.range[length(size.range)])))
  
  # vector of states of interest
  R <- c()
  for(g in 1:length(G.int)){
    R <- c(R, seq(min.c, max.c, 1) + S*(g-1))
  }
  
  Mtilde <- 1 - colSums(Utilde)  # probability of death
  Mtilde[R] <- 0  # already in the absorbing state
  
  Uprime <- Utilde
  Uprime[ ,R] <- 0  # those in the states of interest are in the absorbing state
  Uprime[R, ] <- 0
  
  Mprime <- Mtilde  
  Mprime <- rbind(Mprime, colSums(Utilde[R, ])) 
  Mprime[2, R] <- 1
  
  Uprime <- as.matrix(Uprime)
  Ntilde.prime <- ginv(diag(G*S) - Uprime)
  Bprime <- Mprime %*% Ntilde.prime
  
  ### CONDITIONAL MARKOV CHAIN
  Utilde.c <- diag(Bprime[2,]) %*% Uprime %*% ginv(diag(Bprime[2,]))
  
  eta1 <- ginv(diag(G*S) - Utilde.c)
  
  var <- colSums(2 * (eta1 %*% eta1) - eta1) - colSums(eta1) * colSums(eta1)
  
  m1 <- colSums(eta1)
  m1[R] <- 0
  
  return(list(mu = m1, var = var))
}

# get observed size distributions
get.observed.size.dist <- function(size.vec, S, b){
  complete.size <- size.vec[!is.na(size.vec)]
  tmp <- cut(complete.size, breaks = c(0,b, Inf), labels = FALSE, 
             include.lowest = TRUE)
  S.y <- rep(0, S)  # observed summary statistic - the size distribution
  S.y[as.numeric(names(table(tmp)))] <- as.numeric(table(tmp))
  return(S.y)
}

# compare size distributions 
compare.size.dists <- function(obs, model){
  log.like <- dmultinom(obs, prob = model, log = TRUE)
  return(log.like)
}


