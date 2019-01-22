## four parameter survival curve
predict.surv <- function(x, model) {
  
  linear_p <- predict(model, data.frame(Height = x))[ ,1]
  
  return(linear_p)
}

boot_predict_surv <- function(x, params) {
  linear_p <- exp(params[1] + params[2] * x + params[3] * x^2)
  return(linear_p/(1+linear_p))
  
}

# probability of reproduction
Reproduction.prob <- function(x,params) {
  u = exp(params$prob.repro.int + params$prob.repro.slope * x)
  return(u / (1 + u))
}
# Fecundity for SB1 and SB2
SB.Go <- function(x, params, SB.Year) {
  if(SB.Year == 2){
    Reproduction.prob(x,params) * 
    exp(params$seed.int + params$seed.slope * x) *
    params$E2
  }else if(SB.Year == 3){
    Reproduction.prob(x,params) * 
    exp(params$seed.int + params$seed.slope * x) *
    params$E3
  }
}

# Seeds emerging from seedbank
SB.Emerge <- function(x,params) {
    dnorm(x, mean = params$recruit.size.mean, sd = params$recruit.size.sd)
}

# Growth using truncated density distributions to correct for 
# eviction. This is mean for use with a GAM to predict the mean
grow.prob.trunc <- function(size2,
                            size1, 
                            model,
                            lower.size.bound,
                            upper.size.bound) {
  sse <- sum(resid(model)^2)
  
  sdhat <- sqrt(sse / df.residual(model))
  
  Gdata <- data.frame(Height = size1)
  size2.pred <- predict(model, newdata = Gdata, type ='response')
  return(dtrunc(size2,
                spec = 'norm',
                a = lower.size.bound,
                b = upper.size.bound,
                mean = size2.pred,
                sd = sdhat))
}

# Correct for eviction in the IPM script rather than the function itself. This is
# meant for use w/ a GAM
grow.prob <- function(size2, size1, model) {
  sse <- sum(resid(model)^2)
  
  sdhat <- sqrt(sse / df.residual(model))
  
  Gdata <- data.frame(Height = size1)
  size2.pred <- predict(model, newdata = Gdata, type ='response')
  return(dnorm(size2,
               mean = size2.pred,
               sd = sdhat))
}
