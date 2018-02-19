## four parameter survival curve
predict.surv <- function(x, model) {
  pred <- predict(model, data.frame(Ht14 = x),
                      type = 'response')
  return(pred)
}

# Fecundity
SB.Go <- function(x, params, repro.model, fec.model) {
    predict(repro.model, data.frame(Ht15 = x), type = 'response') *
    predict(fec.model, data.frame(Ht = x), type='response') *
    params$germ 
}

# Seeds exiting
SB.Emerge <- function(x, params) {
  params$est.prob * 
  dnorm(x, params$recruit.size.mean, params$recruit.size.sd)
}

grow.prob<-function(size2, size1,
                    model) {
  
  sd.mod <- sd(residuals(model))
  Gdata <- data.frame(Ht14 = size1)
  size2.pred <- predict(model,
                        newdata = Gdata,
                        type = 'response')
  return(dnorm(size2,
               mean = size2.pred,
               sd = sd.mod))
}
