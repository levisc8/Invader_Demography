## four parameter survival curve
predict.surv <- function(x, model) {
  pred <- predict(model, data.frame(Height = x),
                  type = 'response')
  return(pred)
}

# Fecundity
SB.Go <- function(x, params, repro.model, fec.model) {
  predict(repro.model, data.frame(Height = x), type = 'response') *
    predict(fec.model, data.frame(Height = x), type = 'response') *
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
  Gdata <- data.frame(Height = size1)
  size2.pred <- predict(model,
                        newdata = Gdata,
                        type = 'response')
  return(dnorm(size2,
               mean = size2.pred,
               sd = sd.mod))
}

## -----------test rounding error


# ## four parameter survival curve
# predict.surv <- function(x, int, slope, slope2) {
#   pred <- 1/(1 + exp(-(int + slope * x + slope2 * x^2)))
#   return(pred)
# }
# 
# # Fecundity
# SB.Go <- function(x, f_r_int, f_r_slope, f_s_int, f_s_slope, germ) {
#   
#   1/(1 + exp(-(f_r_int + f_r_slope * x))) *
#   exp(f_s_int + f_s_slope * x) *
#   germ 
# }
# 
# # Seeds exiting
# SB.Emerge <- function(x, est.prob, f_d_mu, f_d_sd) {
#   est.prob * 
#     dnorm(x, f_d_mu, f_d_sd)
# }
# 
# grow.prob<-function(size2, size1,
#                     int, slope, gsd) {
#   
#   size2.pred <- int + slope * size1
#   
#   return(dnorm(size2,
#                mean = size2.pred,
#                sd = gsd))
# }
