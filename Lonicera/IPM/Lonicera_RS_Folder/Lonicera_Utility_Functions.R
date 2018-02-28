# Lonicera IPM utility functions

GrowFun <- function(Size, SizeNext, Model) { 
  ModSD <- sd(residuals(Model))
  NewSizeMean <- predict(Model, data.frame(Plant_Height12 = Size), type = 'response')
  out <- dnorm(SizeNext, mean = NewSizeMean, sd = ModSD)
  return(out)
}
SurvFun <- function(Size, Params) {
  if(length(Params) == 3){
    linear_p <- (Params[1] + Params[2] * Size + Params[3] * Size^2)
  } else {
    linear_p <- (Params[1] + Params[2] * Size)
  }
  
  return(1 / (1 + exp(-linear_p)))
}

FecFun <- function(Size, SizeNext, Models) {
  PRepMod <- Models$PRep
  SeedMod <- Models$Seeds
  EstProb <- Models$EstProb
  SizeDist <- Models$SizeDist
  
  out <- predict(PRepMod, data.frame(Plant_Height13 = Size), type = 'response') * 
    predict(SeedMod, data.frame(Plant_Height13 = Size), type = 'response') * 
    dnorm(SizeNext, mean = SizeDist$Mean, sd = SizeDist$SD) * 
    EstProb
  
  return(out)
}