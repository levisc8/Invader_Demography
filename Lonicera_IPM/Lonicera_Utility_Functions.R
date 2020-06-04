# Lonicera IPM utility functions

GrowFun <- function(SizeNext, Size, Model) { 
  ModSD <- sd(residuals(Model))
  NewSizeMean <- predict(Model, data.frame(Height = Size), type = 'response')
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

FecFun <- function(SizeNext, Size, Models) {
  PRepMod <- Models$PRep
  SeedMod <- Models$Seeds
  EstProb <- Models$EstProb
  SizeDist <- Models$SizeDist
  
  out <- predict(PRepMod, data.frame(HeightNext = Size), type = 'response') * 
    predict(SeedMod, data.frame(HeightNext = Size), type = 'response') * 
    dnorm(SizeNext, mean = SizeDist$Mean, sd = SizeDist$SD) * 
    EstProb
  
  return(out)
}