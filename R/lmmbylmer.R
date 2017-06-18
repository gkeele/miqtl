lmmbylmer <- function(formula, data, REML, weights){
  if(is.null(weights)){
    lmer.fit <- lmer(formula=formula, data=data, REML=REML)
  }
  else{
    lmer.fit <- lmer(formula=formula, data=data, REML=REML, weights=weights)
  }
  return(lmer.fit)
}