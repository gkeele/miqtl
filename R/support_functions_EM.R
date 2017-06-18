make.EM.alt.formula <- function(formula, X){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- "y ~"
  this.formula <- as.formula(paste(this.formula.string, paste(gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE), collapse=" + ")))
  return(this.formula)
}
rotate.full.to.add.data <- function(diplotype.data){
  rotate.diplotypes.to.dosages <- straineff.mapping.matrix()
  dosage.data <- as.matrix(diplotype.data[,-1]) %*% rotate.diplotypes.to.dosages
  dosage.data <- data.frame(SUBJECT.NAME=diplotype.data$SUBJECT.NAME, dosage.data)
  names(dosage.data)[-1] <- LETTERS[1:8]
  return(dosage.data)
}
get.weights.from.diplotype.prob.data <- function(diplotype.data){
  probs.diplotype <- as.matrix(diplotype.data[,-1])
  weights <- as.vector(t(probs.diplotype))
  return(weights)
}
get.expanded.X <- function(diplotype.data, coefficients=NULL){
  if(!is.null(coefficients)){
    diplotype.data <- diplotype.data[,c("SUBJECT.NAME", names(coefficients))]
  }
  
  n <- nrow(diplotype.data)
  p <- ncol(diplotype.data) - 1
  X <- matrix(rep(diag(p), n), ncol=p, byrow=TRUE)
  colnames(X) <- names(diplotype.data)[-1]
  
  expand.subjects <- paste(rep(diplotype.data[,1], each=p), rep(1:p, p), sep=".")
  expand.data <- data.frame(SUBJECT.NAME=expand.subjects, X)
  return(expand.data)
}
get.expanded.y <- function(y, num.groups){
  y <- rep(y, each=num.groups)
  return(y)
}
apply.contrast <- function(X){
  C <- create.treatment.contrast.matrix(num.par=ncol(X), remove.column.index=ncol(X))
  new.X <- cbind(rep(1, nrow(X)), X %*% C)
  colnames(new.X) <- c("(Intercept)", colnames(X)[-ncol(X)])
  return(new.X)
}
create.treatment.contrast.matrix <- function(num.par, remove.column.index){
  C <- contr.treatment(n=num.par, base=remove.column.index)
  return(C)
}
E.step <- function(start.weights, y, X, coefficients, sigma2, p=36){
  mu <- X %*% coefficients
  pheno.prob.matrix <- matrix(dnorm(y, mean=mu, sd=sqrt(sigma2)), ncol=p, byrow=TRUE)
  dip.prob.matrix <- matrix(start.weights, ncol=p, byrow=TRUE)
  product.matrix <- dip.prob.matrix * pheno.prob.matrix
  scale.matrix <- matrix(rep(rowSums(product.matrix), each=p), ncol=p, byrow=TRUE)
  weight.matrix <- product.matrix/scale.matrix
  weights <- as.vector(t(weight.matrix))
  return(weights)
}
M.step <- function(weights, y, X, p=36){
  fit <- lm.wfit(x=X, y=y, w=weights)
  # fit$residuals appear to be not be weighted
  #fit$sigma2 <- sum((fit$residuals^2))/sum(weights)
  fit$sigma2 <- sum(fit$residuals * weights * fit$residuals)/sum(weights)
  return(fit)
}
remove.weights <- function(weights, previous.coefficients, coefficients){
  weight.matrix <- matrix(weights, ncol=length(previous.coefficients), byrow=TRUE)
  remove.col <- is.na(coefficients[names(previous.coefficients)])
  if(any(remove.col)){
    weight.matrix <- weight.matrix[,-which(remove.col)]
  }
  weights <- as.vector(t(weight.matrix))
  return(weights)
}
reshape.X <- function(expanded.X, previous.coefficients, coefficients, n){
  remove.col.logical <- is.na(coefficients[names(previous.coefficients)])
  remove.col <- which(remove.col.logical)
  if(any(remove.col.logical)){
    remove.row <- NULL
    for(i in 1:length(remove.col)){
      remove.row <- c(remove.row, (0:(n-1))*length(previous.coefficients) + remove.col[i])
    }
    remove.row <- sort(remove.row)
    expanded.X <- expanded.X[-remove.row, -remove.col]
  }
  return(expanded.X)
}
get.coef.dif <- function(coefficients, previous.coefficients){
  coef.dif <- coefficients[names(previous.coefficients)] - previous.coefficients
  coef.dif <- coef.dif[!(is.na(coef.dif) | is.nan(coef.dif))]
  return(coef.dif)
}
get.mixture.likelihood <- function(start.weights, expanded.X, expanded.y, coefficients, sigma2, p=36){
  #coefficients[is.na(coefficients)] <- 0
  mu <- expanded.X %*% coefficients
  #p <- length(coefficients)
  pheno.prob.matrix <- matrix(dnorm(expanded.y, mean=mu, sd=sqrt(sigma2)), ncol=p, byrow=TRUE)
  dip.prob.matrix <- matrix(start.weights, ncol=p, byrow=TRUE)
  product.matrix <- dip.prob.matrix * pheno.prob.matrix
  logLik <- sum(log(rowSums(product.matrix)))
  return(logLik)
}
pvalue.per.locus.im <- function(im.logLik, im.rank, fit0){
  stat <- -2*(fit0$logLik - im.logLik)
  stat <- ifelse(stat < 0, 0, stat)
  pvalue <- pchisq(q=stat, df=im.rank-fit0$rank, lower.tail=FALSE)
  return(pvalue)
}