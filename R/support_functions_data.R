make.processed.data <- function(formula, data, cache.subjects, K, pheno.id, geno.id){
  all.variables <- all.vars(formula)
  covariates <- all.variables[-1]
  lh.formula.string <- unlist(strsplit(Reduce(paste, deparse(formula)), split="~"))[1]
  lh.formula.string <- gsub("[[:space:]]", "", lh.formula.string)
  covariates <- unique(c(covariates, pheno.id, geno.id))
  formula.string <- paste(lh.formula.string,
                          paste(covariates, collapse="+"),
                          sep="~")
  data <- model.frame(formula(formula.string), data=data)
  names(data) <- c("y", covariates)
  # Selecting those in both data and cache
  include.subjects <- intersect(unique(as.character(data[,geno.id])), cache.subjects)
  data <- data[as.character(data[,geno.id]) %in% include.subjects,]
  #matching <- match(x=as.character(data[,geno.id]), table=include.subjects)
  #data <- data[matching,]
  if(!is.null(K)){
    ### TODO: further selection based on whether K does not have an individual in cach or data
    K <- K[include.subjects, include.subjects]
  }
  if(length(covariates) > 0){
    covariate.matrix <- matrix(NA, nrow=nrow(data), ncol=length(covariates))
    for(i in 1:length(covariates)){
      if(is.factor(data[,covariates[i]])){
        factor.counts <- table(data[,covariates[i]])
        data[,covariates[i]] <- gdata::reorder.factor(x=data[,covariates[i]], new.order=names(sort(factor.counts[factor.counts != 0], decreasing=TRUE)))
      }
    }
  }
  return(list(data=data, K=K))
}

make.simple.augment.K <- function(K, augment.n){
  if(!is.null(K)){
    original.K.names <- colnames(K)
    K <- as.matrix(bdiag(K, diag(augment.n)))
    rownames(K) <- colnames(K) <- c(original.K.names, paste0("augment.obs", 1:augment.n))
  }
  return(K)
}

make.simple.augment.data <- function(data, formula, augment.n){
  real.y.names <- data$SUBJECT.NAME
  all.variables <- all.vars(formula)
  covariates <- all.variables[-1]
  augment.y <- rep(mean(data$y), augment.n)
  augment.y.names <- paste0("augment.obs", 1:augment.n)
  covariate.matrix <- NULL
  if(length(all.variables) > 1){
    covariate.matrix <- matrix(NA, nrow=augment.n, ncol=length(covariates))
    for(i in 1:length(covariates)){
      if(is.factor(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(levels(data[,covariates[i]])[1], augment.n)
      }
      if(is.numeric(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(mean(data[,covariates[i]]), augment.n)
      }
    }
  }
  if(is.null(covariate.matrix)){
    partial.augment.data <- data.frame(augment.y, augment.y.names)
    names(partial.augment.data) <- c("y", "SUBJECT.NAME")
  }
  if(!is.null(covariate.matrix)){
    partial.augment.data <- data.frame(augment.y, covariate.matrix, augment.y.names)
    names(partial.augment.data) <- c("y", covariates, "SUBJECT.NAME")
  }
  data <- rbind(data, partial.augment.data)
  return(data)
}

make.augment.weights <- function(data, weights, augment.n, added.data.points){
  if(added.data.points == augment.n & is.null(weights)){
    weights <- NULL
  }
  else if(added.data.points != augment.n & is.null(weights)){
    weights <- c(rep(1, nrow(data) - augment.n), rep(added.data.points/augment.n, augment.n))
    names(weights) <- as.character(data$SUBJECT.NAME)
  }
  else if(!is.null(weights)){
    weights <- c(weights, rep(added.data.points/augment.n, augment.n))
    names(weights) <- as.character(data$SUBJECT.NAME)
  }
  return(weights)
}

make.full.null.augment.K <- function(K, augment.n, original.n){
  if(!is.null(K)){
    original.K.names <- colnames(K)
    K <- as.matrix(bdiag(K, diag(augment.n)))
    K[-(1:original.n), 1:original.n] <- K[1:original.n, -(1:original.n)] <- 0.5
    K[-(1:original.n), -(1:original.n)][K[-(1:original.n), -(1:original.n)] == 0] <- 0.5
    rownames(K) <- colnames(K) <- c(original.K.names, paste0("augment.obs", 1:augment.n))
  }
  return(K)
}

make.full.null.augment.data <- function(formula, data, no.augment.K, use.par, brute,
                                        original.n, augment.n, weights){
  all.variables <- all.vars(formula)
  covariates <- all.variables[-1]
  
  null.formula.no.augment <- make.null.formula(formula=formula, is.augmented=FALSE)
  fit0.no.augment <- lmmbygls(formula=null.formula.no.augment, data=data, covariates=covariates, K=no.augment.K,
                              use.par=use.par, brute=brute, null.test=TRUE)
  set.seed(seed)
  y.null.hat <- predict.lmmbygls(fit0.no.augment=fit0.no.augment, original.n=original.n, augment.n=augment.n, 
                                 covariates=covariates, weights=weights)
  real.y.names <- data$SUBJECT.NAME
  augment.y.names <- paste0("augment.obs", 1:augment.n)
  
  covariate.matrix <- NULL
  if(!is.null(covariates)){
    covariate.matrix <- matrix(NA, nrow=augment.n, ncol=length(covariates))
    for(i in 1:length(covariates)){
      if(is.factor(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(levels(data[,covariates[i]])[1], augment.n)
      }
      if(is.numeric(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(mean(data[,covariates[i]]), augment.n)
      }
    }
  }
  partial.augment.data <- data.frame(y.null.hat, covariate.matrix, augment.y.names)
  names(partial.augment.data) <- c(outcome, covariates, "SUBJECT.NAME")
  data <- rbind(data, partial.augment.data)
  data <- cbind(data, data.frame(augment.indicator=c(rep(0, original.n), rep(1, augment.n))))
  return(data)
}
