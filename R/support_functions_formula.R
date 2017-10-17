## Formula manipulation functions
make.null.formula <- function(formula, do.augment){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  this.formula <- as.formula(ifelse(do.augment, paste0(this.formula.string, " + augment.indicator"), this.formula.string))
  return(this.formula)
}
make.alt.formula <- function(formula, X, do.augment){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  this.formula.string <- ifelse(do.augment, paste0(this.formula.string, " + augment.indicator"), this.formula.string)
  this.formula <- as.formula(paste(this.formula.string, paste(gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE), collapse=" + "), sep=" + "))
  return(this.formula)
}
make.snp.null.formula <- function(formula, condition.loci, X.list, model){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  if(!is.null(condition.loci)){
    for(i in 1:length(condition.loci)){
      if(model == "additive"){
        this.formula.string <- paste(this.formula.string, paste("cond_SNP", i, sep="_"), sep=" + ")
      }
      else if(model == "full"){
        this.formula.string <- paste(this.formula.string, c(paste("cond_SNP", i, "aa", sep="_"), paste("cond_SNP", i, "Aa", sep="_")), sep=" + ")
      }
    }
  }
  return(as.formula(this.formula))
}
make.snp.alt.formula <- function(formula, model){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  if(model == "additive"){
    this.formula <- as.formula(paste(this.formula.string, "SNP", sep=" + "))
  }
  else if(model == "full"){
    this.formula <- as.formula(paste(this.formula.string, "SNP_aa", "SNP_Aa", sep=" + "))
  }
  return(this.formula)
}
remove.whitespace.formula <- function(formula){
  formula.string <- paste0(Reduce(paste, deparse(formula)))
  formula.string <- gsub("[[:space:]]", "", formula.string)
  return(as.formula(formula.string))
}
check.for.lmer.formula <- function(formula){
  formula.string <- paste0(Reduce(paste, deparse(formula)))
  formula.string <- gsub("[[:space:]]", "", formula.string)
  use.lmer <- grepl(pattern="\\([a-zA-Z0-9\\.]+\\|[a-zA-Z0-9\\.]+\\)", x=formula.string, perl=TRUE)
  return(use.lmer)
}
process.random.formula <- function(geno.id){
  random.formula <- as.formula(paste("~", paste(geno.id, "1", sep=" - ")))
  return(random.formula)
}
