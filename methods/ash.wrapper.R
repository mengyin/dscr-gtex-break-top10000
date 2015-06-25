# input: betahat - estimate of effect size
#        sebetahat - estimated sd of betahat
ash.wrapper <- function(input,args){
  library(ashr)
  library(qvalue)

  ash.params <- choose.ash.parameters(input, args)

  fit <- ash(ash.params$betahat, ash.params$sebetahat,df = ash.params$df,
            method = "fdr",
            mixcompdist = args$mixcompdist)
  #qvalue = fit$qvalue
  # return estimates of beta, and q-values for testing beta=0
  #return(list(qvalue = qvalue, qvalue.fsr = qval.from.lfdr(fit$lfsr)))
  return(list(fit=fit))
}
