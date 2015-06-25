# joint: vash+mixash/ash
# Using R package "vash" & "mixash"/"ash"
# input: betahat - estimate of effect size
#        sebetahat - estimated sd of betahat
jointash.wrapper = function(input,args){
  library(vash)
  library(qvalue)

  ash.params <- choose.ash.parameters(input, args)

  va = vash(ash.params$sebetahat,df=ash.params$df,singlecomp=args$singlecomp)
  if (args$singlecomp==TRUE){
    library(ashr)
    sebetahat.new = c(sqrt(va$PosteriorRate/va$PosteriorShape))
    df.new = 2*va$PosteriorShape[1]
    pilik = 1  
    fit = ash(ash.params$betahat,sebetahat.new,df=df.new,
                 method="fdr",mixcompdist=args$mixcompdist)
  }else{
    library(mixash)
    sebetahat.new = sqrt(va$PosteriorRate/va$PosteriorShape)
    df.new = 2*va$PosteriorShape
    pilik = va$PosteriorPi   
    fit = mixash(ash.params$betahat,sebetahat.new,df=df.new,pilik=pilik,
                 method="fdr",mixcompdist=args$mixcompdist)
  }
  
  qvalue = fit$qvalue
  # return estimates of variance, and q-values for testing beta=0
  return(list(fit=fit, va=va))
}
