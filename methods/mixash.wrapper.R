# vash+mixash
# Using R package "vash" & "mixash"
# input: betahat - estimate of effect size
#        sebetahat - estimated sd of betahat
mixash.wrapper = function(input,args){
  library(mixash)
  library(vash)
  library(qvalue)
  va = vash(input$sebetahat,df=input$df,singlecomp=args$singlecomp)
  
  sebetahat.mix = sqrt(va$PosteriorRate/va$PosteriorShape)
  df = 2*va$PosteriorShape
  pilik = va$PosteriorPi
  
  fit = mixash(input$betahat,sebetahat.mix,df=df,pilik=pilik,
               method="fdr",
               mixcompdist=args$mixcompdist)
  qvalue = fit$qvalue
  # return estimates of variance, and q-values for testing beta=0
  return(list(mean.est = fit$PosteriorMean,
              qvalue = qvalue))
}