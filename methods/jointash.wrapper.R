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

choose.ash.parameters <- function(input, args) {
  if (args$transform == "voom") {
    betahat <- input$betahat.voom
    sebetahat <- input$sebetahat.voom
    df <- input$df.voom
  } else if (args$transform == "quasibinom") {
    betahat <- input$betahat.qb
    sebetahat <- input$sebetahat.qb
    df <- input$df.qb
  } else if (args$transform == "Myrna+quasibinom") {
    betahat <- input$betahat.Myrnaqb
    sebetahat <- input$sebetahat.Myrnaqb
    df <- input$df.Myrnaqb
  } else if (args$transform == "Myrnaoff+quasibinom") {
    betahat <- input$betahat.Myrnaoffqb
    sebetahat <- input$sebetahat.Myrnaoffqb
    df <- input$df.Myrnaoffqb
  } else if (args$transform == "RUV+quasibinom") {
    betahat <- input$betahat.RUVqb
    sebetahat <- input$sebetahat.RUVqb
    df <- input$df.RUVqb
  } else if (args$transform == "SVA+quasibinom") {
    betahat <- input$betahat.SVAqb
    sebetahat <- input$sebetahat.SVAqb
    df <- input$df.SVAqb
  } else if (args$transform == "edgeRglm") {
    betahat <- input$betahat.edgeRglm
    sebetahat <- input$sebetahat.edgeRglm
    df <- input$df.edgeRglm
  } else if (args$transform == "DESeqglm") {
    betahat <- input$betahat.DESeqglm
    sebetahat <- input$sebetahat.DESeqglm
    df <- input$df.DESeqglm
  }
  
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
} 
