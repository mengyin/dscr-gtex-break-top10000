# joint: vash+mixash/ash
# Suppose dispersions' are exchangeable and can use vash to shrink!
# se_g^2 = dispersion_g*scale_g, where scale_g is a fixed number for each gene.
# Using R package "vash" & "mixash"/"ash"
# input: betahat - estimate of effect size
#        sebetahat - estimated sd of betahat
jointash_disp.wrapper = function(input,args){
  library(vash)
  library(qvalue)
  
  if (args$transform=="quasibinom"){
    betahat = input$betahat.qb
    sebetahat = input$sebetahat.qb
    dispersion = input$dispersion.qb
    df = input$df.qb
  }else if (args$transform=="Myrna+quasibinom"){
    betahat = input$betahat.Myrnaqb
    sebetahat = input$sebetahat.Myrnaqb
    dispersion = input$dispersion.Myrnaqb
    df = input$df.Myrnaqb
  }else if (args$transform=="Myrnaoff+quasibinom"){
    betahat = input$betahat.Myrnaoffqb
    sebetahat = input$sebetahat.Myrnaoffqb
    dispersion = input$dispersion.Myrnaoffqb
    df = input$df.Myrnaoffqb
  }else if (args$transform=="RUV+quasibinom"){
    betahat = input$betahat.RUVqb
    sebetahat = input$sebetahat.RUVqb
    dispersion = input$dispersion.RUVqb
    df = input$df.RUVqb
  }else if (args$transform=="SVA+quasibinom"){
    betahat = input$betahat.SVAqb
    sebetahat = input$sebetahat.SVAqb
    dispersion = input$dispersion.RUVqb
    df = input$df.SVAqb
  }
  scale = sebetahat/sqrt(dispersion)
  
  va = vash(sqrt(dispersion),df=df,singlecomp=args$singlecomp)
  if (args$singlecomp==TRUE){
    library(ashr)
    sebetahat.new = scale*c(sqrt(va$PosteriorRate/va$PosteriorShape))
    df.new = 2*va$PosteriorShape[1]
    pilik = 1  
    fit = ash(betahat,sebetahat.new,df=df.new,
              method="fdr",mixcompdist=args$mixcompdist)
  }else{
    library(mixash)
    sebetahat.new = matrix(rep(scale,length(va$PosteriorPi)),ncol=length(va$PosteriorPi))*sqrt(va$PosteriorRate/va$PosteriorShape)
    df.new = 2*va$PosteriorShape
    pilik = va$PosteriorPi   
    fit = mixash(betahat,sebetahat.new,df=df.new,pilik=pilik,
                 method="fdr",mixcompdist=args$mixcompdist)
  }
  
  qvalue = fit$qvalue
  # return estimates of variance, and q-values for testing beta=0
  return(list(fit=fit, va=va))
}
