# input: betahat - estimate of effect size
#        sebetahat - estimated sd of betahat
ash.wrapper = function(input,args){
  library(ashr)
  library(qvalue)

  if (args$transform=="voom"){
    betahat = input$betahat.voom
    sebetahat = input$sebetahat.voom
    df = input$df.voom
  }else if (args$transform=="quasibinom"){
    betahat = input$betahat.qb
    sebetahat = input$sebetahat.qb
    df = input$df.qb
  }else if (args$transform=="RUV+quasibinom"){
    betahat = input$betahat.RUVqb
    sebetahat = input$sebetahat.RUVqb
    df = input$df.RUVqb
  }else if (args$transform=="SVA+quasibinom"){
    betahat = input$betahat.SVAqb
    sebetahat = input$sebetahat.SVAqb
    df = input$df.SVAqb
  }
  
  fit = ash(betahat,sebetahat,df=df,
            method="fdr",
            mixcompdist=args$mixcompdist)
  #qvalue = fit$qvalue
  # return estimates of beta, and q-values for testing beta=0
  #return(list(qvalue = qvalue, qvalue.fsr = qval.from.lfdr(fit$lfsr)))
  return(list(fit=fit))
}