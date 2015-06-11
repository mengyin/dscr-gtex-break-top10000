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
  }else if (args$transform=="Myrna+quasibinom"){
    betahat = input$betahat.Myrnaqb
    sebetahat = input$sebetahat.Myrnaqb
    df = input$df.Myrnaqb
  }else if (args$transform=="Myrnaoff+quasibinom"){
    betahat = input$betahat.Myrnaoffqb
    sebetahat = input$sebetahat.Myrnaoffqb
    df = input$df.Myrnaoffqb
  }else if (args$transform=="RUV+quasibinom"){
    betahat = input$betahat.RUVqb
    sebetahat = input$sebetahat.RUVqb
    df = input$df.RUVqb
  }else if (args$transform=="SVA+quasibinom"){
    betahat = input$betahat.SVAqb
    sebetahat = input$sebetahat.SVAqb
    df = input$df.SVAqb
  }else if (args$transform=="edgeRglm"){
    betahat = input$betahat.edgeRglm
    sebetahat = input$sebetahat.edgeRglm
    df = input$df.edgeRglm
  }else if (args$transform=="DESeqglm"){
    betahat = input$betahat.DESeqglm
    sebetahat = input$sebetahat.DESeqglm
    df = input$df.DESeqglm
  }
  
  fit = ash(betahat,sebetahat,df=df,
            method="fdr",
            mixcompdist=args$mixcompdist)
  #qvalue = fit$qvalue
  # return estimates of beta, and q-values for testing beta=0
  #return(list(qvalue = qvalue, qvalue.fsr = qval.from.lfdr(fit$lfsr)))
  return(list(fit=fit))
}