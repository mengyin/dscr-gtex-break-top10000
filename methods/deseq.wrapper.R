deseq.wrapper = function(input,args){
  library(DESeq)
  library(qvalue)

  if(args$exacttest==TRUE){
    cds = newCountDataSet( input$counts, input$condition )
    cds = estimateSizeFactors( cds )
    if (length(input$condition)<5){
      cds = estimateDispersions( cds, method="pooled")
    }else{
      cds = estimateDispersions( cds )
    }
    pvalue = nbinomTest( cds, 1, 2,  pvals_only = TRUE)
  }else{
    if(is.null(args$pseudocounts)){
      args$pseudocounts = 1
    }
    cds = newCountDataSet(input$counts+args$pseudocounts, input$condition )
    cds = estimateSizeFactors( cds )
    cds = estimateDispersions( cds )
    fit1 = fitNbinomGLMs( cds, count ~ condition )     
    fit0 = fitNbinomGLMs( cds, count ~ 1 )
    pvalue = nbinomGLMTest( fit1, fit0 )
  }
  qvalue = qvalue(pvalue)$qval

  return(list(qvalue = qvalue, pvalue=pvalue))
}