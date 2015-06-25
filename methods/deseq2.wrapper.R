deseq2.wrapper = function(input,args){
  library(DESeq2)
  library(qvalue)
  
  if(is.null(args$pseudocounts)){
    args$pseudocounts = 1
  }
  cond = input$condition
  dds = DESeqDataSetFromMatrix(input$counts+args$pseudocounts, DataFrame(cond), ~cond)
  dds = estimateSizeFactors(dds)
  dds = estimateDispersions(dds)
  dds = nbinomWaldTest(dds)
  res = results(dds)
  pvalue = res$pvalue
  qvalue = rep(NA,length(pvalue))
  qvalue[!is.na(pvalue)] = qvalue(pvalue[!is.na(pvalue)])$qval
  
  
  return(list(qvalue = qvalue, pvalue=pvalue))
}