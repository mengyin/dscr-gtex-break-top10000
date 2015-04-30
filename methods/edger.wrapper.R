edger.wrapper = function(input,args){
  library(edgeR)
  library(qvalue)
  if(args$exacttest==TRUE){
    y = DGEList(counts=input$counts,group=input$condition)
    y = calcNormFactors(y)
    y = estimateCommonDisp(y)
    y = estimateTagwiseDisp(y)
    et = exactTest(y)
    pvalue = et$table$PValue
    qvalue = qvalue(pvalue)$qval
  }else{
    if(is.null(args$pseudocounts)){
      args$pseudocounts = 1
    }
    if(is.null(args$RUV)){
      args$RUV = FALSE
    }
    if(is.null(args$SVA)){
      args$SVA = FALSE
    }
    if(args$RUV==TRUE){
      design = model.matrix(~input$condition+input$W.RUV)
    }else if(args$SVA==TRUE){
      if (sum(input$W.SVA!=0)){
        design = model.matrix(~input$condition+input$W.SVA) 
      }else{
        design = model.matrix(~input$condition) 
      }        
    }else{
      design = model.matrix(~input$condition)
    }
    
    y = DGEList(counts=input$counts+args$pseudocounts, group=input$condition)
    y = calcNormFactors(y)
    y = estimateGLMCommonDisp(y,design)
    y = estimateGLMTrendedDisp(y,design)
    y = estimateGLMTagwiseDisp(y,design)
    fit = glmFit(y,design)
    lrt = glmLRT(fit,coef=2)
    pvalue = lrt$table$PValue
    qvalue = qvalue(pvalue)$qval
  }
  
  return(list(qvalue = qvalue, pvalue=pvalue))
}
