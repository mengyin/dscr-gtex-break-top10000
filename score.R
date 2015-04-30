#This file should define your score function

score = function(data, output){
  if (class(output)=="list"){
    qvalue = output$qvalue
    FDR = mean(qvalue<=0.05)
    if (is.null(output$qvalue.lfsr)){
      FDR.lfsr = FDR
    }else{
      FDR.lfsr = mean(output$qvalue.lfsr<=0.05)
    }
    res = c(FDR, FDR.lfsr)
    names(res) = c("FDR_005","FDR.fsr_005")
  }else{
    res = c(NA, NA)
    names(res) = c("FDR_005","FDR.fsr_005")
  }
  return(res)
}
