library(dscr)

dsc_gtex = new.dsc("gtex","dsc-gtex-files")
source("scenarios_test.R")
source("methods_test.R")
source("score.R")

jointash2qval_est =function(output){
  library(ashr)
  if (class(output)=="list"){
    qvalue=output$fit$qvalue
    qvalue.fsr = qval.from.lfdr(output$fit$lfsr)
    #qvalue.fsra = qval.from.lfdr(output$fit$lfsra)
    return(list(qvalue=qvalue, qvalue.fsr=qvalue.fsr))
  }else{
    return(list(qvalue=NA, qvalue.fsr=NA))
  }
} 
addOutputParser(dsc_gtex,"jointash2qval",jointash2qval_est,"jointash_output","qval_output")
addOutputParser(dsc_gtex,"ash2qval",jointash2qval_est,"ash_output","qval_output")

addScore(dsc_gtex,score,name="score",outputtype="qval_output")

res=run_dsc(dsc_gtex)

# save(res,file="res.Rdata")

# reg = dsc2BE(dsc_gtex, "BE_gtex", 5)

