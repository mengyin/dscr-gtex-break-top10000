M.list = listMethods(dsc_gtex)
S.list = listScenarios(dsc_gtex)
Onames=getOutputtypes(dsc_gtex)

# Check initial output
for (i in 1:length(S.list)){
  for (j in 1:length(M.list)){
    for (k in 1:50){
      filename = paste0("./dsc-gtex-files/output/",S.list[i],"/",
                        M.list[j],"/",Onames[[M.list[j]]],"/output.",k,".rds")
      temp = try(readRDS(filename))
      if (class(temp)=="try-error"){
        print(filename)
        #file.remove(filename)
      }
    }   
  }
}

# Check qval_output
for (i in 1:length(S.list)){
  for (j in 1:length(M.list)){
    for (k in 1:50){
      filename = paste0("./dsc-gtex-files/output/",S.list[i],"/",
                        M.list[j],"/qval_output/output.",k,".rds")
      temp = try(readRDS(filename))
      if (class(temp)=="try-error"){
        print(filename)
        file.remove(filename)
      }
    }   
  }
}