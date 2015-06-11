sourceDir("datamakers")
fileName = 'data_path.txt'
path = gsub("[[:space:]]", "", readChar(fileName, file.info(fileName)$size))

addScenario(dsc_gtex,name="Adipose-Subcutaneous,2vs2", fn=datamaker,
            args=list(tissue="Adipose-Subcutaneous", 
                      Nsamp=2, Ngene=10000, breaksample=TRUE, path=path,
                      voom.normalize=TRUE),
            seed=1)

# addScenario(dsc_gtex,name="Adipose-Subcutaneous,10vs10", fn=datamaker,
#             args=list(tissue="Adipose-Subcutaneous", 
#                       Nsamp=10, Ngene=10000, breaksample=TRUE, path=path,
#                       voom.normalize=TRUE),
#             seed=1:50)
# 
# addScenario(dsc_gtex,name="Adipose-Subcutaneous,50vs50", fn=datamaker,
#             args=list(tissue="Adipose-Subcutaneous",
#                       Nsamp=50, Ngene=10000, breaksample=TRUE, path=path,
#                       voom.normalize=TRUE),
#             seed=1:50)
