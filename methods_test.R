# Methods: 
sourceDir("methods")

# addMethod(dsc_gtex,name="DESeq",fn=deseq.wrapper,outputtype="qval_output",args=list(exacttest=TRUE))
# addMethod(dsc_gtex,name="DESeq.glm",fn=deseq.wrapper,outputtype="qval_output",args=list(exacttest=FALSE))

addMethod(dsc_gtex,name="edgeR",fn=edger.wrapper,outputtype="qval_output",args=list(exacttest=TRUE))
addMethod(dsc_gtex,name="edgeR.glm",fn=edger.wrapper,outputtype="qval_output",args=list(exacttest=FALSE))
addMethod(dsc_gtex,name="RUV+edgeR.glm",fn=edger.wrapper,outputtype="qval_output",
          args=list(exacttest=FALSE,RUV=TRUE))
addMethod(dsc_gtex,name="SVA+edgeR.glm",fn=edger.wrapper,outputtype="qval_output",
          args=list(exacttest=FALSE,SVA=TRUE))


# addMethod(dsc_gtex,name="voom+limma",fn=limma.wrapper,outputtype="qval_output",
#           args=list(transform="voom",robust=FALSE))
# addMethod(dsc_gtex,name="voom+limmaR",fn=limma.wrapper,outputtype="qval_output",
#            args=list(transform="voom",robust=TRUE))

addMethod(dsc_gtex,name="voom+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="voom",singlecomp=TRUE))
#addMethod(dsc_gtex,name="voom+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="voom",singlecomp=FALSE))
addMethod(dsc_gtex,name="qb+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
           args=list(transform="quasibinom",singlecomp=TRUE))
#addMethod(dsc_gtex,name="qb+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="quasibinom",singlecomp=FALSE))
 
#addMethod(dsc_gtex,name="voom+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="voom",singlecomp=FALSE,mixcompdist="uniform"))
#addMethod(dsc_gtex,name="qb+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="quasibinom",singlecomp=FALSE,mixcompdist="uniform"))
 
addMethod(dsc_gtex,name="voom+ash",fn=ash.wrapper,outputtype="ash_output",
           args=list(transform="voom"))
addMethod(dsc_gtex,name="qb+ash",fn=ash.wrapper,outputtype="ash_output",
           args=list(transform="quasibinom"))
addMethod(dsc_gtex,name="RUVqb+ash",fn=ash.wrapper,outputtype="ash_output",
          args=list(transform="RUV+quasibinom"))
addMethod(dsc_gtex,name="SVAqb+ash",fn=ash.wrapper,outputtype="ash_output",
          args=list(transform="SVA+quasibinom"))
addMethod(dsc_gtex,name="Myrnaqb+ash",fn=ash.wrapper,outputtype="ash_output",
          args=list(transform="Myrna+quasibinom"))


addMethod(dsc_gtex,name="RUVqb+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="RUV+quasibinom",singlecomp=TRUE))
#addMethod(dsc_gtex,name="RUVqb+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="RUV+quasibinom",singlecomp=FALSE))
#addMethod(dsc_gtex,name="RUVqb+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="RUV+quasibinom",singlecomp=FALSE,mixcompdist="uniform"))

addMethod(dsc_gtex,name="SVAqb+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="SVA+quasibinom",singlecomp=TRUE))
#addMethod(dsc_gtex,name="SVAqb+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="SVA+quasibinom",singlecomp=FALSE))
#addMethod(dsc_gtex,name="SVAqb+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="SVA+quasibinom",singlecomp=FALSE,mixcompdist="uniform"))

addMethod(dsc_gtex,name="Myrnaqb+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="Myrna+quasibinom",singlecomp=TRUE))
#addMethod(dsc_gtex,name="Myrnaqb+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="Myrna+quasibinom",singlecomp=FALSE))
#addMethod(dsc_gtex,name="Myrnaqb+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="Myrna+quasibinom",singlecomp=FALSE,mixcompdist="uniform"))