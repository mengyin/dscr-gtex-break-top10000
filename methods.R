# Methods: 
sourceDir("methods")

addMethod(dsc_gtex,name="DESeq",fn=deseq.wrapper,outputtype="qval_output",args=list(exacttest=TRUE))
# addMethod(dsc_gtex,name="DESeq.glm",fn=deseq.wrapper,outputtype="qval_output",args=list(exacttest=FALSE))
addMethod(dsc_gtex,name="DESeq2",fn=deseq2.wrapper,outputtype="qval_output")

addMethod(dsc_gtex,name="edgeR",fn=edger.wrapper,outputtype="qval_output",args=list(exacttest=TRUE))
addMethod(dsc_gtex,name="edgeR.glm",fn=edger.wrapper,outputtype="qval_output",args=list(exacttest=FALSE))
addMethod(dsc_gtex,name="RUV+edgeR.glm",fn=edger.wrapper,outputtype="qval_output",
          args=list(exacttest=FALSE,RUV=TRUE))
addMethod(dsc_gtex,name="SVA+edgeR.glm",fn=edger.wrapper,outputtype="qval_output",
          args=list(exacttest=FALSE,SVA=TRUE))


addMethod(dsc_gtex,name="voom+limma",fn=limma.wrapper,outputtype="qval_output",
          args=list(transform="voom",robust=FALSE))
addMethod(dsc_gtex,name="voom+limmaR",fn=limma.wrapper,outputtype="qval_output",
          args=list(transform="voom",robust=TRUE))

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
addMethod(dsc_gtex,name="Myrnaoffqb+ash",fn=ash.wrapper,outputtype="ash_output",
          args=list(transform="Myrnaoff+quasibinom"))
addMethod(dsc_gtex,name="edgeRglm+ash",fn=ash.wrapper,outputtype="ash_output",
          args=list(transform="edgeRglm"))
addMethod(dsc_gtex,name="DESeqglm+ash",fn=ash.wrapper,outputtype="ash_output",
          args=list(transform="DESeqglm"))


addMethod(dsc_gtex,name="voom+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="voom",singlecomp=TRUE))
addMethod(dsc_gtex,name="voom+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="voom",singlecomp=FALSE))
addMethod(dsc_gtex,name="voom+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="voom",singlecomp=FALSE,mixcompdist="uniform"))

addMethod(dsc_gtex,name="qb+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="quasibinom",singlecomp=TRUE))
addMethod(dsc_gtex,name="qb+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="quasibinom",singlecomp=FALSE))
addMethod(dsc_gtex,name="qb+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="quasibinom",singlecomp=FALSE,mixcompdist="uniform"))

addMethod(dsc_gtex,name="RUVqb+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="RUV+quasibinom",singlecomp=TRUE))
# addMethod(dsc_gtex,name="RUVqb+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="RUV+quasibinom",singlecomp=FALSE))
# addMethod(dsc_gtex,name="RUVqb+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="RUV+quasibinom",singlecomp=FALSE,mixcompdist="uniform"))

addMethod(dsc_gtex,name="SVAqb+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="SVA+quasibinom",singlecomp=TRUE))
# addMethod(dsc_gtex,name="SVAqb+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="SVA+quasibinom",singlecomp=FALSE))
# addMethod(dsc_gtex,name="SVAqb+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="SVA+quasibinom",singlecomp=FALSE,mixcompdist="uniform"))

addMethod(dsc_gtex,name="Myrnaqb+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="Myrna+quasibinom",singlecomp=TRUE))
# addMethod(dsc_gtex,name="Myrnaqb+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="Myrna+quasibinom",singlecomp=FALSE))
# addMethod(dsc_gtex,name="Myrnaqb+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="Myrna+quasibinom",singlecomp=FALSE,mixcompdist="uniform"))

addMethod(dsc_gtex,name="qb+jointash_disp.single.u",fn=jointash_disp.wrapper,outputtype="jointash_output",
          args=list(transform="quasibinom",singlecomp=TRUE))
# addMethod(dsc_gtex,name="qb+jointash_disp.mix.n",fn=jointash_disp.wrapper,outputtype="jointash_output",
#           args=list(transform="quasibinom",singlecomp=FALSE,mixcompdist="normal"))
# addMethod(dsc_gtex,name="qb+jointash_disp.mix.u",fn=jointash_disp.wrapper,outputtype="jointash_output",
#           args=list(transform="quasibinom",singlecomp=FALSE,mixcompdist="uniform"))

addMethod(dsc_gtex,name="Myrnaoffqb+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="Myrnaoff+quasibinom",singlecomp=TRUE))
# addMethod(dsc_gtex,name="Myrnaoffqb+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="Myrnaoff+quasibinom",singlecomp=FALSE,mixcompdist="normal"))
# addMethod(dsc_gtex,name="Myrnaoffqb+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="Myrnaoff+quasibinom",singlecomp=FALSE,mixcompdist="uniform"))

addMethod(dsc_gtex,name="edgeRglm+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="edgeRglm",singlecomp=TRUE))
# addMethod(dsc_gtex,name="edgeRglm+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="edgeRglm",singlecomp=FALSE,mixcompdist="normal"))
# addMethod(dsc_gtex,name="edgeRglm+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="edgeRglm",singlecomp=FALSE,mixcompdist="uniform"))

addMethod(dsc_gtex,name="DESeqglm+jointash.single.u",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="DESeqglm",singlecomp=TRUE))
# addMethod(dsc_gtex,name="DESeqglm+jointash.mix.n",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="DESeqglm",singlecomp=FALSE,mixcompdist="normal"))
# addMethod(dsc_gtex,name="DESeqglm+jointash.mix.u",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="DESeqglm",singlecomp=FALSE,mixcompdist="uniform"))
