datamaker = function(args){  
  tissue = args$tissue
  
  rawdata1 = read.table(paste0(path,"/dscr-gtex/",tissue[1],".txt"),header=TRUE)
    
  Nsamp = args$Nsamp # Number of samples in each of the 2 groups
  
  if (is.null(args$poisthin)){
    args$poisthin = FALSE
  }
  if (is.null(args$breaksample)){
    args$breaksample = FALSE
  }
  if (is.null(args$nullpi)){
    if (args$poisthin==TRUE){
      args$nullpi = 0.9
    }else if (length(args$tissue)==1){
      args$nullpi = 1
    }else if (length(args$tissue)>1){
      args$nullpi = 0
    }
  }
  
  if (length(tissue)>1){
    rawdata2 = read.table(paste0(path,"/dscr-gtex/",tissue[2],".txt"),header=TRUE)
    if (dim(rawdata1)[2]<Nsamp | dim(rawdata2)[2]<Nsamp){
      stop("Not enough samples in the raw dataset!")
    }
    
    if (args$nullpi==0){
      # All genes are alternatives
      temp1 = selectsample(rawdata1, Nsamp, args$breaksample)
      counts1 = temp1$counts
      subsample1 = temp1$subsample
      rm(temp1)      
      temp2 = selectsample(rawdata2, Nsamp, args$breaksample)
      counts2 = temp2$counts
      subsample2 = temp2$subsample
      rm(temp2)
      
      counts = cbind(counts1, counts2)
      subsample = cbind(subsample1, subsample2)
    }else{
      # Some genes are nulls, some are alternatives
      temp1 = selectsample(rawdata1, 2*Nsamp, args$breaksample)
      counts1 = temp1$counts
      subsample1 = temp1$subsample
      rm(temp1)      
      temp2 = selectsample(rawdata2, Nsamp, args$breaksample)
      counts2 = temp2$counts
      subsample2 = temp2$subsample
      rm(temp2)
      counts = cbind(counts1, counts2)
      subsample = cbind(subsample1, subsample2)
    }
  }else{
    if (dim(rawdata1)[2] < 2*Nsamp){
      stop("Not enough samples in the raw dataset!")
    }
    
    temp = selectsample(rawdata1, 2*Nsamp, args$breaksample)
    counts = temp$counts
    subsample = temp$subsample
    rm(temp)
  } 
  
  # Remove genes without any reads
  subsample = subsample[apply(counts,1,sum)>0,]
  counts = counts[apply(counts,1,sum)>0,]
  
  # Take the top Ngene high-expressed genes
  Ngene = args$Ngene 
  if (!is.null(Ngene)){
    Ngene = min(Ngene, dim(counts)[1])
    subsample = subsample[sort(order(rowSums(counts),decreasing=TRUE)[1:Ngene]),]
    counts = counts[sort(order(rowSums(counts),decreasing=TRUE)[1:Ngene]),]
  }
  Ngene = dim(counts)[1]
  
  # Model's design: Nsamp samples for group A and Nsamp samples for group B
  condition = factor(rep(1:2,each=Nsamp))
  design = model.matrix(~condition)
  
  # Ground truth of null hypotheses: beta_g=0
  null = rep(0,Ngene)
  null[sample(Ngene, round(Ngene*args$nullpi))] = 1  
  
  # Poison thinning
  if (args$poisthin==TRUE){
    if (is.null(args$log2foldmean)){
      args$log2foldmean = 0
    }
    if (is.null(args$log2foldsd)){
      args$log2foldsd = 1
    }
     
    log2foldchanges = rnorm(sum(!null), mean=args$log2foldmean, sd=args$log2foldsd)
    foldchanges = 2^log2foldchanges
    
    # thin group A
    counts[which(!null)[log2foldchanges>0],1:Nsamp]=matrix(rbinom(sum(log2foldchanges>0)*Nsamp, 
                                                      size=c(as.matrix(counts[which(!null)[log2foldchanges>0],1:Nsamp])),
                                                      prob=rep(1/foldchanges[log2foldchanges>0],Nsamp)),ncol=Nsamp)
    # thin group B
    counts[which(!null)[log2foldchanges<0],(Nsamp+1):(2*Nsamp)]=matrix(rbinom(sum(log2foldchanges<0)*Nsamp, 
                                                                  size=c(as.matrix(counts[which(!null)[log2foldchanges<0],(Nsamp+1):(2*Nsamp)])),
                                                                  prob=rep(foldchanges[log2foldchanges<0],Nsamp)),ncol=Nsamp)
    
  }else if(args$nullpi<1 & args$nullpi>0){
    newcounts = matrix(rep(0, Ngene*2*Nsamp),nrow=Ngene)
    newcounts[as.logical(null),] = counts[as.logical(null),1:(2*Nsamp)]
    newcounts[!null,] = counts[!null,c(1:Nsamp,(2*Nsamp+1):(3*Nsamp))]
    counts = newcounts
    newsubsample = matrix(rep(0, Ngene*2*Nsamp),nrow=Ngene)
    newsubsample[as.logical(null),] = subsample[as.logical(null),1:(2*Nsamp)]
    newsubsample[!null,] = subsample[!null,c(1:Nsamp,(2*Nsamp+1):(3*Nsamp))]
    subsample = newsubsample
    rm(newcounts); rm(newsubsample);
  }
  
  
  # Voom transformation
  library(edgeR)
  dgecounts = DGEList(counts=counts,group=condition)
  if (args$voom.normalize == TRUE){
      dgecounts = calcNormFactors(dgecounts)
  }
  library(limma)
  v = voom(dgecounts,design,plot=FALSE)
  zdat.voom = apply(cbind(v$E,v$weights),1,wls.wrapper,g=condition)
  betahat.voom = zdat.voom[1,]
  sebetahat.voom = zdat.voom[2,]
  df.voom = length(condition)-2
  
  # Quasi-binomial glm
  zdat.qb = counts.associate(counts, condition)
  betahat.qb = zdat.qb[3,]
  sebetahat.qb = zdat.qb[4,]
  dispersion.qb = zdat.qb[5,]
  df.qb = length(condition)-2
  
  # Myrna & Quasi-binomial glm
  # Use log(75th quantile of samples' counts) as covariate
  W.Myrna = apply(counts,2,function(x) log(quantile(x[x>0],0.75)))
  
  zdat.Myrnaqb = counts.associate(counts, condition, W=W.Myrna)
  betahat.Myrnaqb = zdat.Myrnaqb[3,]
  sebetahat.Myrnaqb = zdat.Myrnaqb[4,]
  dispersion.Myrnaqb = zdat.Myrnaqb[5,]
  df.Myrnaqb = length(condition)-2-1
  # Use log(75th quantile of samples' counts) as offset
  offset.Myrnaoff = apply(counts,2,function(x) quantile(x[x>0],0.75))
  zdat.Myrnaoffqb = counts.associate(counts, condition, W=NULL, offset=offset.Myrnaoff)
  betahat.Myrnaoffqb = zdat.Myrnaoffqb[3,]
  sebetahat.Myrnaoffqb = zdat.Myrnaoffqb[4,]
  dispersion.Myrnaoffqb = zdat.Myrnaoffqb[5,]
  df.Myrnaoffqb = length(condition)-2
  
  # RUV & quasi-binomial glm
  library(RUVSeq)
  seq = newSeqExpressionSet(as.matrix(counts))
  controls = rownames(seq)
  differences = matrix(data=c(1:Nsamp, (Nsamp+1):(2*Nsamp)), byrow=TRUE, nrow=2)
  if (is.null(args$RUV.k)){
    args$RUV.k = round(log2(Nsamp))
  }
  seqRUVs = RUVs(seq, controls, k=args$RUV.k, differences)
  zdat.RUVqb = counts.associate(counts, condition, W=pData(seqRUVs)$W_1)
  betahat.RUVqb = zdat.RUVqb[3,]
  sebetahat.RUVqb = zdat.RUVqb[4,]
  dispersion.RUVqb = zdat.RUVqb[5,]
  df.RUVqb = length(condition)-2-args$RUV.k
  
  # SVA & quasi-binomial glm
  library(sva)
  mod1 = model.matrix(~condition)
  mod0 = cbind(mod1[,1])
  
  if (args$nullpi>0){
    svseq = svaseq(counts,mod1,mod0,control=null)
  }else{
    svseq = svaseq(counts,mod1,mod0)
  }
  
  if(svseq$n.sv>0){
    zdat.SVAqb = counts.associate(counts, condition, W=svseq$sv)
  }else{
    zdat.SVAqb = zdat.qb
  }
  betahat.SVAqb = zdat.SVAqb[3,]
  sebetahat.SVAqb = zdat.SVAqb[4,]
  dispersion.SVAqb = zdat.SVAqb[5,]
  df.SVAqb = length(condition)-2-svseq$n.sv
  
  # Get sebetahat from edgeR.glm (infer from betahat & pval)
  if (is.null(args$pseudocounts)){
    args$pseudocounts = 1
  }
  y = DGEList(counts=counts+args$pseudocounts, group=condition)
  y = calcNormFactors(y)
  y = estimateGLMCommonDisp(y,design)
  y = estimateGLMTrendedDisp(y,design)
  y = estimateGLMTagwiseDisp(y,design)
  fit = glmFit(y,design)
  lrt = glmLRT(fit,coef=2)
  betahat.edgeRglm = fit$coef[,2]
  df.edgeRglm = length(condition)-2
  tscore = qt(1-lrt$table$PValue/2,df=df.edgeRglm)
  sebetahat.edgeRglm = abs(betahat.edgeRglm/tscore)
  
  # Get sebetahat from DESeq.glm (infer from betahat & pval)
  library(DESeq)
  cds = newCountDataSet(counts+args$pseudocounts, condition )
  cds = estimateSizeFactors( cds )
  cds = try(estimateDispersions( cds ),silent=TRUE)
  if (class(cds)=="try-error"){
    betahat.DESeqglm = NA
    sebetahat.DESeqglm = NA
    df.DESeqglm = length(condition)-2
  }else{
    fit1 = fitNbinomGLMs( cds, count ~ condition )     
    fit0 = fitNbinomGLMs( cds, count ~ 1 )
    betahat.DESeqglm = fit1[,2]
    df.DESeqglm = length(condition)-2
    tscore = qt(1-nbinomGLMTest(fit1,fit0)/2,df=df.DESeqglm)
    sebetahat.DESeqglm = abs(betahat.DESeqglm/tscore)
  }
  
  
  
  # meta data
  meta = list(tissue=tissue, subsample=subsample, null=null, 
              args=args)
  
  # input data
  input = list(counts=counts, condition=condition,
               v=v, betahat.voom=betahat.voom, sebetahat.voom=sebetahat.voom, df.voom=df.voom,
               betahat.qb=betahat.qb, sebetahat.qb=sebetahat.qb, df.qb=df.qb, dispersion.qb=dispersion.qb,
               betahat.RUVqb=betahat.RUVqb, sebetahat.RUVqb=sebetahat.RUVqb, dispersion.RUVqb=dispersion.RUVqb, df.RUVqb=df.RUVqb, W.RUV=pData(seqRUVs)$W_1,
               betahat.SVAqb=betahat.SVAqb, sebetahat.SVAqb=sebetahat.SVAqb, dispersion.SVAqb=dispersion.SVAqb, df.SVAqb=df.SVAqb, W.SVA=svseq$sv,
               betahat.Myrnaqb=betahat.Myrnaqb, sebetahat.Myrnaqb=sebetahat.Myrnaqb, dispersion.Myrnaqb=dispersion.Myrnaqb, df.Myrnaqb=df.Myrnaqb, W.Myrna=W.Myrna,
               betahat.Myrnaoffqb=betahat.Myrnaoffqb, sebetahat.Myrnaoffqb=sebetahat.Myrnaoffqb, dispersion.Myrnaoffqb=dispersion.Myrnaoffqb, df.Myrnaoffqb=df.Myrnaoffqb, offset.Myrnaoff=offset.Myrnaoff,
               betahat.edgeRglm=betahat.edgeRglm, sebetahat.edgeRglm=sebetahat.edgeRglm,df.edgeRglm=df.edgeRglm,
               betahat.DESeqglm=betahat.DESeqglm, sebetahat.DESeqglm=sebetahat.DESeqglm,df.DESeqglm=df.DESeqglm)
  
  data = list(meta=meta,input=input)
  return(data)
}

# Weighted least squares regression
# g: formula
# ynweights: matrix of response y and corresponding weights
wls.wrapper = function(ynweights,g,...){
  y = ynweights[1:(length(ynweights)/2)]
  weights = ynweights[(length(ynweights)/2+1):length(ynweights)]
  y.wls = lm(y~g,weights=weights,...)
  
  # slope estimate & standard error
  c = as.vector(t(summary(y.wls)$coeff[2,1:2]))
  return(c)
}

# counts is a ngene (or nwindow) by nsample matrix of counts (eg RNAseq)
# g is an nsample vector of group memberships/covariate
# looks for association between rows and covariate
counts.associate=function(counts, g,pseudocounts=1,W=NULL,offset=NULL){
  y.counts=t(as.matrix(counts)) 
  col.sum = apply(y.counts, 2, sum)
  y.counts=y.counts[,col.sum>0] #remove 0 columns
  y.counts = y.counts+pseudocounts 
  if (is.null(offset)){
    offset = apply(y.counts,1,sum)
  }
  
  y.prop = y.counts/apply(y.counts,1,sum) # compute proportions
  zdat = apply(y.prop,2,glm.binomial.wrapper,g=g,W=W,weights=offset,epsilon=1e-6)  #estimate effect sizes and standard errors
  #zdat.ash = ash(zdat[3,],zdat[4,],df=2,method='fdr') #shrink the estimated effects
  #return(list(zdat=zdat,zdat.ash=zdat.ash))
  return(zdat)
}

glm.binomial.wrapper = function(y,g,W=NULL,...){
  if (is.null(W)){
    y.glm=safe.quasibinomial.glm(y~g,...)
  }else{
    y.glm=safe.quasibinomial.glm(y~g+W,...)
  }
  
  return(c(get.coeff(y.glm),summary(y.glm)$dispersion))
}


#fill NAs with 0s (or other value)
fill.nas=function(x,t=0){
  x[is.na(x)]=t
  return(x)
}

#get estimates and standard errors from a glm object
#return NAs if not converged
get.coeff=function(x.glm){
  c=as.vector(t(summary(x.glm)$coeff[,1:2]))  
  if(x.glm$conv){return(c)} else {return(rep(NA,length(c)))}
}

# use glm to fit quasibinomial, but don't allow for underdispersion!
safe.quasibinomial.glm=function(formula,forcebin=FALSE,...){
  if(forcebin){
    fm=glm(formula,family=binomial,...)
  } else{
    fm = glm(formula,family=quasibinomial,...)
    if(is.na(summary(fm)$dispersion) | summary(fm)$dispersion<1){
      fm=glm(formula,family=binomial,...)
    }
  }
  return(fm)
}

# randomly subsample data for each gene
# gene: a vector of reads for one gene
# Nsamp: # of samples wanted
sampleingene = function(gene, Nsamp){
  sample = sample(length(gene),Nsamp)
  return(c(gene[sample],sample))
}

# randomly select samples
# counts: full count matrix
# Nsamp: # of samples wanted
# breaksample: flag, if select different samples for each gene
selectsample = function(counts, Nsamp, breaksample){
  if (breaksample==FALSE){
    subsample = sample(1:dim(counts)[2],Nsamp)
    counts = counts[,subsample]
    subsample = t(matrix(rep(subsample, dim(counts)[1]), ncol=dim(counts)[1]))
  }else{
    temp = t(apply(counts, 1, sampleingene, Nsamp=Nsamp))
    counts = temp[,1:Nsamp]
    subsample = temp[,(Nsamp+1):(2*Nsamp)]
  }
  return(list(counts=counts, subsample=subsample))
}
