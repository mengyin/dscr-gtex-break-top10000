library(edgeR)
library(limma)
library(RUVSeq)
library(sva)
library(DESeq)

datamaker = function(args){  
  dfargs = default_datamaker_args(args)
  rawdata1 = read.table(paste0(path,"/dscr-gtex/",dfargs$tissue[1],".txt"),header=TRUE)
  
  if (length(dfargs$tissue)>1){
    rawdata2 = read.table(paste0(path,"/dscr-gtex/",dfargs$tissue[2],".txt"),header=TRUE)
    if (dim(rawdata1)[2]<dfargs$Nsamp | dim(rawdata2)[2]<dfargs$Nsamp){
      stop("Not enough samples in the raw dataset!")
    }
    
    if (dfargs$nullpi==0){
      # All genes are alternatives
      temp1 = selectsample(rawdata1, dfargs$Nsamp, dfargs$breaksample)
      counts1 = temp1$counts
      subsample1 = temp1$subsample
      rm(temp1)      
      temp2 = selectsample(rawdata2, dfargs$Nsamp, dfargs$breaksample)
      counts2 = temp2$counts
      subsample2 = temp2$subsample
      rm(temp2)
      
      counts = cbind(counts1, counts2)
      subsample = cbind(subsample1, subsample2)
    }else{
      # Some genes are nulls, some are alternatives
      temp1 = selectsample(rawdata1, 2*dfargs$Nsamp, dfargs$breaksample)
      counts1 = temp1$counts
      subsample1 = temp1$subsample
      rm(temp1)      
      temp2 = selectsample(rawdata2, dfargs$Nsamp, dfargs$breaksample)
      counts2 = temp2$counts
      subsample2 = temp2$subsample
      rm(temp2)
      counts = cbind(counts1, counts2)
      subsample = cbind(subsample1, subsample2)
    }
  }else{
    if (dim(rawdata1)[2] < 2*dfargs$Nsamp){
      stop("Not enough samples in the raw dataset!")
    }
    
    temp = selectsample(rawdata1, 2*dfargs$Nsamp, dfargs$breaksample)
    counts = temp$counts
    subsample = temp$subsample
    rm(temp)
  } 
  
  # Remove genes without any reads
  subsample = subsample[apply(counts,1,sum)>0,]
  counts = counts[apply(counts,1,sum)>0,]
  
  # Take the top Ngene high-expressed genes
  if (!is.null(dfargs$Ngene)){
    dfargs$Ngene = min(dfargs$Ngene, dim(counts)[1])
    subsample = subsample[sort(order(rowSums(counts),decreasing=TRUE)[1:dfargs$Ngene]),]
    counts = counts[sort(order(rowSums(counts),decreasing=TRUE)[1:dfargs$Ngene]),]
  }
  dfargs$Ngene = dim(counts)[1]
  
  # Model's design: Nsamp samples for group A and Nsamp samples for group B
  condition = factor(rep(1:2,each=dfargs$Nsamp))
  design = model.matrix(~condition)
  
  # Ground truth of null hypotheses: beta_g=0
  null = rep(0,dfargs$Ngene)
  null[sample(dfargs$Ngene, round(dfargs$Ngene*dfargs$nullpi))] = 1  
  
  # Poisson thinning (optional)
  counts = pois_thinning(counts, dfargs, null)
  
  # Mix null and alternative genes from different samples (optional)
  mix = mix_sample(counts, dfargs, null, subsample)
  counts = mix$counts
  subsample = mix$subsample
  
  # Voom transformation
  voom = voom_transform(counts, condition)
  
  # Quasi-binomial glm
  qb = quasi_binom(counts, condition)
  
  # Myrna & Quasi-binomial glm
  # Use log(75th quantile of samples' counts) as covariate
  W.Myrna = apply(counts,2,function(x) log(quantile(x[x>0],0.75)))
  Myrnaqb = quasi_binom(counts, condition, W=W.Myrna)
  # Use log(75th quantile of samples' counts) as offset
  offset.Myrnaoff = apply(counts,2,function(x) quantile(x[x>0],0.75))
  Myrnaoffqb = quasi_binom(counts, condition, W=NULL, offset=offset.Myrnaoff)
  
  # RUV & quasi-binomial glm
  W.RUV = RUV_factor(counts, dfargs, null)
  RUVqb = quasi_binom(counts, condition, W=W.RUV)
  
  # SVA & quasi-binomial glm
  W.SVA = SVA_factor(counts, condition, dfargs, null)
  SVAqb = quasi_binom(counts, condition, W=W.SVA)
  
  # Get sebetahat from edgeR.glm (infer from betahat & pval)
  edgeRglm = edgeR_glmest(counts, condition, dfargs)
  
  # Get sebetahat from DESeq.glm (infer from betahat & pval)
  DESeqglm = DESeq_glmest(counts, condition, dfargs)
  
  # meta data
  meta = list(subsample=subsample, null=null, 
              dfargs=dfargs)
  
  # input data
  input = list(counts=counts, condition=condition,
               v=voom$v, betahat.voom=voom$betahat, sebetahat.voom=voom$sebetahat, df.voom=voom$df,
               betahat.qb=qb$betahat, sebetahat.qb=qb$sebetahat, df.qb=qb$df, dispersion.qb=qb$dispersion,
               betahat.RUVqb=RUVqb$betahat, sebetahat.RUVqb=RUVqb$sebetahat, dispersion.RUVqb=RUVqb$dispersion, df.RUVqb=RUVqb$df, W.RUV=W.RUV,
               betahat.SVAqb=SVAqb$betahat, sebetahat.SVAqb=SVAqb$sebetahat, dispersion.SVAqb=SVAqb$dispersion, df.SVAqb=SVAqb$df, W.SVA=W.SVA,
               betahat.Myrnaqb=Myrnaqb$betahat, sebetahat.Myrnaqb=Myrnaqb$sebetahat, dispersion.Myrnaqb=Myrnaqb$dispersion, df.Myrnaqb=Myrnaqb$df, W.Myrna=W.Myrna,
               betahat.Myrnaoffqb=Myrnaoffqb$betahat, sebetahat.Myrnaoffqb=Myrnaoffqb$sebetahat, dispersion.Myrnaoffqb=Myrnaoffqb$dispersion, df.Myrnaoffqb=Myrnaoffqb$df, offset.Myrnaoff=offset.Myrnaoff,
               betahat.edgeRglm=edgeRglm$betahat, sebetahat.edgeRglm=edgeRglm$sebetahat, df.edgeRglm=edgeRglm$df,
               betahat.DESeqglm=DESeqglm$betahat, sebetahat.DESeqglm=DESeqglm$sebetahat, df.DESeqglm=DESeqglm$df)
  
  data = list(meta=meta,input=input)
  return(data)
}

# Set default arguments for datamaker function
default_datamaker_args = function(args){
  # poisthin: flag of Poisson thinning
  if (is.null(args$poisthin)){
    args$poisthin = FALSE
  }
  
  # log2foldmean, log2foldsd: Poisson thinning params
  if (args$poisthin==TRUE){
    if (is.null(args$log2foldmean)){
      args$log2foldmean = 0
    }
    if (is.null(args$log2foldsd)){
      args$log2foldsd = 1
    }
  }
  
  # breaksample: flag of each gene randomly select samples
  if (is.null(args$breaksample)){
    args$breaksample = FALSE
  }
  
  
  # nullpi: proportion of null genes
  if (is.null(args$nullpi)){
    if (args$poisthin==TRUE){
      args$nullpi = 0.9
    }else if (length(args$tissue)==1){
      args$nullpi = 1
    }else if (length(args$tissue)>1){
      args$nullpi = 0
    }
  }
  
  # RUV.k: number of surrogate variables for SVA
  if (is.null(args$RUV.k)){
    args$RUV.k = round(log2(args$Nsamp))
  }
  
  # pseudocounts: add pseudocounts to count matrix
  if (is.null(args$pseudocounts)){
    args$pseudocounts = 1
  }
  
  return(args)
}

# Poisson thinning
pois_thinning = function(counts, args, null){
  if (args$poisthin==TRUE){ 
    log2foldchanges = rnorm(sum(!null), mean=args$log2foldmean, sd=args$log2foldsd)
    foldchanges = 2^log2foldchanges
    
    # thin group A
    counts[which(!null)[log2foldchanges>0],1:args$Nsamp]=matrix(rbinom(sum(log2foldchanges>0)*args$Nsamp, 
                                                                       size=c(as.matrix(counts[which(!null)[log2foldchanges>0],1:args$Nsamp])),
                                                                       prob=rep(1/foldchanges[log2foldchanges>0],Nsamp)),ncol=args$Nsamp)
    # thin group B
    counts[which(!null)[log2foldchanges<0],(args$Nsamp+1):(2*args$Nsamp)]=matrix(rbinom(sum(log2foldchanges<0)*args$Nsamp, 
                                                                                        size=c(as.matrix(counts[which(!null)[log2foldchanges<0],
                                                                                                                (args$Nsamp+1):(2*args$Nsamp)])),
                                                                                        prob=rep(foldchanges[log2foldchanges<0],args$Nsamp)),
                                                                                 ncol=args$Nsamp)
    
  }
  return(counts)
}

# Mix null and alternative genes from different samples
mix_sample = function(counts, args, null, subsample){
  if(args$nullpi<1 & args$nullpi>0){
    newcounts = matrix(rep(0, args$Ngene*2*args$Nsamp),nrow=args$Ngene)
    newcounts[as.logical(null),] = counts[as.logical(null),1:(2*args$Nsamp)]
    newcounts[!null,] = counts[!null,c(1:args$Nsamp,(2*args$Nsamp+1):(3*args$Nsamp))]
    counts = newcounts
    newsubsample = matrix(rep(0, args$Ngene*2*args$Nsamp),nrow=args$Ngene)
    newsubsample[as.logical(null),] = subsample[as.logical(null),1:(2*args$Nsamp)]
    newsubsample[!null,] = subsample[!null,c(1:args$Nsamp,(2*args$Nsamp+1):(3*args$Nsamp))]
    subsample = newsubsample
    rm(newcounts); rm(newsubsample);
  }
  return(list(counts=counts, subsample=subsample))
}

# Voom transformation
voom_transform = function(counts, condition){
  
  dgecounts = calcNormFactors(DGEList(counts=counts,group=condition))
  
  design = model.matrix(~condition)
  v = voom(dgecounts,design,plot=FALSE)
  zdat.voom = apply(cbind(v$E,v$weights),1,wls.wrapper,g=condition)
  betahat.voom = zdat.voom[1,]
  sebetahat.voom = zdat.voom[2,]
  df.voom = length(condition)-2
  
  return(list(betahat=betahat.voom, sebetahat=sebetahat.voom, df=df.voom, v=v))
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

# Quasi-binomial glm
quasi_binom = function(counts, condition, W=NULL, offset=NULL){
  zdat.qb = counts.associate(counts, condition, W=W, offset=offset)
  betahat = zdat.qb[3,]
  sebetahat = zdat.qb[4,]
  dispersion = zdat.qb[5,]
  df = length(condition)-2-!is.null(W)
  return(list(betahat = betahat, sebetahat = sebetahat,
         df = df, dispersion = dispersion))
}

# counts is a ngene (or nwindow) by nsample matrix of counts (eg RNAseq)
# g is an nsample vector of group memberships/covariate
# looks for association between rows and covariate
counts.associate=function(counts, g, W=NULL, offset=NULL, pseudocounts=1){
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

# Use RUV to estimate confounding factor
RUV_factor = function(counts, args, null){
  seq = newSeqExpressionSet(as.matrix(counts[as.logical(null),]))
  if (sum(null)>0){
    controls = rownames(seq)
    differences = matrix(data=c(1:args$Nsamp, (args$Nsamp+1):(2*args$Nsamp)), byrow=TRUE, nrow=2)
    seqRUVs = RUVs(seq, controls, k=args$RUV.k, differences)
    return(W=pData(seqRUVs)$W_1)
  }else{
    return(W=NULL)
  }
}

# Use SVA to estimate confounding factor
SVA_factor = function(counts, condition, args, null){
  mod1 = model.matrix(~condition)
  mod0 = cbind(mod1[,1])
  
  if (args$nullpi>0){
    svseq = svaseq(counts,mod1,mod0,control=null)
  }else{
    svseq = svaseq(counts,mod1,mod0)
  }
  
  if(svseq$n.sv>0){
    return(W = svseq$sv)
  }else{
    return(W = NULL)
  }
}

# Get sebetahat from edgeR.glm (infer from betahat & pval)
edgeR_glmest = function(counts, condition, args){
  design = model.matrix(~condition)
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
  
  return(list(betahat=betahat.edgeRglm, sebetahat=sebetahat.edgeRglm,
              df=df.edgeRglm)) 
}

# Get sebetahat from DESeq.glm (infer from betahat & pval)
DESeq_glmest = function(counts, condition, args){
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
  return(list(betahat=betahat.DESeqglm, sebetahat=sebetahat.DESeqglm,
              df=df.DESeqglm)) 
}