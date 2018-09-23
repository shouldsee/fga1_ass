calc_BatchGroupCellMeans <-function(sim.groups,expr_mat = 'counts')
{
  # sim.groups<-sim.groups.0
  cD<-colData(sim.groups)
  # ct <- counts(sim.groups)
  ct <- assays(sim.groups)[[expr_mat]]
  # ASS$
  # i = 1
  # sapply(1:nlevels(gp_fac),function(i)
  # data.frame(sim.groups)
  gp_vct  = (cD$Group)
  bch_vct = (cD$Batch)
  gp_fac  = factor(cD$Group)
  bch_fac = factor(cD$Batch)
  bch_gp <- cbind(bch_vct,gp_vct)
  lvl = data.frame(unique(bch_gp))
  lvl$id <- 1:nrow(lvl)
  # lvl
  # levels(lvl[,1]) <-levels(bch_fac)
  # levels(lvl[,2]) <-levels(gp_fac)
  # as.factor(lvl)
  # bch_gp
  # levels(bch_gp[,1])
  # sapply(1:ncol(lvl),function(i){levels(bch_gp[,i])[lvl[,i]] })
  # lvl[,1] <-levels
  # levels(bch_gp)
  ROW=lvl[1,]
  # lvl
  # ROW[2]
  # ROW
  
  GroupMeans<- apply( lvl,MARGIN= 1,function(ROW)
  {
    # fac <-levels(gp_fac)[i]
    bch = ROW[1]
    gp = ROW[2]
    idx =((gp_fac== gp) & (bch_fac==bch))
    print(sum(idx))
    rowSums(ct[,idx])%>%
      as.integer
    
    # rowMedians(ct[,idx])%>%as.integer
  }
  )
  
  # SingleCellExperiment(GroupMeans)
  # ?SingleCellExperiment
  # GroupMeans%>%head
  # bch_gp%>%distinct
  # bch_gp == lvl
  
  # lvl
  # idx
  idxdf<-as.data.frame(bch_gp) %>% left_join(lvl)
  # idxdf<-as.data.frame(bch_gp) %>% left_join(lvl) %>%select('id')%>%unlist
  idx<-idxdf[['id']]
  BatchGroupCellMeans=GroupMeans[,idx]
  colnames(BatchGroupCellMeans) <- cD$Cell 
  # ASS <- assays(sim.groups)
  # ASS$BatchGroupCellMeans = BatchGroupCellMeans
  # assays(sim.groups) <- ASS
  assays(sim.groups)$BatchGroupCellMeans <- BatchGroupCellMeans
  assays(sim.groups)$logBatchGroupCellMeans <- log(assays(sim.groups)$BatchGroupCellMeans + 1)
  # sim.groups.0 <- sim.groups
  return(sim.groups)
}


#### RORC wrappers (ROC evaluation)

DE_Quality_AUC <- function(pVals) {
  pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                   names(pVals) %in% GroundTruth$notDE]
  truth <- rep(1, times = length(pVals));
  truth[names(pVals) %in% GroundTruth$DE] = 0;
  pred <- ROCR::prediction(pVals, truth)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  # ROCR::plot(perf)
  (aucObj <- ROCR::performance(pred, "auc"))
  # print(aucObj@y.values[[1]])
  return(aucObj@y.values[[1]])
  # return(perf)
}

DE_Quality_ROC <- function(pVals) {
  pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                   names(pVals) %in% GroundTruth$notDE]
  truth <- rep(1, times = length(pVals));
  truth[names(pVals) %in% GroundTruth$DE] = 0;
  pred <- ROCR::prediction(pVals, truth)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  # ROCR::plot(perf)
  # (aucObj <- ROCR::performance(pred, "auc"))
  # print(aucObj@y.values[[1]])
  return(perf)
}

DE_Quality_PR <- function(pVals) {
  pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                   names(pVals) %in% GroundTruth$notDE]
  truth <- rep(1, times = length(pVals));
  truth[names(pVals) %in% GroundTruth$DE] = 0;
  pred <- ROCR::prediction(pVals, truth)
  perf <- ROCR::performance(pred,  "prec","rec")
  # perf <- ROCR::performance(pred, "tpr", "prec")
  # ROCR::plot(perf)
  (prObj <- ROCR::performance(pred, "prbe"))
  print(prObj@y.values[[1]])
  # return(aucObj@y.values[[1]])
  return(perf)
}
DE_Quality_PRBE <- function(pVals) {
  pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                   names(pVals) %in% GroundTruth$notDE]
  truth <- rep(1, times = length(pVals));
  truth[names(pVals) %in% GroundTruth$DE] = 0;
  pred <- ROCR::prediction(pVals, truth)
  perf <- ROCR::performance(pred,  "prec","rec")
  # perf <- ROCR::performance(pred, "tpr", "prec")
  # ROCR::plot(perf)
  (prObj <- ROCR::performance(pred, "prbe"))
  print(prObj@y.values)
  return(tail(prObj@y.values[[1]],1))
  
}