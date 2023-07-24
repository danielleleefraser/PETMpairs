

pairs_calc<-function(occ.matrix){
  library(cooccur)
  results<-list()
  cooccur.matrix<-cooccur(mat=occ.matrix,spp_names = TRUE)
  results[[1]]<-cooccur.matrix$results[cooccur.matrix$results$p_gt<=0.05,]
  results[[2]]<-cooccur.matrix$results[cooccur.matrix$results$p_lt<=0.05,]
  temp<-cooccur.matrix$results[cooccur.matrix$results$p_gt>0.05,]
  results[[3]]<-temp[temp$p_lt>0.05,]
  names(results)<-c("Aggregated pairs","Segregated pairs","Random pairs")
  return(results)
}


compare_climpref<-function(pairs_results,climprefs,calc.means=c("TRUE","FALSE")){
  pairs_clim<-matrix(ncol=12,nrow=1)
  diffs_pos_V1<-c()
  diffs_neg_V1<-c()
  diffs_rand_V1<-c()
  diffs_pos_V2<-c()
  diffs_neg_V2<-c()
  diffs_rand_V2<-c()
  positive_sig<-pairs_results[[1]]
  negative_sig<-pairs_results[[2]]
  random_sig<-pairs_results[[3]]
  for(i in 1:nrow(positive_sig)){
    focal<-positive_sig[i,]
    first_spec<-climprefs[focal$sp1_name,]
    second_spec<-climprefs[focal$sp2_name,]
    diffs_pos_V1[i]<-abs(first_spec$V1-second_spec$V1)
    diffs_pos_V2[i]<-abs(first_spec$V2-second_spec$V2)
  }
  for(i in 1:nrow(negative_sig)){
    focal<-negative_sig[i,]
    first_spec<-climprefs[focal$sp1_name,]
    second_spec<-climprefs[focal$sp2_name,]
    diffs_neg_V1[i]<-abs(first_spec$V1-second_spec$V1)
    diffs_neg_V2[i]<-abs(first_spec$V2-second_spec$V2)
  }
  for(i in 1:nrow(random_sig)){
    focal<-random_sig[i,]
    first_spec<-climprefs[focal$sp1_name,]
    second_spec<-climprefs[focal$sp2_name,]
    diffs_rand_V1[i]<-abs(first_spec$V1-second_spec$V1)
    diffs_rand_V2[i]<-abs(first_spec$V2-second_spec$V2)
  }
  if(calc.means=="TRUE"){ # add V1 and V2 component here
    
    pairs_clim[1,1]<-mean(na.omit(diffs_pos_V1))
    pairs_clim[1,2]<-mean(na.omit(diffs_neg_V1))
    pairs_clim[1,3]<-mean(na.omit(diffs_rand_V1))
    pairs_clim[1,4]<-sd(na.omit(diffs_pos_V1))/sqrt(length(diffs_pos_V1))
    pairs_clim[1,5]<-sd(na.omit(diffs_neg_V1))/sqrt(length(diffs_neg_V1))
    pairs_clim[1,6]<-sd(na.omit(diffs_rand_V1))/sqrt(length(diffs_rand_V1))
    
    pairs_clim[1,7]<-mean(na.omit(diffs_pos_V2))
    pairs_clim[1,8]<-mean(na.omit(diffs_neg_V2))
    pairs_clim[1,9]<-mean(na.omit(diffs_rand_V2))
    pairs_clim[1,10]<-sd(na.omit(diffs_pos_V2))/sqrt(length(diffs_pos_V2))
    pairs_clim[1,11]<-sd(na.omit(diffs_neg_V2))/sqrt(length(diffs_neg_V2))
    pairs_clim[1,12]<-sd(na.omit(diffs_rand_V2))/sqrt(length(diffs_rand_V2))
    
    return(pairs_clim)
  }
  if(calc.means=="FALSE"){
    results.list<-list(diffs_pos_V1,diffs_neg_V1,diffs_rand_V1,diffs_pos_V2,diffs_neg_V2,diffs_rand_V2)
    names(results.list)<-c("agg_v1","seg_v1","rand_v1","agg_v2","seg_v2","rand_v2")
    return(results.list)
  }
}


Summarize.null.clims<-function(null.list){
  
  clarkfork_three_aggV1<-c()
  clarkfork_three_segV1<-c()
  clarkfork_three_randV1<-c()
  
  Was_zero_aggV1<-c()
  Was_zero_segV1<-c()
  Was_zero_randV1<-c()
  
  Was_one_aggV1<-c()
  Was_one_segV1<-c()
  Was_one_randV1<-c()
  
  clarkfork_three_aggV2<-c()
  clarkfork_three_segV2<-c()
  clarkfork_three_randV2<-c()
  
  Was_zero_aggV2<-c()
  Was_zero_segV2<-c()
  Was_zero_randV2<-c()
  
  Was_one_aggV2<-c()
  Was_one_segV2<-c()
  Was_one_randV2<-c()
  
  results.final<-matrix(nrow=3,ncol=30)
  colnames(results.final)<-c("MeanNullaggV1","SENullaggV1","SDNullaggV1","UpperNullaggV1","LowerNullaggV1",
                             "MeanNullsegV1","SENullsegV1","SDNullsegV1","UpperNullsegV1","LowerNullsegV1",
                             "MeanNullrandV1","SENullrandV1","SDNullrandV1","UpperNullrandV1","LowerNullrandV1",
                             "MeanNullaggV2","SENullaggV2","SDNullaggV2","UpperNullaggV2","LowerNullaggV2",
                             "MeanNullsegV2","SENullsegV2","SDNullsegV2","UpperNullsegV2","LowerNullsegV2",
                             "MeanNullrandV2","SENullrandV2","SDNullrandV2","UpperNullrandV2","LowerNullrandV2")# double because it needs V1 and V2 etc.
  
  temp.list<-null.list[[1]] # each temp is an iteration
  
  clarkfork_three_aggV1<-temp.list[[1]][[1]]
  clarkfork_three_segV1<-temp.list[[1]][[2]]
  clarkfork_three_randV1<-temp.list[[1]][[3]]
  
  clarkfork_three_aggV2<-temp.list[[1]][[4]]
  clarkfork_three_segV2<-temp.list[[1]][[5]]
  clarkfork_three_randV2<-temp.list[[1]][[6]]
  
  Was_zero_aggV1<-temp.list[[2]][[1]]
  Was_zero_segV1<-temp.list[[2]][[2]]
  Was_zero_randV1<-temp.list[[2]][[3]]
  
  Was_zero_aggV2<-temp.list[[2]][[4]]
  Was_zero_segV2<-temp.list[[2]][[5]]
  Was_zero_randV2<-temp.list[[2]][[6]]
  
  Was_one_aggV1<-temp.list[[3]][[1]]
  Was_one_segV1<-temp.list[[3]][[2]]
  Was_one_randV1<-temp.list[[3]][[3]]
  
  Was_one_aggV2<-temp.list[[3]][[4]]
  Was_one_segV2<-temp.list[[3]][[5]]
  Was_one_randV2<-temp.list[[3]][[6]]
  
  for(i in 2:length(null.list)){
    
    temp.list<-null.list[[i]]
    
    clarkfork_three_aggV1 <- c(clarkfork_three_aggV1, temp.list[[1]][[1]])
    clarkfork_three_segV1 <- c(clarkfork_three_segV1, temp.list[[1]][[2]])
    clarkfork_three_randV1 <-c(clarkfork_three_randV1, temp.list[[1]][[3]])
    
    clarkfork_three_aggV2 <- c(clarkfork_three_aggV2, temp.list[[1]][[4]])
    clarkfork_three_segV2 <- c(clarkfork_three_segV2, temp.list[[1]][[5]])
    clarkfork_three_randV2 <-c(clarkfork_three_randV2, temp.list[[1]][[6]])
    
    Was_zero_aggV1 <- c(Was_zero_aggV1, temp.list[[2]][[1]])
    Was_zero_segV1 <- c(Was_zero_segV1, temp.list[[2]][[2]])
    Was_zero_randV1 <-c(Was_zero_randV1, temp.list[[2]][[3]])
    
    Was_zero_aggV2 <- c(Was_zero_aggV2, temp.list[[2]][[4]])
    Was_zero_segV2 <- c(Was_zero_segV2, temp.list[[2]][[5]])
    Was_zero_randV2 <-c(Was_zero_randV2, temp.list[[2]][[6]])
    
    Was_one_aggV1 <- c(Was_one_aggV1, temp.list[[3]][[1]])
    Was_one_segV1 <- c(Was_one_segV1, temp.list[[3]][[2]])
    Was_one_randV1 <-c(Was_one_randV1, temp.list[[3]][[3]])
    
    Was_one_aggV2 <- c(Was_one_aggV2, temp.list[[3]][[4]])
    Was_one_segV2 <- c(Was_one_segV2, temp.list[[3]][[5]])
    Was_one_randV2 <-c(Was_one_randV2, temp.list[[3]][[6]])
    
  }
  #Clarkfork3 V1
  results.final[1,1]<-mean(clarkfork_three_aggV1,na.rm = T)
  results.final[1,2]<-sd(clarkfork_three_aggV1,na.rm = T)/sqrt(length(clarkfork_three_aggV1)) #
  results.final[1,3]<-sd(clarkfork_three_aggV1,na.rm = T)
  results.final[1,4]<-(quantile(clarkfork_three_aggV1, c(0.025), type = 1,na.rm=TRUE))
  results.final[1,5]<-(quantile(clarkfork_three_aggV1, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[1,6]<-mean(clarkfork_three_segV1,na.rm = T)
  results.final[1,7]<-sd(clarkfork_three_segV1,na.rm = T)/sqrt(length(clarkfork_three_segV1)) #
  results.final[1,8]<-sd(clarkfork_three_segV1,na.rm = T)
  results.final[1,9]<-(quantile(clarkfork_three_segV1, c(0.025), type = 1,na.rm=TRUE))
  results.final[1,10]<-(quantile(clarkfork_three_segV1, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[1,11]<-mean(clarkfork_three_randV1,na.rm = T)
  results.final[1,12]<-sd(clarkfork_three_randV1,na.rm = T)/sqrt(length(clarkfork_three_randV1)) #
  results.final[1,13]<-sd(clarkfork_three_randV1,na.rm = T)
  results.final[1,14]<-(quantile(clarkfork_three_randV1, c(0.025), type = 1,na.rm=TRUE))
  results.final[1,15]<-(quantile(clarkfork_three_randV1, c(0.975), type = 1,na.rm=TRUE))
  
  #Clarkfork3 V2
  results.final[1,16]<-mean(clarkfork_three_aggV2,na.rm = T)
  results.final[1,17]<-sd(clarkfork_three_aggV2,na.rm = T)/sqrt(length(clarkfork_three_aggV2)) #
  results.final[1,18]<-sd(clarkfork_three_aggV2,na.rm = T)
  results.final[1,19]<-(quantile(clarkfork_three_aggV2, c(0.025), type = 1,na.rm=TRUE))
  results.final[1,20]<-(quantile(clarkfork_three_aggV2, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[1,21]<-mean(clarkfork_three_segV2,na.rm = T)
  results.final[1,22]<-sd(clarkfork_three_segV2,na.rm = T)/sqrt(length(clarkfork_three_segV2)) #
  results.final[1,23]<-sd(clarkfork_three_segV2,na.rm = T)
  results.final[1,24]<-(quantile(clarkfork_three_segV2, c(0.025), type = 1,na.rm=TRUE))
  results.final[1,25]<-(quantile(clarkfork_three_segV2, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[1,26]<-mean(clarkfork_three_randV2,na.rm = T)
  results.final[1,27]<-sd(clarkfork_three_randV2,na.rm = T)/sqrt(length(clarkfork_three_randV2)) #
  results.final[1,28]<-sd(clarkfork_three_randV2,na.rm = T)
  results.final[1,29]<-(quantile(clarkfork_three_randV2, c(0.025), type = 1,na.rm=TRUE))
  results.final[1,30]<-(quantile(clarkfork_three_randV2, c(0.975), type = 1,na.rm=TRUE))
  
  #Wasatch0 V1 Was_zero_
  results.final[2,1]<-mean(Was_zero_aggV1,na.rm = T)
  results.final[2,2]<-sd(Was_zero_aggV1,na.rm = T)/sqrt(length(Was_zero_aggV1)) #
  results.final[2,3]<-sd(Was_zero_aggV1,na.rm = T)
  results.final[2,4]<-(quantile(Was_zero_aggV1, c(0.025), type = 1,na.rm=TRUE))
  results.final[2,5]<-(quantile(Was_zero_aggV1, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[2,6]<-mean(Was_zero_segV1,na.rm = T)
  results.final[2,7]<-sd(Was_zero_segV1,na.rm = T)/sqrt(length(Was_zero_segV1)) #
  results.final[2,8]<-sd(Was_zero_segV1,na.rm = T)
  results.final[2,9]<-(quantile(Was_zero_segV1, c(0.025), type = 1,na.rm=TRUE))
  results.final[2,10]<-(quantile(Was_zero_segV1, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[2,11]<-mean(Was_zero_randV1,na.rm = T)
  results.final[2,12]<-sd(Was_zero_randV1,na.rm = T)/sqrt(length(Was_zero_randV1)) #
  results.final[2,13]<-sd(Was_zero_randV1,na.rm = T)
  results.final[2,14]<-(quantile(Was_zero_randV1, c(0.025), type = 1,na.rm=TRUE))
  results.final[2,15]<-(quantile(Was_zero_randV1, c(0.975), type = 1,na.rm=TRUE))
  
  #Wasatch0 V2 Was_zero_
  results.final[2,16]<-mean(Was_zero_aggV2,na.rm = T)
  results.final[2,17]<-sd(Was_zero_aggV2,na.rm = T)/sqrt(length(Was_zero_aggV2)) #
  results.final[2,18]<-sd(Was_zero_aggV2,na.rm = T)
  results.final[2,19]<-(quantile(Was_zero_aggV2, c(0.025), type = 1,na.rm=TRUE))
  results.final[2,20]<-(quantile(Was_zero_aggV2, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[2,21]<-mean(Was_zero_segV2,na.rm = T)
  results.final[2,22]<-sd(Was_zero_segV2,na.rm = T)/sqrt(length(Was_zero_segV2)) #
  results.final[2,23]<-sd(Was_zero_segV2,na.rm = T)
  results.final[2,24]<-(quantile(Was_zero_segV2, c(0.025), type = 1,na.rm=TRUE))
  results.final[2,25]<-(quantile(Was_zero_segV2, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[2,26]<-mean(Was_zero_randV2,na.rm = T)
  results.final[2,27]<-sd(Was_zero_randV2,na.rm = T)/sqrt(length(Was_zero_randV2)) #
  results.final[2,28]<-sd(Was_zero_randV2,na.rm = T)
  results.final[2,29]<-(quantile(Was_zero_randV2, c(0.025), type = 1,na.rm=TRUE))
  results.final[2,30]<-(quantile(Was_zero_randV2, c(0.975), type = 1,na.rm=TRUE))
  
  #Wasatch one Was_one_
  results.final[3,1]<-mean(Was_one_aggV1,na.rm = T)
  results.final[3,2]<-sd(Was_one_aggV1,na.rm = T)/sqrt(length(Was_one_aggV1)) #
  results.final[3,3]<-sd(Was_one_aggV1,na.rm = T)
  results.final[3,4]<-(quantile(Was_one_aggV1, c(0.025), type = 1,na.rm=TRUE))
  results.final[3,5]<-(quantile(Was_one_aggV1, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[3,6]<-mean(Was_one_segV1,na.rm = T)
  results.final[3,7]<-sd(Was_one_segV1,na.rm = T)/sqrt(length(Was_one_segV1)) #
  results.final[3,8]<-sd(Was_one_segV1,na.rm = T)
  results.final[3,9]<-(quantile(Was_one_segV1, c(0.025), type = 1,na.rm=TRUE))
  results.final[3,10]<-(quantile(Was_one_segV1, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[3,11]<-mean(Was_one_randV1,na.rm = T)
  results.final[3,12]<-sd(Was_one_randV1,na.rm = T)/sqrt(length(Was_one_randV1)) #
  results.final[3,13]<-sd(Was_one_randV1,na.rm = T)
  results.final[3,14]<-(quantile(Was_one_randV1, c(0.025), type = 1,na.rm=TRUE))
  results.final[3,15]<-(quantile(Was_one_randV1, c(0.975), type = 1,na.rm=TRUE))
  
  #Wasatch0 V2 Was_zero_
  results.final[3,16]<-mean(Was_one_aggV2,na.rm = T)
  results.final[3,17]<-sd(Was_one_aggV2,na.rm = T)/sqrt(length(Was_one_aggV2)) #
  results.final[3,18]<-sd(Was_one_aggV2,na.rm = T)
  results.final[3,19]<-(quantile(Was_one_aggV2, c(0.025), type = 1,na.rm=TRUE))
  results.final[3,20]<-(quantile(Was_one_aggV2, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[3,21]<-mean(Was_one_segV2,na.rm = T)
  results.final[3,22]<-sd(Was_one_segV2,na.rm = T)/sqrt(length(Was_one_segV2)) #
  results.final[3,23]<-sd(Was_one_segV2,na.rm = T)
  results.final[3,24]<-(quantile(Was_one_segV2, c(0.025), type = 1,na.rm=TRUE))
  results.final[3,25]<-(quantile(Was_one_segV2, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[3,26]<-mean(Was_one_randV2,na.rm = T)
  results.final[3,27]<-sd(Was_one_randV2,na.rm = T)/sqrt(length(Was_one_randV2)) #
  results.final[3,28]<-sd(Was_one_randV2,na.rm = T)
  results.final[3,29]<-(quantile(Was_one_randV2, c(0.025), type = 1,na.rm=TRUE))
  results.final[3,30]<-(quantile(Was_one_randV2, c(0.975), type = 1,na.rm=TRUE))
  
  return(results.final)
}
  

compare_mass<-function(pairs_results,BMs,calc.means=c("TRUE","FALSE")){
  pairs_BM<-matrix(ncol=6,nrow=1)
  diffs_pos<-c()
  diffs_neg<-c()
  diffs_rand<-c()
  positive_sig<-pairs_results[[1]]
  negative_sig<-pairs_results[[2]]
  random_sig<-pairs_results[[3]]
  for(i in 1:nrow(positive_sig)){
    focal<-positive_sig[i,]
    first_spec<-BMs[focal$sp1_name,]
    second_spec<-BMs[focal$sp2_name,]
    diffs_pos[i]<-abs(first_spec$ln_mass-second_spec$ln_mass)
  }
  for(i in 1:nrow(negative_sig)){
    focal<-negative_sig[i,]
    first_spec<-BMs[focal$sp1_name,]
    second_spec<-BMs[focal$sp2_name,]
    diffs_neg[i]<-abs(first_spec$ln_mass-second_spec$ln_mass)
  }
  for(i in 1:nrow(random_sig)){
    focal<-random_sig[i,]
    first_spec<-BMs[focal$sp1_name,]
    second_spec<-BMs[focal$sp2_name,]
    diffs_rand[i]<-abs(first_spec$ln_mass-second_spec$ln_mass)
  }
  if(calc.means=="TRUE"){
    pairs_BM[1,1]<-mean(na.omit(diffs_pos))
    pairs_BM[1,2]<-mean(na.omit(diffs_neg))
    pairs_BM[1,3]<-mean(na.omit(diffs_rand))
    pairs_BM[1,4]<-sd(na.omit(diffs_pos))/sqrt(length(diffs_pos))
    pairs_BM[1,5]<-sd(na.omit(diffs_neg))/sqrt(length(diffs_neg))
    pairs_BM[1,6]<-sd(na.omit(diffs_rand))/sqrt(length(diffs_rand))
    return(pairs_BM)
  }
  if(calc.means=="FALSE"){
    results.list<-list(diffs_pos,diffs_neg,diffs_rand)
    names(results.list)<-c("agg","seg","rand")
    return(results.list)
  }
}


# compare locomotion

compare_loco<-function(pairs_results,loco,calc.means=c("TRUE","FALSE")){
  pairs_loco<-matrix(ncol=6,nrow=1)
  diffs_pos<-c()
  diffs_neg<-c()
  diffs_rand<-c()
  positive_sig<-pairs_results[[1]]
  negative_sig<-pairs_results[[2]]
  random_sig<-pairs_results[[3]]
  for(i in 1:nrow(positive_sig)){
    focal<-positive_sig[i,]
    first_spec<-loco[focal$sp1_name,]
    second_spec<-loco[focal$sp2_name,]
    diffs_pos[i]<-abs(first_spec$npose-second_spec$npose)
  }
  for(i in 1:nrow(negative_sig)){
    focal<-negative_sig[i,]
    first_spec<-loco[focal$sp1_name,]
    second_spec<-loco[focal$sp2_name,]
    diffs_neg[i]<-abs(first_spec$npose-second_spec$npose)
  }
  for(i in 1:nrow(random_sig)){
    focal<-random_sig[i,]
    first_spec<-loco[focal$sp1_name,]
    second_spec<-loco[focal$sp2_name,]
    diffs_rand[i]<-abs(first_spec$npose-second_spec$npose)
  }
  if(calc.means=="TRUE"){
    pairs_loco[1,1]<-mean(na.omit(diffs_pos))
    pairs_loco[1,2]<-mean(na.omit(diffs_neg))
    pairs_loco[1,3]<-mean(na.omit(diffs_rand))
    pairs_loco[1,4]<-sd(na.omit(diffs_pos))/sqrt(length(diffs_pos))
    pairs_loco[1,5]<-sd(na.omit(diffs_neg))/sqrt(length(diffs_neg))
    pairs_loco[1,6]<-sd(na.omit(diffs_rand))/sqrt(length(diffs_rand))
    return(pairs_loco)
  }
  if(calc.means=="FALSE"){
    results.list<-list(diffs_pos,diffs_neg,diffs_rand)
    names(results.list)<-c("agg","seg","rand")
    return(results.list)
  }
}


Summarize.null.morph<-function(null.list){ # summarizes just for the three time periods but not across all the iterations
  results.final<-matrix(nrow=3,ncol=15)
  
  clarkfork_three_agg<-c()
  clarkfork_three_seg<-c()
  clarkfork_three_rand<-c()
  
  Was_zero_agg<-c()
  Was_zero_seg<-c()
  Was_zero_rand<-c()
  
  Was_one_agg<-c()
  Was_one_seg<-c()
  Was_one_rand<-c()
  
  colnames(results.final)<-c("MeanNullagg","SENullagg","SSNullagg","LowerNullagg","UpperNullagg","MeanNullseg","SENullseg","SDNullseg","LowerNullseg","UpperNullseg",
                             "MeanNullrand","SENullrand","SDNullrand","LowerNullrand","UpperNullrand")
  temp.list<-null.list[[1]] 
  
  clarkfork_three_agg<-temp.list[[1]][[1]]
  clarkfork_three_seg<-temp.list[[1]][[2]]
  clarkfork_three_rand<-temp.list[[1]][[3]]
  
  Was_zero_agg<-temp.list[[2]][[1]]
  Was_zero_seg<-temp.list[[2]][[2]]
  Was_zero_rand<-temp.list[[2]][[3]]
  
  Was_one_agg<-temp.list[[3]][[1]]
  Was_one_seg<-temp.list[[3]][[2]]
  Was_one_rand<-temp.list[[3]][[3]]
  
  for(i in 2:length(null.list)) {
    
    temp.list <- null.list[[i]]
    
    clarkfork_three_agg <- c(clarkfork_three_agg, temp.list[[1]][[1]])
    clarkfork_three_seg <- c(clarkfork_three_seg, temp.list[[1]][[2]])
    clarkfork_three_rand <-c(clarkfork_three_rand, temp.list[[1]][[3]])
    
    Was_zero_agg <- c(Was_zero_agg, temp.list[[2]][[1]])
    Was_zero_seg <- c(Was_zero_seg, temp.list[[2]][[2]])
    Was_zero_rand <-c(Was_zero_rand, temp.list[[2]][[3]])
    
    Was_one_agg <- c(Was_one_agg, temp.list[[3]][[1]])
    Was_one_seg <- c(Was_one_seg, temp.list[[3]][[2]])
    Was_one_rand <-c(Was_one_rand, temp.list[[3]][[3]])
    
  }
  #Clarkfork3
  results.final[1,1]<-mean(clarkfork_three_agg,na.rm = T)
  results.final[1,2]<-sd(clarkfork_three_agg,na.rm = T)/sqrt(length(clarkfork_three_agg)) #
  results.final[1,3]<-sd(clarkfork_three_agg,na.rm = T)
  results.final[1,4]<-(quantile(clarkfork_three_agg, c(0.025), type = 1,na.rm=TRUE))
  results.final[1,5]<-(quantile(clarkfork_three_agg, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[1,6]<-mean(clarkfork_three_seg,na.rm = T)
  results.final[1,7]<-sd(clarkfork_three_seg,na.rm = T)/sqrt(length(clarkfork_three_seg)) #
  results.final[1,8]<-sd(clarkfork_three_seg,na.rm = T)
  results.final[1,9]<-(quantile(clarkfork_three_seg, c(0.025), type = 1,na.rm=TRUE))
  results.final[1,10]<-(quantile(clarkfork_three_seg, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[1,11]<-mean(clarkfork_three_rand,na.rm = T)
  results.final[1,12]<-sd(clarkfork_three_rand,na.rm = T)/sqrt(length(clarkfork_three_rand)) #
  results.final[1,13]<-sd(clarkfork_three_rand,na.rm = T)
  results.final[1,14]<-(quantile(clarkfork_three_rand, c(0.025), type = 1,na.rm=TRUE))
  results.final[1,15]<-(quantile(clarkfork_three_rand, c(0.975), type = 1,na.rm=TRUE))
  
  #Wasatch 0
  results.final[2,1]<-mean(Was_zero_agg,na.rm = T)
  results.final[2,2]<-sd(Was_zero_agg,na.rm = T)/sqrt(length(Was_zero_agg)) #
  results.final[2,3]<-sd(Was_zero_agg,na.rm = T)
  results.final[2,4]<-(quantile(Was_zero_agg, c(0.025), type = 1,na.rm=TRUE))
  results.final[2,5]<-(quantile(Was_zero_agg, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[2,6]<-mean(Was_zero_seg,na.rm = T)
  results.final[2,7]<-sd(Was_zero_seg,na.rm = T)/sqrt(length(Was_zero_seg)) #
  results.final[2,8]<-sd(Was_zero_seg,na.rm = T)
  results.final[2,9]<-(quantile(Was_zero_seg, c(0.025), type = 1,na.rm=TRUE))
  results.final[2,10]<-(quantile(Was_zero_seg, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[2,11]<-mean(Was_zero_rand,na.rm = T)
  results.final[2,12]<-sd(Was_zero_rand,na.rm = T)/sqrt(length(Was_zero_rand)) #
  results.final[2,13]<-sd(Was_zero_rand,na.rm = T)
  results.final[2,14]<-(quantile(Was_zero_rand, c(0.025), type = 1,na.rm=TRUE))
  results.final[2,15]<-(quantile(Was_zero_rand, c(0.975), type = 1,na.rm=TRUE))
  
  #Wasatch one
  results.final[3,1]<-mean(Was_one_agg,na.rm = T)
  results.final[3,2]<-sd(Was_one_agg,na.rm = T)/sqrt(length(Was_one_agg)) #
  results.final[3,3]<-sd(Was_one_agg,na.rm = T)
  results.final[3,4]<-(quantile(Was_one_agg, c(0.025), type = 1,na.rm=TRUE))
  results.final[3,5]<-(quantile(Was_one_agg, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[3,6]<-mean(Was_one_seg,na.rm = T)
  results.final[3,7]<-sd(Was_one_seg,na.rm = T)/sqrt(length(Was_one_seg)) #
  results.final[3,8]<-sd(Was_one_seg,na.rm = T)
  results.final[3,9]<-(quantile(Was_one_seg, c(0.025), type = 1,na.rm=TRUE))
  results.final[3,10]<-(quantile(Was_one_seg, c(0.975), type = 1,na.rm=TRUE))
  
  results.final[3,11]<-mean(Was_one_rand,na.rm = T)
  results.final[3,12]<-sd(Was_one_rand,na.rm = T)/sqrt(length(Was_one_rand)) #
  results.final[3,13]<-sd(Was_one_rand,na.rm = T)
  results.final[3,14]<-(quantile(Was_one_rand, c(0.025), type = 1,na.rm=TRUE))
  results.final[3,15]<-(quantile(Was_one_rand, c(0.975), type = 1,na.rm=TRUE))
  return(results.final)
}



####