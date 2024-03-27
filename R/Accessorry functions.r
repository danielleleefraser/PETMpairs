

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
    
    colnames(pairs_clim)<-c("Mean V1 agg","Mean V1 seg","Mean V1 rand",
                            "SD V1 agg","SD V1 seg","SD V1 rand",
                            "Mean V2 agg","Mean V2 seg","Mean V2 rand",
                            "SD V2 agg","SD V2 seg","SD V2 rand")
    
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
    colnames(pairs_BM)<-c("Mean agg","Mean seg","Mean rand","SD agg","SD seg","SD rand")
    return(pairs_BM)
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
  
  colnames(results.final)<-c("MeanNullagg","SENullagg","SDNullagg","LowerNullagg","UpperNullagg","MeanNullseg","SENullseg","SDNullseg","LowerNullseg","UpperNullseg",
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


Shuffletaxa<-function(occur,bins,clim_prefs,BMs,loco,niter,summarize=c("TRUE","FALSE")){
  library(EcoSimR)#must download archived version
  clim_list<-list()
  BM_list<-list()
  loco_list<-list()
  for(f in 1:niter){
    clim_comp<-list()
    BM_comp<-list()
    loco_comp<-list()
    randommatrix<-sim3(occur)
    rownames(randommatrix)<-rownames(occur)
    for(b in 1:length(bins)){
      matches_sites<-intersect(bins[[b]],rownames(randommatrix))
      matrix1<-randommatrix[matches_sites,]
      matrix1<-t(matrix1)
      sig_pairs<-pairs_calc(matrix1)
      clim_comp[[b]]<-compare_climpref(sig_pairs,clim_prefs,calc.means = FALSE)
      BM_comp[[b]]<-compare_mass(sig_pairs,BMs,calc.means = FALSE) # needs to be lists
      loco_comp[[b]]<-compare_loco(sig_pairs,loco,calc.means = FALSE)
    }
    clim_list[[f]]<-clim_comp #1
    BM_list[[f]]<-BM_comp #2
    loco_list[[f]]<-loco_comp #3
  }
  if(summarize=="TRUE"){
    BM_null<-Summarize.null.morph(BM_list)# has to work with a list.
    Loco_null<-Summarize.null.morph(loco_list)
    Clim_null<-Summarize.null.clims(clim_list)
    all_null<-list(BM_null,Loco_null,Clim_null)
  }
  if(summarize=="FALSE"){
    BM_null<-Bind_null_morph(BM_list)
    Clim_null<-Bind_null_clim(clim_list) # have to edit climt summary function
    all_null<-list(BM_null,Clim_null)
  }
  return(all_null)
}



Bind_null_morph<-function(null.list){
  
  clarkfork_three_agg<-c()#1 (BM)
  clarkfork_three_seg<-c()
  clarkfork_three_rand<-c()
  
  Was_zero_agg<-c() #BM
  Was_zero_seg<-c()
  Was_zero_rand<-c()
  
  Was_one_agg<-c() #BM
  Was_one_seg<-c()
  Was_one_rand<-c()
  
  temp.list<-null.list[[1]] 
  
  clarkfork_three_agg<-c(temp.list[[1]][[1]])
  clarkfork_three_seg<-c(temp.list[[1]][[2]])
  clarkfork_three_rand<-c(temp.list[[1]][[3]])
  
  Was_zero_agg<-c(temp.list[[2]][[1]])
  Was_zero_seg<-c(temp.list[[2]][[2]])
  Was_zero_rand<-c(temp.list[[2]][[3]])
  
  Was_one_agg<-c(temp.list[[3]][[1]])
  Was_one_seg<-c(temp.list[[3]][[2]])
  Was_one_rand<-c(temp.list[[3]][[3]])
  
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
  
  clarkfork_three_agg<-data.frame(clarkfork_three_agg)
  colnames(clarkfork_three_agg)<-c("n")
  clarkfork_three_rand<-data.frame(clarkfork_three_rand)
  colnames(clarkfork_three_rand)<-c("n")
  clarkfork_three_seg<-data.frame(clarkfork_three_seg)
  colnames(clarkfork_three_seg)<-c("n")
  
  clarkforknames<-c(rep("agg",nrow(clarkfork_three_agg)),rep("rand",nrow(clarkfork_three_rand)),
                    rep("seg",nrow(clarkfork_three_seg)))
  clarkfork<-data.frame(rbind(clarkfork_three_agg,clarkfork_three_rand,clarkfork_three_seg))
  clarkfork<-cbind(clarkforknames,clarkfork)
  
  Was_zero_agg<-data.frame(Was_zero_agg)
  colnames(Was_zero_agg)<-c("n")
  Was_zero_rand<-data.frame(Was_zero_rand)
  colnames(Was_zero_rand)<-c("n")
  Was_zero_seg<-data.frame(Was_zero_seg)
  colnames(Was_zero_seg)<-c("n")
  
  Was_zeronames<-c(rep("agg",nrow(Was_zero_agg)),rep("rand",nrow(Was_zero_rand)),
                   rep("seg",nrow(Was_zero_seg)))
  Was_zero<-data.frame(rbind(Was_zero_agg,Was_zero_rand,Was_zero_seg))
  Was_zero<-cbind(Was_zeronames,Was_zero)
  
  Was_one_agg<-data.frame(Was_one_agg)
  colnames(Was_one_agg)<-c("n")
  Was_one_rand<-data.frame(Was_one_rand)
  colnames(Was_one_rand)<-c("n")
  Was_one_seg<-data.frame(Was_one_seg)
  colnames(Was_one_seg)<-c("n")
  
  Was_onenames<-c(rep("agg",nrow(Was_one_agg)),rep("rand",nrow(Was_one_rand)),
                  rep("seg",nrow(Was_one_seg)))
  Was_one<-data.frame(rbind(Was_one_agg,Was_one_rand,Was_one_seg))
  Was_one<-cbind(Was_onenames,Was_one)
  
  results.final<-list(clarkfork,Was_one,Was_zero)
  return(results.final)
}


Bind_null_clim<-function(null.list){
  
  clarkfork_three_agg_v1<-c()#1 (BM)
  clarkfork_three_seg_v1<-c()
  clarkfork_three_rand_v1<-c()
  
  Was_zero_agg_v1<-c() #BM
  Was_zero_seg_v1<-c()
  Was_zero_rand_v1<-c()
  
  Was_one_agg_v1<-c() #BM
  Was_one_seg_v1<-c()
  Was_one_rand_v1<-c()
  
  clarkfork_three_agg_v2<-c()#1 (BM)
  clarkfork_three_seg_v2<-c()
  clarkfork_three_rand_v2<-c()
  
  Was_zero_agg_v2<-c() #BM
  Was_zero_seg_v2<-c()
  Was_zero_rand_v2<-c()
  
  Was_one_agg_v2<-c() #BM
  Was_one_seg_v2<-c()
  Was_one_rand_v2<-c()
  
  temp.list.v1<-null.list[[1]] 
  
  clarkfork_three_agg_v1<-c(temp.list.v1[[1]][[1]])
  clarkfork_three_seg_v1<-c(temp.list.v1[[1]][[2]])
  clarkfork_three_rand_v1<-c(temp.list.v1[[1]][[3]])
  
  Was_zero_agg_v1<-c(temp.list.v1[[2]][[1]])
  Was_zero_seg_v1<-c(temp.list.v1[[2]][[2]])
  Was_zero_rand_v1<-c(temp.list.v1[[2]][[3]])
  
  Was_one_agg_v1<-c(temp.list.v1[[3]][[1]])
  Was_one_seg_v1<-c(temp.list.v1[[3]][[2]])
  Was_one_rand_v1<-c(temp.list.v1[[3]][[3]])
  
  temp.list.v2<-null.list[[1]] 
  
  clarkfork_three_agg_v2<-c(temp.list.v2[[1]][[1]])
  clarkfork_three_seg_v2<-c(temp.list.v2[[1]][[2]])
  clarkfork_three_rand_v2<-c(temp.list.v2[[1]][[3]])
  
  Was_zero_agg_v2<-c(temp.list.v2[[2]][[1]])
  Was_zero_seg_v2<-c(temp.list.v2[[2]][[2]])
  Was_zero_rand_v2<-c(temp.list.v2[[2]][[3]])
  
  Was_one_agg_v2<-c(temp.list.v2[[3]][[1]])
  Was_one_seg_v2<-c(temp.list.v2[[3]][[2]])
  Was_one_rand_v2<-c(temp.list.v2[[3]][[3]])
  
  for(i in 2:length(null.list)) {
    
    temp.list_v1 <- null.list[[i]]
    
    clarkfork_three_agg_v1 <- c(clarkfork_three_agg_v1, temp.list_v1[[1]][[1]])
    clarkfork_three_seg_v1 <- c(clarkfork_three_seg_v1, temp.list_v1[[1]][[2]])
    clarkfork_three_rand_v1 <-c(clarkfork_three_rand_v1, temp.list_v1[[1]][[3]])
    
    Was_zero_agg_v1 <- c(Was_zero_agg_v1, temp.list_v1[[2]][[1]])
    Was_zero_seg_v1 <- c(Was_zero_seg_v1, temp.list_v1[[2]][[2]])
    Was_zero_rand_v1 <-c(Was_zero_rand_v1, temp.list_v1[[2]][[3]])
    
    Was_one_agg_v1 <- c(Was_one_agg_v1, temp.list_v1[[3]][[1]])
    Was_one_seg_v1 <- c(Was_one_seg_v1, temp.list_v1[[3]][[2]])
    Was_one_rand_v1 <-c(Was_one_rand_v1, temp.list_v1[[3]][[3]])
    
    temp.list_v2 <- null.list[[i]]
    
    clarkfork_three_agg_v2 <- c(clarkfork_three_agg_v2, temp.list_v2[[1]][[1]])
    clarkfork_three_seg_v2 <- c(clarkfork_three_seg_v2, temp.list_v2[[1]][[2]])
    clarkfork_three_rand_v2 <-c(clarkfork_three_rand_v2, temp.list_v2[[1]][[3]])
    
    Was_zero_agg_v2 <- c(Was_zero_agg_v2, temp.list_v2[[2]][[1]])
    Was_zero_seg_v2 <- c(Was_zero_seg_v2, temp.list_v2[[2]][[2]])
    Was_zero_rand_v2 <-c(Was_zero_rand_v2, temp.list_v2[[2]][[3]])
    
    Was_one_agg_v2 <- c(Was_one_agg_v2, temp.list_v2[[3]][[1]])
    Was_one_seg_v2 <- c(Was_one_seg_v2, temp.list_v2[[3]][[2]])
    Was_one_rand_v2 <-c(Was_one_rand_v2, temp.list_v2[[3]][[3]])
    
  }
  
  clarkfork_three_agg_v1<-data.frame(clarkfork_three_agg_v1)
  colnames(clarkfork_three_agg_v1)<-c("n")
  clarkfork_three_rand_v1<-data.frame(clarkfork_three_rand_v1)
  colnames(clarkfork_three_rand_v1)<-c("n")
  clarkfork_three_seg_v1<-data.frame(clarkfork_three_seg_v1)
  colnames(clarkfork_three_seg_v1)<-c("n")
  
  clarkforknames_v1<-c(rep("agg",nrow(clarkfork_three_agg_v1)),rep("rand",nrow(clarkfork_three_rand_v1)),
                       rep("seg",nrow(clarkfork_three_seg_v1)))
  clarkfork_v1<-data.frame(rbind(clarkfork_three_agg_v1,clarkfork_three_rand_v1,clarkfork_three_seg_v1))
  clarkfork_v1<-cbind(clarkforknames_v1,clarkfork_v1)
  
  clarkfork_three_agg_v2<-data.frame(clarkfork_three_agg_v2)
  colnames(clarkfork_three_agg_v2)<-c("n")
  clarkfork_three_rand_v2<-data.frame(clarkfork_three_rand_v2)
  colnames(clarkfork_three_rand_v2)<-c("n")
  clarkfork_three_seg_v2<-data.frame(clarkfork_three_seg_v2)
  colnames(clarkfork_three_seg_v2)<-c("n")
  
  clarkforknames_v2<-c(rep("agg",nrow(clarkfork_three_agg_v2)),rep("rand",nrow(clarkfork_three_rand_v2)),
                       rep("seg",nrow(clarkfork_three_seg_v2)))
  clarkfork_v2<-data.frame(rbind(clarkfork_three_agg_v2,clarkfork_three_rand_v2,clarkfork_three_seg_v2))
  clarkfork_v2<-cbind(clarkforknames_v2,clarkfork_v2)
  
  Was_zero_agg_v1<-data.frame(Was_zero_agg_v1)
  colnames(Was_zero_agg_v1)<-c("n")
  Was_zero_rand_v1<-data.frame(Was_zero_rand_v1)
  colnames(Was_zero_rand_v1)<-c("n")
  Was_zero_seg_v1<-data.frame(Was_zero_seg_v1)
  colnames(Was_zero_seg_v1)<-c("n")
  
  Was_zeronames_v1<-c(rep("agg",nrow(Was_zero_agg_v1)),rep("rand",nrow(Was_zero_rand_v1)),
                      rep("seg",nrow(Was_zero_seg_v1)))
  Was_zero_v1<-data.frame(rbind(Was_zero_agg_v1,Was_zero_rand_v1,Was_zero_seg_v1))
  Was_zero_v1<-cbind(Was_zeronames_v1,Was_zero_v1)
  
  Was_zero_agg_v2<-data.frame(Was_zero_agg_v2)
  colnames(Was_zero_agg_v2)<-c("n")
  Was_zero_rand_v2<-data.frame(Was_zero_rand_v2)
  colnames(Was_zero_rand_v2)<-c("n")
  Was_zero_seg_v2<-data.frame(Was_zero_seg_v2)
  colnames(Was_zero_seg_v2)<-c("n")
  
  Was_zeronames_v2<-c(rep("agg",nrow(Was_zero_agg_v2)),rep("rand",nrow(Was_zero_rand_v2)),
                      rep("seg",nrow(Was_zero_seg_v2)))
  Was_zero_v2<-data.frame(rbind(Was_zero_agg_v2,Was_zero_rand_v2,Was_zero_seg_v2))
  Was_zero_v2<-cbind(Was_zeronames_v2,Was_zero_v2)
  
  Was_one_agg_v1<-data.frame(Was_one_agg_v1)
  colnames(Was_one_agg_v1)<-c("n")
  Was_one_rand_v1<-data.frame(Was_one_rand_v1)
  colnames(Was_one_rand_v1)<-c("n")
  Was_one_seg_v1<-data.frame(Was_one_seg_v1)
  colnames(Was_one_seg_v1)<-c("n")
  
  Was_onenames_v1<-c(rep("agg",nrow(Was_one_agg_v1)),rep("rand",nrow(Was_one_rand_v1)),
                     rep("seg",nrow(Was_one_seg_v1)))
  Was_one_v1<-data.frame(rbind(Was_one_agg_v1,Was_one_rand_v1,Was_one_seg_v1))
  Was_one_v1<-cbind(Was_onenames_v1,Was_one_v1)
  
  Was_one_agg_v2<-data.frame(Was_one_agg_v2)
  colnames(Was_one_agg_v2)<-c("n")
  Was_one_rand_v2<-data.frame(Was_one_rand_v2)
  colnames(Was_one_rand_v2)<-c("n")
  Was_one_seg_v2<-data.frame(Was_one_seg_v2)
  colnames(Was_one_seg_v2)<-c("n")
  
  Was_onenames_v2<-c(rep("agg",nrow(Was_one_agg_v2)),rep("rand",nrow(Was_one_rand_v2)),
                     rep("seg",nrow(Was_one_seg_v2)))
  Was_one_v2<-data.frame(rbind(Was_one_agg_v2,Was_one_rand_v2,Was_one_seg_v2))
  Was_one_v2<-cbind(Was_onenames_v2,Was_one_v2)
  
  results.final.v1<-list(clarkfork_v1,Was_one_v1,Was_zero_v1)
  results.final.v2<-list(clarkfork_v2,Was_one_v1,Was_zero_v2)
  results.final<-list(results.final.v1,results.final.v2)
  return(results.final)
}

####