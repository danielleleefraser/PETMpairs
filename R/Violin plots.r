

#assemble the data from cooccur analysis.r

clark_agg<-Clark3clim_comp$agg_v1
clark_agg<-data.frame(clark_agg)
clark1_agg<-cbind("Clark3",clark_agg)
colnames(clark1_agg)<-c("NALMA","Aggs")
Was0_agg<-Was0clim_comp$agg_v1
Was0_agg<-data.frame(Was0_agg)
Was0_agg<-cbind("Was0",Was0_agg)
colnames(Was0_agg)<-c("NALMA","Aggs")
Was1_agg<-Was1clim_comp$agg_v1
Was1_agg<-data.frame(Was1_agg)
Was1_agg<-cbind("Was1",Was1_agg)
colnames(Was1_agg)<-c("NALMA","Aggs")
aggregates<-rbind(clark1_agg,Was0_agg,Was1_agg)

clark_agg<-Clark3clim_comp$agg_v2
clark_agg<-data.frame(clark_agg)
clark1_agg<-cbind("Clark3",clark_agg)
colnames(clark1_agg)<-c("NALMA","Aggs")
Was0_agg<-Was0clim_comp$agg_v2
Was0_agg<-data.frame(Was0_agg)
Was0_agg<-cbind("Was0",Was0_agg)
colnames(Was0_agg)<-c("NALMA","Aggs")
Was1_agg<-Was1clim_comp$agg_v2
Was1_agg<-data.frame(Was1_agg)
Was1_agg<-cbind("Was1",Was1_agg)
colnames(Was1_agg)<-c("NALMA","Aggs")
aggregates<-rbind(clark1_agg,Was0_agg,Was1_agg)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

library(ggplot2)
p<-ggplot(aggregates, aes(x=NALMA, y=Aggs)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  #ylim(0.0,0.4)+
  geom_violin()+ 
  stat_summary(fun=median, geom="point", size=2, color="red")+
  stat_summary(fun.data=data_summary)+
  geom_jitter(shape=16, position=position_jitter(0.02))

# Body mass

clark_agg<-Clark3bm_comp$agg
clark_agg<-data.frame(clark_agg)
clark1_agg<-cbind("Clarkfork",clark_agg)
colnames(clark1_agg)<-c("NALMA","Aggs")
Was0_agg<-Was0bm_comp$agg
Was0_agg<-data.frame(Was0_agg)
Was0_agg<-cbind("Wasatch0",Was0_agg)
colnames(Was0_agg)<-c("NALMA","Aggs")
Was1_agg<-Was1bm_comp$agg
Was1_agg<-data.frame(Was1_agg)
Was1_agg<-cbind("Wasatch1",Was1_agg)
colnames(Was1_agg)<-c("NALMA","Aggs")
aggregates<-rbind(clark1_agg,Was0_agg,Was1_agg)

clark_seg<-Clark3bm_comp$seg
clark_seg<-data.frame(clark_seg)
clark1_seg<-cbind("Clarkfork",clark_seg)
colnames(clark1_seg)<-c("NALMA","segs")
Was0_seg<-Was0bm_comp$seg
Was0_seg<-data.frame(Was0_seg)
Was0_seg<-cbind("Wasatch0",Was0_seg)
colnames(Was0_seg)<-c("NALMA","segs")
Was1_seg<-Was1bm_comp$seg
Was1_seg<-data.frame(Was1_seg)
Was1_seg<-cbind("Wasatch1",Was1_seg)
colnames(Was1_seg)<-c("NALMA","segs")
segregates<-rbind(clark1_seg,Was0_seg,Was1_seg)

clark_rand<-Clark3bm_comp$rand
clark_rand<-data.frame(clark_rand)
clark1_rand<-cbind("Clarkfork",clark_rand)
colnames(clark1_rand)<-c("NALMA","rands")
Was0_rand<-Was0bm_comp$rand
Was0_rand<-data.frame(Was0_rand)
Was0_rand<-cbind("Wasatch0",Was0_rand)
colnames(Was0_rand)<-c("NALMA","rands")
Was1_rand<-Was1bm_comp$rand
Was1_rand<-data.frame(Was1_rand)
Was1_rand<-cbind("Wasatch1",Was1_rand)
colnames(Was1_rand)<-c("NALMA","rands")
randregates<-rbind(clark1_rand,Was0_rand,Was1_rand)

library(ggplot2)

p<-ggplot(aggregates, aes(x=NALMA, y=Aggs)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  ylim(0.0,0.4)+
  geom_violin()+ 
  stat_summary(fun=median, geom="point", size=2, color="red")+
  stat_summary(fun.data=data_summary)+
  geom_jitter(shape=16, position=position_jitter(0.02))



##
