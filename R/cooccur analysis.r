

### time vs NALMA co-occurrence analysis

library(here)

petm_occur<-read.csv(here("Data/petm_occur.csv"),header=T)

source("Acessory functions.r")

# Clarkfork3

library(cooccur)

localities<-c("SC-8",
              "SC-9",
              "SC-10",
              "SC-11",
              "SC-21",
              "SC-22",
              "SC-23",
              "SC-24",
              "SC-25",
              "SC-28",
              "SC-29",
              "SC-48",
              "SC-49",
              "SC-50",
              "SC-51",
              "SC-52",
              "SC-55",
              "SC-56",
              "SC-57",
              "SC-59",
              "SC-60",
              "SC-70",
              "SC-71",
              "SC-72",
              "SC-73",
              "SC-75",
              "SC-76",
              "SC-77",
              "SC-80",
              "SC-81",
              "SC-90",
              "SC-100",
              "SC-101",
              "SC-102",
              "SC-103",
              "SC-105",
              "SC-106",
              "SC-107",
              "SC-138",
              "SC-140",
              "SC-149",
              "SC-150",
              "SC-152",
              "SC-153",
              "SC-154",
              "SC-155",
              "SC-158",
              "SC-159",
              "SC-162",
              "SC-163",
              "SC-164",
              "SC-175",
              "SC-176",
              "SC-183",
              "SC-184",
              "SC-202",
              "SC-203",
              "SC-204",
              "SC-214",
              "SC-230",
              "SC-231",
              "SC-233",
              "SC-234",
              "SC-235",
              "SC-141",
              "SC-350",
              "SC-404",
              "SC-343")

# Create the occurrence matrix

clark3<-petm_occur[petm_occur$collection_name==localities[1],]

for(i in 2:length(localities)){
  clark_temp<-petm_occur[petm_occur$collection_name==localities[i],]
  clark3<-rbind(clark3,clark_temp)
}

table<-xtabs(~clark3$collection_name+clark3$accepted_name)
table<-apply(table,c(1,2),function(x) {ifelse(any(x>0),1,0)})
clark3_occur<-t(table)

# cooccur analyses

cooccur.clark3<-pairs_calc(clark3_occur)

clim_prefs<-read.csv(here("Data/Mammal climate prefs new.csv"),header=T,row.names = 1)
BMs<-read.csv(here("Data/All masses combined Feb 2017.csv"),header=T,row.names=1)
loco<-read.csv(here("Data/LocomotorData.csv"),header=T,row.names=1)

# Compare climate preferences of pairs
Clark3clim_comp<-compare_climpref(cooccur.clark3,clim_prefs,calc.means=TRUE)
colnames(Clark3clim_comp)<-c("DiffposNMDS1","DiffnegNMDS1","DiffrandNMDS1","sdposNMDS1","sdnegNMDS1","sdrandNMDS1","DiffposNMDS2",
                             "DiffnegNMDS2","DiffrandNMDS2","sdposNMDS2","sdnegNMDS2","sdrandNMDS2")

# Compare masses of pairs
Clark3BM_comp<-compare_mass(cooccur.clark3,BMs,calc.means = TRUE)
colnames(Clark3BM_comp)<-c("DiffposBM","DiffnegBM","DiffrandBM","sefposBM","senegBM","serandBM")

# Compare locomotor mode of pairs 
Clark3loco_comp<-compare_loco(cooccur.clark3,loco,calc.means = TRUE)
colnames(Clark3loco_comp)<-c("Diffposloco","Diffnegloco","Diffrandloco","sdfposloco","sdnegloco","sdrandloco")

# Repeat for other two NALMAs

# Was0

localities<-c("WW-71",
              "WW-72",
              "WW-73",
              "WW-74",
              "WW-75",
              "WW-76",
              "WW-77",
              "WW-78",
              "WW-79",
              "WW-80",
              "WW-83",
              "WW-84",
              "WW-85",
              "WW-86",
              "WW-87",
              "WW-88",
              "WW-90",
              "WW-91",
              "WW-96",
              "WW-97",
              "WW-98",
              "WW-99",
              "WW-101",
              "WW-112",
              "WW-113",
              "WW-114",
              "WW-115",
              "WW-117",
              "WW-118",
              "WW-119",
              "WW-125",
              "WW-126",
              "WW-128",
              "WW-171",
              "WW-172",
              "WW-184",
              "WW-186",
              "WW-191",
              "WW-192",
              "SC-67",
              "SC-69",
              "SC-121",
              "SC-308",
              "SC-345",
              "SC-349",
              "SC-351",
              "SC-405",
              "SC-79")

was0<-petm_occur[petm_occur$collection_name==localities[1],]

for(i in 2:length(localities)){
  was0_temp<-petm_occur[petm_occur$collection_name==localities[i],]
  was0<-rbind(was0,was0_temp)
}

table<-xtabs(~was0$collection_name+was0$accepted_name)
table<-apply(table,c(1,2),function(x) {ifelse(any(x>0),1,0)})
Was0_occur<-t(table)

cooccur.Was0<-pairs_calc(Was0_occur)

Was0clim_comp<-compare_climpref(cooccur.Was0,clim_prefs,calc.means=TRUE)
colnames(Was0clim_comp)<-c("DiffposNMDS1","DiffnegNMDS1","DiffrandNMDS1","sdposNMDS1","sdnegNMDS1","sdrandNMDS1","DiffposNMDS2",
                             "DiffnegNMDS2","DiffrandNMDS2","sdposNMDS2","sdnegNMDS2","sdrandNMDS2")


Was0BM_comp<-compare_mass(cooccur.Was0,BMs,calc.means = TRUE)
colnames(Was0BM_comp)<-c("DiffposBM","DiffnegBM","DiffrandBM","sefposBM","senegBM","serandBM")

Was0loco_comp<-compare_loco(cooccur.Was0,loco,calc.means = TRUE)
colnames(Was0loco_comp)<-c("Diffposloco","Diffnegloco","Diffrandloco","sdfposloco","sdnegloco","sdrandloco")

# WASOne

localities<-c("WW-89",
              "WW-173",
              "WW-178",
              "SC-40",
              "SC-182",
              "SC-4",
              "SC-6",
              "SC-40",
              "SC-89",
              "SC-123",
              "SC-253",
              "SC-295",
              "SC-128",
              "SC-192",
              "SC-68",
              "SC-122",
              "SC-206",
              "SC-142",
              "SC-44",
              "SC-17",
              "SC-18",
              "SC-16",
              "SC-37",
              "SC-210",
              "SC-54",
              "SC-2",
              "SC-12",
              "SC-6",
              "SC-4",
              "SC-129",
              "SC-47",
              "SC-46")

wasone<-petm_occur[petm_occur$collection_name==localities[1],]

for(i in 2:length(localities)){
  wasone_temp<-petm_occur[petm_occur$collection_name==localities[i],]
  wasone<-rbind(wasone,wasone_temp)
}

table<-xtabs(~wasone$collection_name+wasone$accepted_name)
table<-apply(table,c(1,2),function(x) {ifelse(any(x>0),1,0)})
Was1_occur<-t(table)

cooccur.Was1<-pairs_calc(Was1_occur)

Was1clim_comp<-compare_climpref(cooccur.Was1,clim_prefs,calc.means=TRUE)
colnames(Was1clim_comp)<-c("DiffposNMDS1","DiffnegNMDS1","DiffrandNMDS1","sdposNMDS1","sdnegNMDS1","sdrandNMDS1","DiffposNMDS2",
                           "DiffnegNMDS2","DiffrandNMDS2","sdposNMDS2","sdnegNMDS2","sdrandNMDS2")


Was1BM_comp<-compare_mass(cooccur.Was1,BMs,calc.means = TRUE)
colnames(Was1BM_comp)<-c("DiffposBM","DiffnegBM","DiffrandBM","sefposBM","senegBM","serandBM")

Was1loco_comp<-compare_loco(cooccur.Was1,loco,calc.means = TRUE)
colnames(Was1loco_comp)<-c("Diffposloco","Diffnegloco","Diffrandloco","sdfposloco","sdnegloco","sdrandloco")


# Plot changes in numbers of pairs

Clark3.nums<-c(nrow(cooccur.clark3[[1]]),nrow(cooccur.clark3[[2]]),nrow(cooccur.clark3[[3]]))
Was0.nums<-c(nrow(cooccur.Was0[[1]]),nrow(cooccur.Was0[[2]]),nrow(cooccur.Was0[[3]]))
Was1.nums<-c(nrow(cooccur.Was1[[1]]),nrow(cooccur.Was1[[2]]),nrow(cooccur.Was1[[3]]))

ForPlot<-rbind(Clark3.nums,Was0.nums,Was1.nums)
colnames(ForPlot)<-c("Agg","Seg","Rand")
Ages<-c(1:3)
ForPlot<-cbind(Ages,ForPlot)
ForPlot<-data.frame(ForPlot)
Numtot<-rowSums(ForPlot[,2:3])
ForPlot<-cbind(ForPlot,Numtot)

cols <- c("Agg"="tan4","Seg"="tan1","Rand"="darkorange3","Total"="burlywood")

library(ggplot2)
ggplot(ForPlot, aes(x=Ages,y=Agg))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  #scale_x_reverse(breaks=c(3,2,1),labels=c("","25 ka-20 ka","20 ka-15 ka","15 ka-10 ka","10 ka - 5 ka","5 ka - modern","Modern"))+  
  #ylim(0.0,4)+
  xlab("NALMA")+
  ylab("Number of Pairs")+
  geom_point(aes(x=Ages,y=Agg,colour="Agg"))+
  geom_line(aes(x=Ages,y=Agg,colour="Agg"),size=2)+
  geom_point(aes(x=Ages,y=Seg,colour="Seg"))+
  geom_line(aes(x=Ages,y=Seg,colour="Seg"),size=2)+
  geom_point(aes(x=Ages,y=Rand,colour="Rand"))+
  geom_line(aes(x=Ages,y=Rand,colour="Rand"),size=2)+
  geom_point(aes(x=Ages,y=Numtot,colour="Total"))+
  geom_line(aes(x=Ages,y=Numtot,colour="Total"),size=2)+
  scale_colour_manual(name="Type of Pair",values=cols, guide = guide_legend(fill = NULL,colour = NULL)) + 
  theme_classic()


##