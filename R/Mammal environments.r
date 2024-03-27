
library(here)

mammals<-read.csv(here("Data/Mammal lat longs.csv",header=T),row.names=2)
plants<-read.csv(here("Data/Plant_latlong_new.csv"),header=T,row.names=1)

mam_plant<-rbind(mammals,plants)

#Find distances among all the sites

library(fossil)

#do for each NALMA
NALMAs<-unique(mam_plant$NALMA)
colnames(mam_plant)[2:3]<-c("longitude","latitude")

temp<-mam_plant[mam_plant$NALMA==NALMAs[1],]
temp[,1]<-NULL
dists_mp<-earth.dist(temp)
dists_mp<-as.matrix(dists_mp)
colnames(dists_mp)<-rownames(temp)
rownames(dists_mp)<-rownames(temp)
# take out the ones that compare mammals to plants
mamm_temp<-mammals[mammals$NALMA==NALMAs[1],]
dists_mpnew<-dists_mp[rownames(mamm_temp),]
plant_temp<-plants[plants$NALMA==NALMAs[1],]
dists_mpnew<-dists_mpnew[,rownames(plant_temp)]
# for each mammal site find the closest plant site
closest<-matrix(nrow=nrow(dists_mpnew),ncol=1)
rownames(closest)<-rownames(dists_mpnew)
for(j in 1:nrow(dists_mpnew)){
  closest[j,1]<-names(which(dists_mpnew[j,]==min(dists_mpnew[j,])))
}

for(i in 2:3){
  temp<-mam_plant[mam_plant$NALMA==NALMAs[i],]
  temp[,1]<-NULL
  dists_mp<-earth.dist(temp)
  dists_mp<-as.matrix(dists_mp)
  colnames(dists_mp)<-rownames(temp)
  rownames(dists_mp)<-rownames(temp)
  # take out the ones that compare mammals to plants
  mamm_temp<-mammals[mammals$NALMA==NALMAs[i],]
  dists_mpnew<-dists_mp[rownames(mamm_temp),]
  plant_temp<-plants[plants$NALMA==NALMAs[i],]
  dists_mpnew<-dists_mpnew[,rownames(plant_temp)]
 # for each mammal site find the closest plant site
  closest_temp<-matrix(nrow=nrow(dists_mpnew),ncol=1)
  rownames(closest_temp)<-rownames(dists_mpnew)
  for(j in 1:nrow(dists_mpnew)){
    closest_temp[j,1]<-names(which(dists_mpnew[j,]==min(dists_mpnew[j,])))
  }
  closest<-rbind(closest,closest_temp)
}

write.csv(closest,here("Data/Closest mammal and fossil sites_new.csv"))

# Now find the NMDS1 and 2

plant_NMDS<-read.csv(here("Data/Pollen scores new.csv"),header=T,row.names = 1)
mammal_clim<-matrix(nrow=nrow(closest),ncol=2)
rownames(mammal_clim)<-rownames(closest)

for(i in 1:nrow(closest)){
  print(i)
  mammal_clim[i,1]<-plant_NMDS[which(rownames(plant_NMDS)==closest[i,1]),1]
  mammal_clim[i,2]<-plant_NMDS[which(rownames(plant_NMDS)==closest[i,1]),2]
}

write.csv(mammal_clim,here("Data/Mammal site climates new.csv"))

############################################################################################################################

# Now assign each species a climate envelope based on the sites they occur at

library(cooccur)

petm_occur<-read.csv(here("Data/petm_occur.csv"),header=T)

# Clarkfork3

clark3<-c("SC-8",
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

# Was0

Was0<-c("WW-71",
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

# WASOne

Was1<-c("WW-89",
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

localities<-c(clark3,Was0,Was1)

mammal_clim<-read.csv(here("Data/Mammal site climates new.csv"),header=T,row.names=1)
mammal_pref<-matrix(nrow=nrow(all_table),ncol=2)
rownames(mammal_pref)<-rownames(all_table)

for(i in 1:nrow(all_table)){
  temp1<-all_table[i,all_table[i,1:ncol(all_table)]>0]
  if(length(temp1)>1){
    mammal_pref[i,1]<-mean(mammal_clim[names(temp1),][,1])
    mammal_pref[i,2]<-mean(mammal_clim[names(temp1),][,2])
  }
  if(length(temp1)==1){
    site<-which(all_table[i,1:ncol(all_table)]>0)
    mammal_pref[i,1]<-mean(mammal_clim[names(site),][,1])
    mammal_pref[i,2]<-mean(mammal_clim[names(site),][,2])
  }
}

write.csv(mammal_pref,here("Data/Mammal climate prefs new.csv"))


# Assess change in mean climate preference through time

climpref<-read.csv(here("Data/Mammal climate prefs new.csv"),header=T,row.names=1)
BMs<-read.csv(here("Data/All masses combined.csv"),header=T,row.names=1)

# Clarkfork3

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

clark3<-petm_occur[petm_occur$collection_name==localities[1],]

for(i in 2:length(localities)){
  clark_temp<-petm_occur[petm_occur$collection_name==localities[i],]
  clark3<-rbind(clark3,clark_temp)
}

table<-xtabs(~clark3$collection_name+clark3$accepted_name)
table<-apply(table,c(1,2),function(x) {ifelse(any(x>0),1,0)})
clark3_occur<-t(table)

# Mean climate preference

clark3_clims<-climpref[rownames(clark3_occur),]
clark3_BM<-BMs[rownames(clark3_occur),]
# make new object to put means in

mean_clims<-matrix(nrow=3,ncol=6)
colnames(mean_clims)<-c("MeanV1","SDV1","SEV1","MeanV2","SDV2","SEV2")
mean_clims[1,1]<-mean(clark3_clims$V1)
mean_clims[1,2]<-sd(clark3_clims$V1)
mean_clims[1,3]<-sd(clark3_clims$V1)/sqrt(nrow(clark3_clims))
mean_clims[1,4]<-mean(clark3_clims$V2)
mean_clims[1,5]<-sd(clark3_clims$V2)
mean_clims[1,6]<-sd(clark3_clims$V2)/sqrt(nrow(clark3_clims))

mean_BM<-matrix(nrow=3,ncol=3) # is this in grams? yes, check on methods for masses being combined
colnames(mean_BM)<-c("MeanBM","SDBM","SEBM")
mean_BM[1,1]<-mean(na.omit(clark3_BM$ln_mass))
mean_BM[1,2]<-sd(na.omit(clark3_BM$ln_mass))
mean_BM[1,3]<-sd(na.omit(clark3_BM$ln_mass))/sqrt(nrow(clark3_BM))


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

# Mean climate preference

Was0_clims<-climpref[rownames(Was0_occur),]
Was0_BM<-BMs[rownames(Was0_occur),]
# make new object to put means in

mean_clims[2,1]<-mean(Was0_clims$V1)
mean_clims[2,2]<-sd(Was0_clims$V1)
mean_clims[2,3]<-sd(Was0_clims$V1)/sqrt(nrow(Was0_clims))
mean_clims[2,4]<-mean(Was0_clims$V2)
mean_clims[2,5]<-sd(Was0_clims$V2)
mean_clims[2,6]<-sd(Was0_clims$V2)/sqrt(nrow(Was0_clims))


colnames(mean_BM)<-c("MeanBM","SDBM","SEBM")
mean_BM[2,1]<-mean(na.omit(Was0_BM$ln_mass))
mean_BM[2,2]<-sd(na.omit(Was0_BM$ln_mass))
mean_BM[2,3]<-sd(na.omit(Was0_BM$ln_mass))/sqrt(nrow(Was0_BM))


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

Was1_clims<-climpref[rownames(Was1_occur),]
Was1_BM<-BMs[rownames(Was1_occur),]
# make new object to put means in

mean_clims[3,1]<-mean(Was1_clims$V1)
mean_clims[3,2]<-sd(Was1_clims$V1)
mean_clims[3,3]<-sd(Was1_clims$V1)/sqrt(nrow(Was1_clims))
mean_clims[3,4]<-mean(Was1_clims$V2)
mean_clims[3,5]<-sd(Was1_clims$V2)
mean_clims[3,6]<-sd(Was1_clims$V2)/sqrt(nrow(Was1_clims))


colnames(mean_BM)<-c("MeanBM","SDBM","SEBM")
mean_BM[3,1]<-mean(na.omit(Was1_BM$ln_mass))
mean_BM[3,2]<-sd(na.omit(Was1_BM$ln_mass))
mean_BM[3,3]<-sd(na.omit(Was1_BM$ln_mass))/sqrt(nrow(Was1_BM))



#### Without newly appearing taxa during the PETM

# Clarkfork3

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

clark3<-petm_occur[petm_occur$collection_name==localities[1],]

for(i in 2:length(localities)){
  clark_temp<-petm_occur[petm_occur$collection_name==localities[i],]
  clark3<-rbind(clark3,clark_temp)
}

table<-xtabs(~clark3$collection_name+clark3$accepted_name)
table<-apply(table,c(1,2),function(x) {ifelse(any(x>0),1,0)})
clark3_occur<-t(table)

# Mean climate preference

clark3_clims<-climpref[rownames(clark3_occur),]
clark3_clims<-cbind(clark3_clims,"Clark3")
colnames(clark3_clims)<-c("V1","V2","NALMA")
clark3_BM<-BMs[rownames(clark3_occur),]
clark3_BM<-cbind(clark3_BM,"Clark3")
colnames(clark3_BM)<-c("V1","V2","NALMA")
# make new object to put means in

mean_clims<-matrix(nrow=3,ncol=6)
colnames(mean_clims)<-c("MeanV1","SDV1","SEV1","MeanV2","SDV2","SEV2")
mean_clims[1,1]<-mean(clark3_clims$V1)
mean_clims[1,2]<-sd(clark3_clims$V1)
mean_clims[1,3]<-sd(clark3_clims$V1)/sqrt(nrow(clark3_clims))
mean_clims[1,4]<-mean(clark3_clims$V2)
mean_clims[1,5]<-sd(clark3_clims$V2)
mean_clims[1,6]<-sd(clark3_clims$V2)/sqrt(nrow(clark3_clims))

mean_BM<-matrix(nrow=3,ncol=3) # is this in grams? yes, check on methods for masses being combined
colnames(mean_BM)<-c("MeanBM","SDBM","SEBM")
mean_BM[1,1]<-mean(na.omit(clark3_BM$ln_mass))
mean_BM[1,2]<-sd(na.omit(clark3_BM$ln_mass))
mean_BM[1,3]<-sd(na.omit(clark3_BM$ln_mass))/sqrt(nrow(clark3_BM))


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
#table<-t(table)
#table[1:nrow(table),rownames(firstocc)]<-0
Was0_occur<-t(table)

# First appearances in Was0

diffs<-setdiff(rownames(Was0_occur),rownames(clark3_occur))
diffs2<-setdiff(rownames(Was0_occur),diffs)
#Was0_occur<-Was0_occur[diffs2,] # Taxa that do not first appear in Wasatch 0, selection of taxa with wider
# climate pref
Was0_occur<-Was0_occur[diffs,] # first appearances only


# Mean climate preference

Was0_clims<-climpref[rownames(Was0_occur),]
Was0_clims<-cbind(Was0_clims,"Was0")
colnames(Was0_clims)<-c("V1","V2","NALMA")
Was0_BM<-BMs[rownames(Was0_occur),]
Was0_BM<-cbind(Was0_BM,"Was0")
colnames(Was0_BM)<-c("V1","V2","NALMA")
# make new object to put means in

mean_clims[2,1]<-mean(Was0_clims$V1)
mean_clims[2,2]<-sd(Was0_clims$V1)
mean_clims[2,3]<-sd(Was0_clims$V1)/sqrt(nrow(Was0_clims))
mean_clims[2,4]<-mean(Was0_clims$V2)
mean_clims[2,5]<-sd(Was0_clims$V2)
mean_clims[2,6]<-sd(Was0_clims$V2)/sqrt(nrow(Was0_clims))


colnames(mean_BM)<-c("MeanBM","SDBM","SEBM")
mean_BM[2,1]<-mean(na.omit(Was0_BM$ln_mass))
mean_BM[2,2]<-sd(na.omit(Was0_BM$ln_mass))
mean_BM[2,3]<-sd(na.omit(Was0_BM$ln_mass))/sqrt(nrow(Was0_BM))

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

diffs2<-setdiff(rownames(Was1_occur),diffs)
Was1_occur<-Was1_occur[diffs2,]

Was1_clims<-climpref[rownames(Was1_occur),]
Was1_clims<-cbind(Was1_clims,"Was1")
colnames(Was1_clims)<-c("V1","V2","NALMA")
Was1_BM<-BMs[rownames(Was1_occur),]
Was1_BM<-cbind(Was1_BM,"Was1")
colnames(Was1_BM)<-c("V1","V2","NALMA")
# make new object to put means in

mean_clims[3,1]<-mean(Was1_clims$V1)
mean_clims[3,2]<-sd(Was1_clims$V1)
mean_clims[3,3]<-sd(Was1_clims$V1)/sqrt(nrow(Was1_clims))
mean_clims[3,4]<-mean(Was1_clims$V2)
mean_clims[3,5]<-sd(Was1_clims$V2)
mean_clims[3,6]<-sd(Was1_clims$V2)/sqrt(nrow(Was1_clims))


colnames(mean_BM)<-c("MeanBM","SDBM","SEBM")
mean_BM[3,1]<-mean(na.omit(Was1_BM$ln_mass))
mean_BM[3,2]<-sd(na.omit(Was1_BM$ln_mass))
mean_BM[3,3]<-sd(na.omit(Was1_BM$ln_mass))/sqrt(nrow(Was1_BM))

#ran similar code for plots with newly appearing species in was0

forPlot<-rbind(clark3_clims,Was0_clims,Was1_clims)
forPlot<-data.frame(forPlot)
#v2
p<-ggplot(forPlot, aes(x=NALMA, y=V2,color=NALMA,fill=NALMA)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  ylim(-0.2,0.2)+
  theme_classic()+
  geom_violin()+ 
  scale_color_manual(values=c("black", "black", "black"))+
  scale_fill_manual(values=c("salmon", "tan4", "purple"))+
  labs(y = "Mean Climate Preference (NMDS 2)")

p+
  #geom_jitter(data=randregates,aes(x=NALMA, y=n),position = position_jitter(width=0.02))+
  stat_summary(data=forPlot,aes(x=NALMA, y=V2),fun=median, geom="point", size=4, color="orange")+
  stat_summary(data=forPlot,aes(x=NALMA, y=V2),fun.data=data_summary,size=1,color="grey")

# Now do the same for body mass

forPlot<-rbind(clark3_BM,Was0_BM,Was1_BM)
forPlot<-data.frame(forPlot)
#v2
p<-ggplot(forPlot, aes(x=NALMA, y=V2,color=NALMA,fill=NALMA)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  #ylim(-0.2,0.2)+
  theme_classic()+
  geom_violin()+ 
  scale_color_manual(values=c("black", "black", "black"))+
  scale_fill_manual(values=c("salmon", "tan4", "purple"))+
  labs(y = "Mean Body Mass (g)")

p+
  #geom_jitter(data=randregates,aes(x=NALMA, y=n),position = position_jitter(width=0.02))+
  stat_summary(data=forPlot,aes(x=NALMA, y=V2),fun=median, geom="point", size=4, color="orange")+
  stat_summary(data=forPlot,aes(x=NALMA, y=V2),fun.data=data_summary,size=1,color="grey")



###