
# Null models to test for change through time

library(here)
library(cooccur)

petm_all<-read.csv(here("Data/petm_occur.csv"),header=T)

# Clarkfork3


clark3_localities<-c("SC-8",
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

Was0_localities<-c("WW-71",
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


Was1_localities<-c("WW-89",
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

#occurrences?
all_locs<-c(clark3_localities,Was0_localities,Was1_localities)
bins<-list(clark3_localities,Was0_localities,Was1_localities)
matches<-intersect(petm_all$collection_name,all_locs)

occurrences<-petm_all[matches,]
occurrences<-occurrences[rowSums(occurrences)>0,]
occurrences<-occurrences[,colSums(occurrences)>0]

clim_prefs<-read.csv(here("Data/Mammal climate prefs new.csv"),header=T,row.names = 1)
BMs<-read.csv(here("Data/All masses combined.csv"),header=T,row.names=1)

null_results<-Shuffletaxa(occurences,clim_prefs,BMs,loco,niter=1000,summarize="FALSE")

###





