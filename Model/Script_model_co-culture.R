#install package BacArena, sybil, glpkAPI, data.table

library(BacArena)
library(sybil)
library(glpkAPI)
library(data.table)
sybil::SYBIL_SETTINGS("SOLVER","glpkAPI")

#read input file
#medium <-read.csv("E:/Đề tài/GEM đi thi/gapseq-master/dat/media/TSBmed.csv")
pa <- readRDS(".../pa.RDS") #directory file in your computer
sa <- readRDS(".../sa.RDS") #directory file in your computer

#from data to bacterial
bacpa <- Bac(pa)
bacsa <- Bac(sa)

#setting enviroment culture
arena = Arena(n=100, m=100, Lx=0.025, Ly=0.025)

arena <- addOrg(arena, bacpa, amount=1)
arena <- addOrg(arena, bacsa, amount=1)

arena_subs <- fread(".../TSBmed.csv")
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]
arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)


#start simulation
CF_sim <- simEnv(arena,time=13, sec_obj = "mtf")


#show graph
plotGrowthCurve(sim)
plotSpecActivity(sim)[2]
