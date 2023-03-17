#install package BacArena, sybil, glpkAPI, data.table

library(BacArena)
library(sybil)
library(glpkAPI)
library(data.table)
sybil::SYBIL_SETTINGS("SOLVER","glpkAPI")

#read input file
#medium <-read.csv("E:/Đề tài/GEM đi thi/gapseq-master/dat/media/TSBmed.csv")
pa <- readRDS("C:/Khiem/Model BacArena/pa.RDS") #directory file in your computer
sa <- readRDS("C:/Khiem/Model BacArena/sa.RDS") #directory file in your computer

#from data to bacterial
bacpa <- Bac(pa, type="Pseudomonas aeruginosa")
bacsa <- Bac(sa, type="MRSA")

#setting enviroment culture
arena = Arena(n=100, m=100, Lx=0.025, Ly=0.025)

arena <- addOrg(arena, bacpa, amount=1)
arena <- addOrg(arena, bacsa, amount=6)

arena_subs <- fread("C:/Khiem/Model BacArena/TSBmed.csv")
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]
arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)


#start simulation
sim <- simEnv(arena,time=10)

par(mfrow=c(3,2))
evalArena(sim, plot_items = "Population", phencol = F, time = seq(1,9,2), 
          show_legend = FALSE, legend_pos = "bottom")

plotCurves2(sim, legendpos = "right")

#show graph
par(mfrow=c(2,3))
evalArena(sim, show_legend = FALSE, time=seq(1,5,1))

plotCurves2(sim)
plotGrowthCurve(CF_sim)
plotSpecActivity(CF_sim)[2]
