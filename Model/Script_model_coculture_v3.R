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
bacpa <- Bac(pa, type="Pseudomonas aeruginosa", chem = "EX_cpd00007_e0")
bacsa <- Bac(sa, type="MRSA")

#setting enviroment culture
arena = Arena(n=100, m=100, Lx=0.025, Ly=0.025)

arena <- addOrg(arena, bacpa, amount=15)
arena <- addOrg(arena, bacsa, amount=85)

arena_subs <- fread("C:/Khiem/Model BacArena/TSBmed.csv")
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]
arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
arena <- createGradient(arena,smax=75,mediac=c("EX_cpd00007_e0"),
                    position='top',steep=0.5, add=TRUE)
chemotaxis(bacpa,arena,1, "EX_cpd00007_e0", arena@occupyM)


#start simulation
sim <- simEnv(arena,time=18)

par(mfrow=c(1,1))
evalArena(sim, plot_items = "Population", phencol = F, time = seq(0,18,1), 
          show_legend = FALSE, legend_pos = "bottom")


plotCurves2(sim, legendpos = "right")

''
#show graph

