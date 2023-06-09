#install package BacArena, sybil, glpkAPI, data.table

#load library
library(BacArena)
library(sybil)
library(glpkAPI)
library(data.table)
sybil::SYBIL_SETTINGS("SOLVER","glpkAPI")

#input genome metabolic model
pa <- readRDS("/pa.RDS") #directory file in your computer
sa <- readRDS("/sa.RDS") #directory file in your computer

#from data to bacterial
bacpa <- Bac(pa, type="Pseudomonas aeruginosa", chem = "EX_cpd00007_e0")
bacsa <- Bac(sa, type="MRSA")

#setting enviroment culture
arena = Arena(n=100, m=100, Lx=0.025, Ly=0.025)

#add bacteria
arena <- addOrg(arena, bacpa, amount=15)
arena <- addOrg(arena, bacsa, amount=85)

arena_subs <- fread("/TSBmed.csv") #dỉrectory
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]
arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
arena <- createGradient(arena,smax=50,mediac=c("EX_cpd00007_e0"),
                    position='top',steep=0.5, add=TRUE)
chemotaxis(bacpa,arena,1, "EX_cpd00007_e0", arena@occupyM)

#start simulation
sim <- simEnv(arena,time=24)

#show result
evalArena(sim, plot_items = "Population", phencol = F, time = seq(0,18,1), 
          show_legend = FALSE, legend_pos = "bottom")

plotCurves2(sim, legendpos = "right")

#end