library(BacArena)
library(sybil)
library(glpkAPI)
library(data.table)
sybil::SYBIL_SETTINGS("SOLVER","glpkAPI")

#read input file
#medium <-read.csv("E:/Đề tài/GEM đi thi/gapseq-master/dat/media/TSBmed.csv")

sa <- readRDS("C:/Khiem/Model BacArena/sa.RDS") #directory file in your computer

#from data to bacterial

bacsa <- Bac(sa, type="MRSA")

#setting enviroment culture
arena = Arena(n=100, m=100, Lx=0.025, Ly=0.025)


arena <- addOrg(arena, bacsa, amount=20)

arena_subs <- fread("C:/Khiem/Model BacArena/TSBmed.csv")
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]
arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
arena <- createGradient(arena,smax=75,mediac=c("EX_cpd00007_e0"),
                        position='top',steep=0.5, add=TRUE)



#start simulation
sim <- simEnv(arena,time=10)

par(mfrow=c(1,1))
evalArena(sim, shadepalette(endcol = "red", incol="white"), plot_items = "Population", phencol = F, time = seq(0,10,1), 
          show_legend = FALSE, legend_pos = "bottom")


plotCurves2(sim, legendpos = "left")

''
#show graph

