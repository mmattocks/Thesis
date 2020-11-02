using Plots, Images, FileIO, Plots.PlotMeasures

sib36=load("/bench/PhD/Thesis/images/rys/C2 cell 1/Rys sib C2 36k cell 1.jpg")
sib36=sib36[1:2465,:]
sib110=load("/bench/PhD/Thesis/images/rys/C2 cell 1/Rys sib C2 110k cell 1.jpg")
sib110=sib110[1:2465,:]
sib210=load("/bench/PhD/Thesis/images/rys/C2 cell 1/Rys sib C2 210k cell 1.jpg")
sib210=sib210[1:2465,:]

rys36=load("/bench/PhD/Thesis/images/rys/D1 cell 1/Rys mut D1 CMZ 36k.jpg")
rys36=rys36[1:2465,:]
rys110=load("/bench/PhD/Thesis/images/rys/D1 cell 1/Rys mut D1 CMZ 110k.jpg")
rys110=rys110[1:2465,:]
rys210=load("/bench/PhD/Thesis/images/rys/D1 cell 1/Rys mut D1 CMZ 210k.jpg")
rys210=rys210[1:2465,:]


s36p=plot(sib36, axis=nothing, border=:none, margin=0mm)
annotate!([(200,105,Plots.text("SIB",25,:white,:bold))])
s110p=plot(sib110, axis=nothing,  border=:none, margin=0mm)
annotate!([(200,105,Plots.text("SIB",25,:white,:bold))])
s210p=plot(sib210, axis=nothing, border=:none, margin=0mm)
annotate!([(200,105,Plots.text("SIB",25,:white,:bold))])
r36p=plot(rys36, axis=nothing, border=:none, margin=0mm)
annotate!([(200,105,Plots.text("RYS",25,:white,:bold))])
r110p=plot(rys110, axis=nothing,  border=:none, margin=0mm)
annotate!([(200,105,Plots.text("RYS",25,:white,:bold))])
r210p=plot(rys210, axis=nothing, border=:none, margin=0mm)
annotate!([(200,105,Plots.text("RYS",25,:white,:bold))])

combined=plot(s36p, r36p, s110p, r110p, s210p, r210p, layout=grid(3,2),size=(1200,2400))

savefig(combined, "/bench/PhD/Thesis/images/emtogimp.png")