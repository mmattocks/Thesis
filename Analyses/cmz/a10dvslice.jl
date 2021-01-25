using BayesianLinearRegression,CSV,DataFrames,Distributions,CMZNicheSims,GMC_NS, Serialization, Plots, KernelDensityEstimate,KernelDensityEstimatePlotting, Gadfly
a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"

edpath="/bench/PhD/NGS_binaries/CNS/A10/ed"
evpath="/bench/PhD/NGS_binaries/CNS/A10/ev"
edvpath="/bench/PhD/NGS_binaries/CNS/A10/edv"
paths=[edpath,evpath,edvpath]

default(legendfont = (10), guidefont = (12), tickfont = (10))

a10df=DataFrame(CSV.read(a10pth))
a10df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a10df)]

X=Vector{Float64}()
measure_dict=Dict{String,Vector{Vector{Float64}}}()

measure_dict["CMZ Sum"]=Vector{Vector{Float64}}()
measure_dict["Dorsal CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["Ventral CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["Lens circumference"]=Vector{Vector{Float64}}()

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for t_df in groupby(a10df, "Time point (d)")
    if !(Float64(unique(t_df."Time point (d)")[1]) in [180., 360.])

        push!(X,Float64(unique(t_df."Time point (d)")[1]))
        tvec,dvec,vvec,lcirc = [zeros(0) for i in 1:4]

        for i_df in groupby(t_df,"BSR")
            length([skipmissing(i_df."Lens diameter")...])>0 && push!(lcirc,mean(skipmissing(i_df."Lens diameter"))*π)
        end

        for i_df in groupby(t_df,"BSR")
            length([skipmissing(i_df."CMZ Sum")...])>0 && mean(skipmissing(i_df."CMZ Sum")) >0. && push!(tvec,mean(skipmissing(i_df."CMZ Sum")))
            length([skipmissing(i_df."Dorsal CMZ (#)")...])>0 && mean(skipmissing(i_df."Dorsal CMZ (#)")) >0. && push!(dvec,mean(skipmissing(i_df."Dorsal CMZ (#)")))
            length([skipmissing(i_df."Ventral CMZ (#)")...])>0 && mean(skipmissing(i_df."Ventral CMZ (#)")) >0. && push!(vvec,mean(skipmissing(i_df."Ventral CMZ (#)")))
        end

        push!(measure_dict["CMZ Sum"],tvec)
        push!(measure_dict["Dorsal CMZ (#)"],dvec)
        push!(measure_dict["Ventral CMZ (#)"],vvec)
        push!(measure_dict["Lens circumference"],lcirc)
    end
end

obsset=[measure_dict["Dorsal CMZ (#)"],measure_dict["Ventral CMZ (#)"],measure_dict["CMZ Sum"]]

dpopdist=fit(LogNormal,measure_dict["Dorsal CMZ (#)"][1])
vpopdist=fit(LogNormal,measure_dict["Ventral CMZ (#)"][1])
dvpopdist=fit(LogNormal,measure_dict["CMZ Sum"][1])

logLC=vcat([log.(measure_dict["Lens circumference"][n]) for n in 1:length(X)]...)
logX=vcat([[log.(X[n]) for ms in measure_dict["Lens circumference"][n]] for n in 1:length(X)]...)
logX=hcat(ones(length(logX)),logX)

lm=BayesianLinearRegression.fit!(BayesianLinReg(logX,logLC))
w1,w2=posteriorWeights(lm)
factor=10^w1.val
power=w2.val

lm1=Lens_Model(factor,power,14.)
lm2=Lens_Model(factor,power,28.)

mc_its=Int64(1e6)
base_scale=[144.,5.]
base_min=[10.,0.]
time_scale=90
time_min=4.

cycle_prior=[LogNormal(log(20),log(2)),LogNormal(log(.9),log(1.6))]
end_prior=[Uniform(4,90)]

function compose_priors(phases)
    prs=Vector{Vector{Distribution}}()
    bxes=Vector{Matrix{Float64}}()
    for p in phases
        pr=[repeat(cycle_prior,p)...,repeat(end_prior,p-1)...]
        push!(prs,pr)
        bx=hcat([repeat(base_min,p)...,fill(time_min,p-1)...],
        [repeat(base_scale,p)...,fill(time_scale,p-1)...])
        push!(bxes,GMC_NS.to_unit_ball.(bx,pr))
    end
    return prs,bxes
end

prior,box=compose_priors(2)

d_constants=[X,dpopdist,lm1,mc_its,2]
v_constants=[X,vpopdist,lm1,mc_its,2] 
dv_constants=[X,dvpopdist,lm2,mc_its,2]
constants=[d_constants,v_constants,dv_constants]

uds=Vector{Vector{Function}}([[],[liwi_display],[convergence_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display],[ensemble_display]])

gmc_settings=GMC_DEFAULTS
gmc_settings[end]=25000

for (pth,constants,obs) in zip(paths,constants,obsset)
    if isfile(pth*"/ens")
        e=deserialize(pth*"/ens")
    else
        e=Slice_Ensemble(pth,3000,obs, prior..., constants, box..., gmc_settings)
    end

    converge_ensemble!(e,backup=(true,50),upper_displays=uds, lower_displays=lds, disp_rot_its=100, mc_noise=.3, converge_criterion="compression", converge_factor=1.)
end

ed=deserialize(edpath*"/ens")
ev=deserialize(evpath*"/ens")
edv=deserialize(edvpath*"/ens")

edev=measure_evidence(ed)
evev=measure_evidence(ev)
edvev=measure_evidence(edv)

evidence_ratio=edvev-(edev+evev)

println("\\begin{tabular}{|l|l|l|l|l|}")
println("\\hline")
println("{\\bf Total logZ} & {\\bf Dorsal logZ} & {\\bf Ventral logZ} & {\\bf logZR} & {\\bf \$\\sigma\$ Significance}\\\\ \\hline")
println("\\textbf{$(round(edvev,digits=3))} & $(round(edev,digits=3)) & $(round(evev,digits=3)) & $(round(evidence_ratio,digits=3)) & $(round(evidence_ratio.val/evidence_ratio.err,digits=3))\\\\ \\hline")
println("\\end{tabular}")

mapd=deserialize(ed.models[findmax([m.log_Li for m in ed.models])[2]].path)
mapv=deserialize(ev.models[findmax([m.log_Li for m in ev.models])[2]].path)
mapdv=deserialize(edv.models[findmax([m.log_Li for m in edv.models])[2]].path)

println("\\begin{tabular}{|l|l|l|l|}")
println("\\hline")
println("{\\bf Parameter} & {\\bf Total MAP} & {\\bf Dorsal MAP} & {\\bf Ventral MAP}\\ \\hline")
println("Phase 1 \$CT\$ (h) & $(round(mapdv.θ[1], digits=1)) & $(round(mapd.θ[1], digits=1)) & $(round(mapv.θ[1], digits=1))\\\\ \\hline")
println("Phase 1 \$\\epsilon\$ & $(round(mapdv.θ[2], digits=2)) & $(round(mapd.θ[2], digits=2)) & $(round(mapv.θ[2], digits=2))\\\\ \\hline")
println("Phase 2 \$CT\$ (h) & $(round(mapdv.θ[3], digits=1)) & $(round(mapd.θ[3], digits=1)) & $(round(mapv.θ[3], digits=1))\\\\ \\hline")
println("Phase 2 \$\\epsilon\$ & $(round(mapdv.θ[4], digits=2)) & $(round(mapd.θ[4], digits=2)) & $(round(mapv.θ[4], digits=2))\\\\ \\hline")
println("Transition age & $(round(mapdv.θ[5], digits=1)) & $(round(mapd.θ[5], digits=1)) & $(round(mapv.θ[5], digits=1))\\\\ \\hline")
println("\\end{tabular}")

X=ed.constants[1]

plotticks=[10*i for i in 1:9]

catdobs=vcat([ed.obs[t] for t in 1:length(X)]...)
catvobs=vcat([ev.obs[t] for t in 1:length(X)]...)
catdvobs=vcat([edv.obs[t] for t in 1:length(X)]...)

dXs=vcat([[X[n] for i in 1:length(ed.obs[n])] for n in 1:length(X)]...)
vXs=vcat([[X[n] for i in 1:length(ev.obs[n])] for n in 1:length(X)]...)
dvXs=vcat([[X[n] for i in 1:length(edv.obs[n])] for n in 1:length(X)]...)

d_mean=mapd.disp_mat[:,2]
d_upper=mapd.disp_mat[:,3].-mapd.disp_mat[:,2]
d_lower=mapd.disp_mat[:,2].-mapd.disp_mat[:,1]

v_mean=mapv.disp_mat[:,2]
v_upper=mapv.disp_mat[:,3].-mapv.disp_mat[:,2]
v_lower=mapv.disp_mat[:,2].-mapv.disp_mat[:,1]

dv_mean=mapdv.disp_mat[:,2]
dv_upper=mapdv.disp_mat[:,3].-mapdv.disp_mat[:,2]
dv_lower=mapdv.disp_mat[:,2].-mapdv.disp_mat[:,1]

mapd_plt=scatter(dXs,catdobs, marker=:cross, color=:black, markersize=3, label="Dorsal CMZ population data", showaxis=:y, xticks=plotticks, xformatter=_->"", ylabel="Population")
plot!(mapd_plt, X, d_mean, ribbon=(d_lower,d_upper), color=:green, label="Dorsal model population")
annotate!([(8,150,Plots.text("B",18))])

mapv_plt=scatter(vXs,catvobs, marker=:cross, color=:black, markersize=3, label="Ventral CMZ population data", xlabel="Age (dpf)", ylabel="Population", xticks=plotticks)
plot!(mapv_plt, X, v_mean, ribbon=(v_lower,v_upper), color=:darkmagenta, label="Ventral model population")
annotate!([(8,175,Plots.text("C",18))])

mapdv_plt=scatter(dvXs,catdvobs, marker=:cross, color=:black, markersize=3, label="Total CMZ population data", showaxis=:y, xticks=plotticks, xformatter=_->"", ylabel="Population")
plot!(mapdv_plt, X, dv_mean, ribbon=(dv_lower,dv_upper), color=:darkorange, label="Total model population")
annotate!([(8,300,Plots.text("A",18))])

combined_map=Plots.plot(mapdv_plt,mapd_plt,mapv_plt,layout=grid(3,1), size=(800,900),  link=:x)

savefig(combined_map,"/bench/PhD/Thesis/images/cmz/a10dvMAP.png")

kde=posterior_kde(edv)

ph1marg=marginal(kde,[1;2])
ph2marg=marginal(kde,[3;4])
ph12marg=marginal(kde,[1;3])
transmarg=marginal(kde,[5])

ph1mplt=KernelDensityEstimatePlotting.plot(ph1marg;dimLbls=["Phase 1 CT (hr)","Phase 1 ϵ rate"], axis=[0. 144.; 0. 4.])
ph2mplt=KernelDensityEstimatePlotting.plot(ph2marg;dimLbls=["Phase 2 CT (hr)","Phase 2 ϵ rate"], axis=[0. 144.; 0. 4.])
ph12mplt=KernelDensityEstimatePlotting.plot(ph12marg; dimLbls=["Phase 1 CT (hr)","Phase 2 CT (hr)"],axis=[0. 144.; 0. 144.])
tmplt=KernelDensityEstimatePlotting.plotKDE(transmarg,xlbl="Phase transition age (dpf)", points=false, c=["green"])

combined_marg=Gadfly.gridstack([ph1mplt ph12mplt; ph2mplt tmplt])
img=SVG("/bench/PhD/Thesis/images/cmz/a10dvmarginals.svg",24cm,24cm)
draw(img,combined_marg)
