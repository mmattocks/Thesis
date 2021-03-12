using BayesianLinearRegression,CSV,DataFrames,Distributions,CMZNicheSims,GMC_NS, Serialization, Plots
a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"

edpath="/bench/PhD/NGS_binaries/CNS/A10/ed"
evpath="/bench/PhD/NGS_binaries/CNS/A10/ev"
emultipath="/bench/PhD/NGS_binaries/CNS/A10/emulti"
paths=[edpath,evpath,emultipath]

default(legendfont = (10), guidefont = (12), tickfont = (10))

a10df=DataFrame(CSV.read(a10pth))
a10df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a10df)]

X=Vector{Float64}()
measure_dict=Dict{String,Vector{Vector{Float64}}}()

measure_dict["Dorsal CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["Ventral CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["Lens circumference"]=Vector{Vector{Float64}}()

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for t_df in groupby(a10df, "Time point (d)")
    if !(Float64(unique(t_df."Time point (d)")[1]) in [180., 360.])

        push!(X,Float64(unique(t_df."Time point (d)")[1]))
        dvec,vvec,lcirc = [zeros(0) for i in 1:4]

        for i_df in groupby(t_df,"BSR")
            length([skipmissing(i_df."Lens diameter")...])>0 && push!(lcirc,mean(skipmissing(i_df."Lens diameter"))*π)
        end

        for i_df in groupby(t_df,"BSR")

            length([skipmissing(i_df."Dorsal CMZ (#)")...])>0 && mean(skipmissing(i_df."Dorsal CMZ (#)")) >0. && push!(dvec,mean(skipmissing(i_df."Dorsal CMZ (#)")))
            length([skipmissing(i_df."Ventral CMZ (#)")...])>0 && mean(skipmissing(i_df."Ventral CMZ (#)")) >0. && push!(vvec,mean(skipmissing(i_df."Ventral CMZ (#)")))
        end

        push!(measure_dict["Dorsal CMZ (#)"],dvec)
        push!(measure_dict["Ventral CMZ (#)"],vvec)
        push!(measure_dict["Lens circumference"],lcirc)
    end
end

obsset=[measure_dict["Dorsal CMZ (#)"],measure_dict["Ventral CMZ (#)"],[measure_dict["Dorsal CMZ (#)"],measure_dict["Ventral CMZ (#)"]]]

dpopdist=fit(LogNormal,measure_dict["Dorsal CMZ (#)"][1])
vpopdist=fit(LogNormal,measure_dict["Ventral CMZ (#)"][1])

logLC=vcat([log.(measure_dict["Lens circumference"][n]) for n in 1:length(X)]...)
logX=vcat([[log.(X[n]) for ms in measure_dict["Lens circumference"][n]] for n in 1:length(X)]...)
logX=hcat(ones(length(logX)),logX)

lm=BayesianLinearRegression.fit!(BayesianLinReg(logX,logLC))
w1,w2=posteriorWeights(lm)
factor=10^w1.val
power=w2.val

lm=Lens_Model(factor,power,14.)

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

d_constants=[X,dpopdist,lm,mc_its,2]
v_constants=[X,vpopdist,lm,mc_its,2] 
multi_constants=[X,[dpopdist,vpopdist],lm,mc_its,2]
constants=[d_constants,v_constants,multi_constants]

uds=Vector{Vector{Function}}([[],[liwi_display],[convergence_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display],[ensemble_display]])

gmc_settings=GMC_DEFAULTS
gmc_settings[end]=25000

for (pth,constants,obs) in zip(paths,constants,obsset)
    if isfile(pth*"/ens")
        e=deserialize(pth*"/ens")
    else
        if !occursin("emulti",pth)
            e=Slice_Ensemble(pth,3000,obs, prior..., constants, box..., gmc_settings)
        else
            e=MultiSlice_Ensemble(pth,3000,obs, prior..., constants, box..., gmc_settings)
        end
    end

    converge_ensemble!(e,backup=(true,20),upper_displays=uds, lower_displays=lds, disp_rot_its=50, mc_noise=.3, converge_factor=1e-6)
end

ed=deserialize(edpath*"/ens")
ev=deserialize(evpath*"/ens")
em=deserialize(emultipath*"/ens")

edev=measure_evidence(ed)
evev=measure_evidence(ev)
emev=measure_evidence(em)

evidence_ratio=emev-(edev+evev)

println("\\begin{tabular}{|l|l|l|l|l|}")
println("\\hline")
println("{\\bf Combined D/V model logZ} & {\\bf Separate D/V model logZ} & {\\bf logZR} & {\\bf \$\\sigma\$ Significance}\\\\ \\hline")
println("\\textbf{$(round(emev,digits=3))} & $(round(edev+evev,digits=3)) & $(round(evidence_ratio,digits=3)) & $(round(evidence_ratio.val/evidence_ratio.err,digits=3))\\\\ \\hline")
println("\\end{tabular}")

mapd=deserialize(ed.models[findmax([m.log_Li for m in ed.models])[2]].path)
mapv=deserialize(ev.models[findmax([m.log_Li for m in ev.models])[2]].path)
mapm=deserialize(em.models[findmax([m.log_Li for m in em.models])[2]].path)

println("\\begin{tabular}{|l|l|l|l|}")
println("\\hline")
println("{\\bf Parameter} & {\\bf Total MAP} & {\\bf Dorsal MAP} & {\\bf Ventral MAP}\\ \\hline")
println("Phase 1 \$CT\$ (h) & $(round(mapm.θ[1], digits=1)) & $(round(mapd.θ[1], digits=1)) & $(round(mapv.θ[1], digits=1))\\\\ \\hline")
println("Phase 1 \$\\epsilon\$ & $(round(mapm.θ[2], digits=2)) & $(round(mapd.θ[2], digits=2)) & $(round(mapv.θ[2], digits=2))\\\\ \\hline")
println("Phase 2 \$CT\$ (h) & $(round(mapm.θ[3], digits=1)) & $(round(mapd.θ[3], digits=1)) & $(round(mapv.θ[3], digits=1))\\\\ \\hline")
println("Phase 2 \$\\epsilon\$ & $(round(mapm.θ[4], digits=2)) & $(round(mapd.θ[4], digits=2)) & $(round(mapv.θ[4], digits=2))\\\\ \\hline")
println("Transition age & $(round(mapm.θ[5], digits=1)) & $(round(mapd.θ[5], digits=1)) & $(round(mapv.θ[5], digits=1))\\\\ \\hline")
println("\\end{tabular}")

X=ed.constants[1]

plotticks=[10*i for i in 1:9]

catdobs=vcat([ed.obs[t] for t in 1:length(X)]...)
catvobs=vcat([ev.obs[t] for t in 1:length(X)]...)

dXs=vcat([[X[n] for i in 1:length(ed.obs[n])] for n in 1:length(X)]...)
vXs=vcat([[X[n] for i in 1:length(ev.obs[n])] for n in 1:length(X)]...)
dvXs=vcat([[X[n] for i in 1:length(em.obs[n])] for n in 1:length(X)]...)

d_mean=mapd.disp_mat[:,2]
d_upper=mapd.disp_mat[:,3].-mapd.disp_mat[:,2]
d_lower=mapd.disp_mat[:,2].-mapd.disp_mat[:,1]

v_mean=mapv.disp_mat[:,2]
v_upper=mapv.disp_mat[:,3].-mapv.disp_mat[:,2]
v_lower=mapv.disp_mat[:,2].-mapv.disp_mat[:,1]

md_mean=mapm.slices[1].disp_mat[:,2]
md_upper=mapm.slices[1].disp_mat[:,3].-mapm.slices[1].disp_mat[:,2]
md_lower=mapm.slices[1].disp_mat[:,2].-mapm.slices[1].disp_mat[:,1]

mv_mean=mapm.slices[2].disp_mat[:,2]
mv_upper=mapm.slices[2].disp_mat[:,3].-mapm.slices[2].disp_mat[:,2]
mv_lower=mapm.slices[2].disp_mat[:,2].-mapm.slices[2].disp_mat[:,1]

mapmd_plt=scatter(dXs,catdobs, marker=:cross, color=:black, markersize=3, label="Dorsal CMZ population data", showaxis=:y, xticks=plotticks, xformatter=_->"", ylabel="Population")
plot!(mapmd_plt, X, md_mean, ribbon=(md_lower,md_upper), color=:green, label="Combined model D population")
annotate!([(8,200,Plots.text("A",18))])

mapmv_plt=scatter(vXs,catvobs, marker=:cross, color=:black, markersize=3, label="Ventral CMZ population data", showaxis=:y, ylabel="Population", xticks=plotticks, xformatter=_->"",)
plot!(mapmv_plt, X, mv_mean, ribbon=(mv_lower,mv_upper), color=:darkmagenta, label="Combined model V population")
annotate!([(8,145,Plots.text("B",18))])

mapd_plt=scatter(dXs,catdobs, marker=:cross, color=:black, markersize=3, label="Dorsal CMZ population data", xlabel="Age (dpf)", xticks=plotticks, ylabel="Population")
plot!(mapd_plt, X, d_mean, ribbon=(d_lower,d_upper), color=:green, label="Separate model D population")
annotate!([(8,150,Plots.text("C",18))])

mapv_plt=scatter(vXs,catvobs, marker=:cross, color=:black, markersize=3, label="Ventral CMZ population data", xlabel="Age (dpf)", ylabel="Population", xticks=plotticks)
plot!(mapv_plt, X, v_mean, ribbon=(v_lower,v_upper), color=:darkmagenta, label="Separate model V population")
annotate!([(8,200,Plots.text("D",18))])

combined_map=Plots.plot(mapmd_plt,mapmv_plt,mapd_plt,mapv_plt,layout=grid(2,2), size=(1200,900),  link=:x)

savefig(combined_map,"/bench/PhD/Thesis/images/cmz/a10dvMAP.png")

###KDE FIG
kdes=posterior_kde(em,bivar=[1=>2,3=>4])

balx=10.:1.:144.
baly=[2^(24/ct)-1 for ct in balx]

ph1marg=contour(kdes[6].x,kdes[6].y,transpose(kdes[6].density), c=:viridis, xlabel="Phase 1 CT (hr)", ylabel="Phase 1 ϵ rate", colorbar=:false, xlims=[0,144], ylims=[0,5], levels=200, clim=(0,.5), legend=(.8,.3), xticks=24.:24.:144.)
plot!(balx,baly,color=:black,linestyle=:dash,label="Balanced pop.")
scatter!([mapm.θ[1]],[mapm.θ[2]], marker=:+, markersize=15, markerstrokestyle=:bold, color=:magenta, label="MAP")
for i in 1:5
    scatter!([mapm.θ[1]],[mapm.θ[2]], marker=:+, markersize=15, markerstrokestyle=:bold, color=:magenta, label=:none, framestyle=:origin)
end
lens!([15, 45], [.35, 1.6], inset = (1, bbox(0.2, 0.05, 0.75, 0.55)), colorbar=:false, framestyle=:box, subplot=2, xformatter=_->"", yformatter=_->"")

annotate!([(4,4.5,Plots.text("A",18))])

ph1ctmarg=plot(Uniform(10,144), color=:darkmagenta, fill=true, 
fillalpha=.5,xlims=[0,144], ylabel="p",xaxis=false, xformatter=_->"", legend=:none)
plot!(kdes[1].x,kdes[1].density,color=:green, fill=true, fillalpha=.5)
plot!([mapm.θ[1],mapm.θ[1]],[0,maximum(kdes[1].density)],color=:black,linewidth=1)
yflip!(ph1ctmarg)

ymarg=map(x->pdf(cycle_prior[2],x),kdes[2].x)
ph1ermarg=plot(ymarg,kdes[2].x, color=:darkmagenta, fill=true, fillalpha=.5, ylims=[0,5], legend=:none, yaxis=false, yformatter=_->"", xlabel="p")
plot!(kdes[2].density,kdes[2].x,color=:green, fill=true, fillalpha=.5)
plot!([0,maximum(kdes[2].density)],[mapm.θ[2],mapm.θ[2]],color=:black,linewidth=1)
xflip!(ph1ermarg)

marglayout=@layout [ymarg{.1w} bivar{.8h}; empty xmarg]

cbar=heatmap(transpose(ones(101,1).*(0:0.005:.5)), c=:viridis, legend=:none, yticks=:none, xticks=(1:20:101, string.(0:0.1:.5)),xlabel="Probability density")

empty=plot([],fill=true,xaxis=false,yaxis=false,grid=false, color=:darkmagenta, xformatter=_->"", yformatter=_->"", label="Prior")
plot!([],color=:green,fill=true,label="Posterior")
plot!([], color=:black, linewidth=2, label="MAP")
ph1c=plot(ph1ermarg,ph1marg,empty,ph1ctmarg,layout=marglayout, link=:both)

ph2marg=contour(kdes[7].x,kdes[7].y,transpose(kdes[7].density), c=:viridis, xlabel="Phase 2 CT (hr)", ylabel="Phase 2 ϵ rate", colorbar=:false, xlims=[0,144], ylims=[0,5], levels=200, clim=(0,.5), legend=(.8,.3), xticks=24.:24.:144.)
plot!(balx,baly,color=:black,linestyle=:dash,label="Balanced pop.")
scatter!([mapm.θ[3]],[mapm.θ[4]], marker=:+, markersize=15, color=:magenta, label="MAP")
for i in 1:5
    scatter!([mapm.θ[3]],[mapm.θ[4]], marker=:+, markersize=15, color=:magenta, label=:none, framestyle=:origin)
end
lens!([15, 45], [.35, 1.6], inset = (1, bbox(0.2, 0.05, 0.75, 0.55)), colorbar=:false, framestyle=:box, subplot=2, xformatter=_->"", yformatter=_->"")

annotate!(ph2marg,[(4,4.5,Plots.text("B",18))])


ph2ctmarg=plot(Uniform(10,144), color=:darkmagenta, fill=true, 
fillalpha=.5, ylabel="p", xlims=[0,144], xaxis=false, xformatter=_->"", legend=:none)
plot!(kdes[3].x,kdes[3].density,color=:green, fill=true, fillalpha=.5)
plot!([mapm.θ[3],mapm.θ[3]],[0,maximum(kdes[3].density)],color=:black,linewidth=1)
yflip!(ph2ctmarg)

ymarg=map(x->pdf(cycle_prior[2],x),kdes[4].x)
ph2ermarg=plot(ymarg,kdes[4].x, color=:darkmagenta, fill=true, fillalpha=.5, ylims=[0,5], legend=:none, yaxis=false, yformatter=_->"", xlabel="p")
plot!(kdes[4].density,kdes[4].x,color=:green, fill=true, fillalpha=.5)
plot!([0,maximum(kdes[4].density)],[mapm.θ[4],mapm.θ[4]],color=:black,linewidth=1)
xflip!(ph2ermarg)

marglayout2=@layout [ymarg{.1w} bivar{.8h}; empty xmarg]

ph2c=plot(ph2ermarg,ph2marg,deepcopy(empty),ph2ctmarg,layout=marglayout2, link=:both)

pend=plot(end_prior[1], color=:darkmagenta, fill=true, fillalpha=.5, xlims=[0,90], ylabel="Density", xlabel="Phase transition age", xticks=X,  label="Prior", legend=:none)
plot!(kdes[5].x, kdes[5].density, color=:green, fill=true, fillalpha=.5, label="Posterior")
plot!([mapm.θ[5],mapm.θ[5]],[0,maximum(kdes[5].density)],color=:black, label="MAP")
annotate!([(1,.015,Plots.text("C",18))])

combined_layout=@layout[ph1{.45h}
                        cbar{.01h}
                        ph2{.45h}
                        pend{.09h}]

combined_marg=plot(ph1c,cbar,ph2c,pend,layout=combined_layout,size=(1000,1200))

savefig(combined_marg,"/bench/PhD/Thesis/images/cmz/a10dvmarginals.png")