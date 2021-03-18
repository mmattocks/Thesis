using BayesianLinearRegression,CSV,DataFrames,Distributions,CMZNicheSims,GMC_NS, Serialization, Plots
a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"

edpth="/bench/PhD/NGS_binaries/CNS/A10/ed_decay"
evpth="/bench/PhD/NGS_binaries/CNS/A10/ev_decay"
empth="/bench/PhD/NGS_binaries/CNS/A10/emulti_decay"
paths=[edpth,evpth,empth]

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

dobs=measure_dict["Dorsal CMZ (#)"]
vobs=measure_dict["Ventral CMZ (#)"]
mobs=[measure_dict["Dorsal CMZ (#)"],measure_dict["Ventral CMZ (#)"]]

obsset=[dobs,vobs,mobs]

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

ed_constants=[X,dpopdist,lm,mc_its]
ev_constants=[X,vpopdist,lm,mc_its]
em_constants=[X,[dpopdist,vpopdist],lm,mc_its]
constants=[ed_constants,ev_constants,em_constants]

priors=[LogNormal(log(20),log(2)),Uniform(eps(),.03),LogNormal(log(.9),log(1.6))]

box=[10. 144.
     eps() .03
     0. 5.]
box=GMC_NS.to_unit_ball.(box,priors)

uds=Vector{Vector{Function}}([[],[liwi_display],[convergence_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display],[ensemble_display]])

gmc_settings=GMC_DEFAULTS
gmc_settings[end]=100

for (pth,constants,obs) in zip(paths,constants,obsset)
    if isfile(pth*"/ens")
        e=deserialize(pth*"/ens")
    else
        if !occursin("emulti",pth)
            e=Decay_Ensemble(pth,1000,obs, priors, constants, box, gmc_settings)
        else
            e=MultiSlice__Decay_Ensemble(pth,1000,obs, priors, constants, box, gmc_settings)
        end
    end

    converge_ensemble!(e,backup=(true,20),upper_displays=uds, lower_displays=lds, disp_rot_its=50, mc_noise=.3, converge_factor=1e-6)
end

ed=deserialize(edpth*"/ens")
ev=deserialize(evpth*"/ens")
em=deserialize(empth*"/ens")

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
dvXs=vcat([[X[n] for i in 1:length(em.obs[1][n])] for n in 1:length(X)]...)

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

mapmd_plt=scatter(dXs,catdobs, marker=:cross, color=:black, markersize=3, label="Dorsal CMZ population data", showaxis=:y, xformatter=_->"", ylabel="Population")
plot!(mapmd_plt, X, md_mean, ribbon=(md_lower,md_upper), color=:green, label="Combined model D population")
annotate!([(8,200,Plots.text("A",18))])

mapmv_plt=scatter(vXs,catvobs, marker=:cross, color=:black, markersize=3, label="Ventral CMZ population data", showaxis=:y, ylabel="Population", xformatter=_->"",)
plot!(mapmv_plt, X, mv_mean, ribbon=(mv_lower,mv_upper), color=:darkmagenta, label="Combined model V population")
annotate!([(8,145,Plots.text("B",18))])

mapd_plt=scatter(dXs,catdobs, marker=:cross, color=:black, markersize=3, label="Dorsal CMZ population data", xlabel="Age (dpf)", xticks=X, ylabel="Population")
plot!(mapd_plt, X, d_mean, ribbon=(d_lower,d_upper), color=:green, label="Separate model D population")
annotate!([(8,150,Plots.text("C",18))])

mapv_plt=scatter(vXs,catvobs, marker=:cross, color=:black, markersize=3, label="Ventral CMZ population data", xlabel="Age (dpf)", ylabel="Population", xticks=X)
plot!(mapv_plt, X, v_mean, ribbon=(v_lower,v_upper), color=:darkmagenta, label="Separate model V population")
annotate!([(8,250,Plots.text("D",18))])

combined_map=Plots.plot(mapmd_plt,mapmv_plt,mapd_plt,mapv_plt,layout=grid(2,2), size=(1200,900),  link=:x)

savefig(combined_map,"/bench/PhD/Thesis/images/cmz/a10dvMAP.png")