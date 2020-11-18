using CSV,DataFrames,Distributions,StatsBase,GMC_NS,Serialization, BioSimpleStochastic, Random, KernelDensityEstimate, KernelDensityEstimatePlotting, Plots
gr()
default(legendfont = (8,"courier"), guidefont = (10,"courier"), tickfont = (8,"courier"))

Random.seed!(786)

a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"
e2ph="/bench/PhD/NGS_binaries/BSS/A10/e2ph"
e3ph="/bench/PhD/NGS_binaries/BSS/A10/e3ph"

paths=[e3ph]

a10df=DataFrame(CSV.read(a10pth))
a10df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a10df)]

X=Vector{Float64}()
measure_dict=Dict{String,Vector{Vector{Float64}}}()

measure_dict["PopEst"]=Vector{Vector{Float64}}()
measure_dict["VolEst"]=Vector{Vector{Float64}}()
measure_dict["SphEst"]=Vector{Vector{Float64}}()
measure_dict["CircEst"]=Vector{Vector{Float64}}()
measure_dict["ThiEst"]=Vector{Vector{Float64}}()

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for t_df in groupby(a10df, "Time point (d)")
    push!(X,Float64(unique(t_df."Time point (d)")[1]))
    pevec,vevec,tvec,ldvec,spvec,cvec,thivec = [zeros(0) for i in 1:7]

    for i_df in groupby(t_df,"BSR")
        length([skipmissing(i_df."Lens diameter")...])>0 && push!(ldvec,mean(skipmissing(i_df."Lens diameter")))
    end

    for i_df in groupby(t_df,"BSR")
        if length([skipmissing(i_df."CMZ Sum")...])>0 && mean(skipmissing(i_df."CMZ Sum")) > 0
            push!(tvec,mean(skipmissing(i_df."CMZ Sum")))
            if length([skipmissing(i_df."Lens diameter")...])>0
                push!(pevec,mean(skipmissing(i_df."CMZ Sum"))*mean(skipmissing(i_df."Lens diameter"))/14.)
            else
                push!(pevec,mean(skipmissing(i_df."CMZ Sum"))*mean(ldvec)/14.)
            end
        end
        if length([skipmissing(i_df."Retina thickness")...])>0 &&
         length([skipmissing(i_df."IPL")...])>0 &&
         length([skipmissing(i_df."OPL")...])>0 &&
         length([skipmissing(i_df."RPE length")...])>0

            rthi=mean(skipmissing(i_df."Retina thickness"))-mean(skipmissing(i_df."OPL"))-mean(skipmissing(i_df."IPL"))
            push!(thivec,rthi)
            rpel=mean(skipmissing(i_df."RPE length"))
            if length([skipmissing(i_df."Lens diameter")...])>0
                rcirc=rpel+mean(skipmissing(i_df."Lens diameter"))
            else
                rcirc=rpel+mean(ldvec)
            end

            push!(cvec,rcirc)

            or=rcirc/2π
            ir=or-.5*rthi

            ov=(4/3)*π*(or^3)
            iv=(4/3)*π*(ir^3)

            push!(spvec,(ov-iv)*(4/5))
        end
    end

    push!(measure_dict["PopEst"],pevec)
    push!(measure_dict["SphEst"],spvec)
    push!(measure_dict["CircEst"],cvec)
    push!(measure_dict["ThiEst"],thivec)
end

obs=[(measure_dict["PopEst"][t],measure_dict["SphEst"][t]) for t in 1:length(X)]

nucpop3d=2075.6 #from a38 sib_measure_dict NucPop mean
or=mean(measure_dict["CircEst"][1])/2π
ir=or-.5*mean(measure_dict["ThiEst"][1])
sxvol_3d=(((π*or^2)-(π*ir^2))*14)*(4/5)
vol_const=sxvol_3d/nucpop3d
mc_its=Int64(1.e6)

popdist=fit(LogNormal,measure_dict["PopEst"][1])
voldist=fit(LogNormal,measure_dict["SphEst"][1])

ph2_constants=[X, popdist, voldist, vol_const, mc_its, 2]
ph3_constants=[X, popdist, voldist, vol_const, mc_its, 3]
constants=[ph3_constants]

cycle_prior=[Uniform(10.,144.), LogNormal(log(.9),log(2))]
end_prior=[Uniform(4.,359.)]

phase_max=[144.,5.]
phase_min=[10.,0.]
time_max=359
time_min=4.

function compose_priors(phases)
    prs=Vector{Vector{Distribution}}()
    bxes=Vector{Matrix{Float64}}()
    for p in phases
        pr=[repeat(cycle_prior,p)...,repeat(end_prior,p-1)...]
        push!(prs,pr)
        bx=hcat([repeat(phase_min,p)...,fill(time_min,p-1)...],
        [repeat(phase_max,p)...,fill(time_max,p-1)...])
        push!(bxes,GMC_NS.to_unit_ball.(bx,pr))
    end
    return prs,bxes
end

phases=[3]
priors,boxes=compose_priors(phases)

uds=Vector{Vector{Function}}([[],[liwi_display],[convergence_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display],[ensemble_display]])

for (pth,prior,constants,box) in zip(paths,priors,constants,boxes)
    if isfile(pth*"/ens")
        e=deserialize(pth*"/ens")
    else
        e=CMZ_Ensemble(pth,3000,obs, prior, constants, box, GMC_DEFAULTS)
    end

    converge_ensemble!(e,backup=(true,50),upper_displays=uds, lower_displays=lds, disp_rot_its=100, mc_noise=.3, converge_criterion="compression", converge_factor=1.)
end

e2=deserialize(e2ph*"/ens")
e3=deserialize(e3ph*"/ens")

map2=deserialize(e2.models[findmax([m.log_Li for m in e2.models])[2]].path)
map3=deserialize(e3.models[findmax([m.log_Li for m in e3.models])[2]].path)

X=e2.constants[1]
catpobs=vcat([e2.obs[t][1] for t in 1:length(X)]...)
catvobs=vcat([e2.obs[t][2] for t in 1:length(X)]...)
Xs=vcat([[X[n] for i in 1:length(e2.obs[n][1])] for n in 1:length(X)]...)

p2_mean=map2.disp_mat[:,2,1]
p2_upper=map2.disp_mat[:,1,1].-map2.disp_mat[:,2,1]
p2_lower=map2.disp_mat[:,2,1].-map2.disp_mat[:,3,1]

map2_popplt=scatter(Xs,catpobs, marker=:cross, color=:black, markersize=3, label="CMZ population estimate", ylabel="Population", xaxis=nothing, showaxis=:y)
plot!(map2_popplt, X, p2_mean, ribbon=(p2_lower,p2_upper), color=:magenta, label="2-phase model population")
annotate!([(30,1.25e4,Plots.text("A1",18))])

v2_mean=map2.disp_mat[:,2,2]
v2_upper=map2.disp_mat[:,1,2].-map2.disp_mat[:,2,2]
v2_lower=map2.disp_mat[:,2,2].-map2.disp_mat[:,3,2]

map2_volplt=scatter(Xs,catvobs, marker=:cross, color=:black, markersize=3, label="Retinal volume estimate", ylabel="Volume (μm^3)", xlabel="Age (dpf)", legend=:top)
plot!(map2_volplt, X, v2_mean, ribbon=(v2_lower,v2_upper), color=:green, label="2-phase model volume")
annotate!([(30,3.5e8,Plots.text("A2",18))])


p3_mean=map3.disp_mat[:,2,1]
p3_upper=map3.disp_mat[:,1,1].-map3.disp_mat[:,2,1]
p3_lower=map3.disp_mat[:,2,1].-map3.disp_mat[:,3,1]

map3_popplt=scatter(Xs,catpobs, marker=:cross, color=:black, markersize=3, label="CMZ population estimate", ylabel="Population", xaxis=nothing, showaxis=:y)
plot!(map3_popplt, X, p3_mean, ribbon=(p3_lower,p3_upper), color=:magenta, label="3-phase model population")
annotate!([(30,1.5e4,Plots.text("B1",18))])

v3_mean=map3.disp_mat[:,2,2]
v3_upper=map3.disp_mat[:,1,2].-map3.disp_mat[:,2,2]
v3_lower=map3.disp_mat[:,2,2].-map3.disp_mat[:,3,2]

map3_volplt=scatter(Xs,catvobs, marker=:cross, color=:black, markersize=3, label="Retinal volume estimate", ylabel="Volume (μm^3)", xlabel="Age (dpf)", legend=:top)
plot!(map3_volplt, X, v3_mean, ribbon=(v3_lower,v3_upper), color=:green, label="3-phase model volume")
annotate!([(30,3.5e8,Plots.text("B2",18))])

combined_map=Plots.plot(map2_popplt,map3_popplt,map2_volplt,map3_volplt,layout=grid(2,2), size=(900,600))

savefig(combined_map,"/bench/PhD/Thesis/images/cmz/a10pMAP.png")

println("MAP 2ph θ: $(map2.θ)")
println("MAP 3ph θ: $(map3.θ)")

kde2=posterior_kde(e2)
# kde3=posterior_kde(e3)

ph1marg=marginal(kde2,[1;2])
ph2marg=marginal(kde2,[3;4])
ph12marg=marginal(kde2,[1;3])
transmarg=marginal(kde2,[5])

ph1mplt=KernelDensityEstimatePlotting.plot(ph1marg;dimLbls=["RPC cycle length (hr)","CMZ Exit rate"], axis=[0. 144.; 0. 4.])
ph2mplt=KernelDensityEstimatePlotting.plot(ph2marg;dimLbls=["RPC cycle length (hr)","CMZ Exit rate"], axis=[0. 144.; 0. 4.])
ph12mplt=KernelDensityEstimatePlotting.plot(ph12marg; dimLbls=["Phase 1 cycle length (hr)","Phase 2 cycle length (hr)"],axis=[0. 144.; 0. 144.])
tmplt=KernelDensityEstimatePlotting.plotKDE(transmarg,xlbl="Phase transition age (dpf)", points=false, c=["green"])

combined_marg=gridstack([ph1mplt ph12mplt; ph2mplt tmplt])
img=SVG("/bench/PhD/Thesis/images/cmz/a10pmarginals.svg",24cm,24cm)
draw(img,combined_marg)