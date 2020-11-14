using CSV,DataFrames,Distributions,StatsBase,GMC_NS,Serialization, BioSimpleStochastic, Random, KernelDensityEstimate, KernelDensityEstimatePlotting

Random.seed!(786)

a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"
e2ph="/bench/PhD/NGS_binaries/BSS/A10/e2ph"
e3ph="/bench/PhD/NGS_binaries/BSS/A10/e3ph"

paths=[e2ph,e3ph]

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
constants=ph2_constants,ph3_constants

cycle_prior=[Uniform(10.,144.), Uniform(0.,6)]
end_prior=[Uniform(4.,359.)]

phases=[2,3]
phase_max=[144.,6.]
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

kde2=posterior_kde(e2)
kde3=posterior_kde(e3)

