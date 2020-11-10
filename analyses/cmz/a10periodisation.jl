using CSV,DataFrames,Distributions,NGRefTools,StatsBase,Plots,GMC_NS,Serialization, Measurements, BioSimpleStochastic, Random

Random.seed!(786)

a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"
e2ph="/bench/PhD/NGS_binaries/BSS/A10/e2ph"
e3ph="/bench/PhD/NGS_binaries/BSS/A10/e3ph"
e4ph="/bench/PhD/NGS_binaries/BSS/A10/e4ph"
e5ph="/bench/PhD/NGS_binaries/BSS/A10/e5ph"
e6ph="/bench/PhD/NGS_binaries/BSS/A10/e6ph"
e7ph="/bench/PhD/NGS_binaries/BSS/A10/e7ph"
e10ph="/bench/PhD/NGS_binaries/BSS/A10/e10ph"

paths=[e2ph,e3ph,e4ph]
#,e5ph,e6ph,e7ph,e10ph]

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

            push!(vevec, rthi * (rpel/2)^2 *π)

            push!(spvec,(ov-iv)*(4/5))
        end
    end

    push!(measure_dict["PopEst"],pevec)
    push!(measure_dict["VolEst"],vevec)
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
#noise_dist=Normal(0,.01)
base_scale=[144.,5.]
base_min=[10.,0.]
volconst_scale=5000.
vol_min=10.
time_scale=360
time_min=3.

popdist=fit(LogNormal,measure_dict["PopEst"][1])
voldist=fit(LogNormal,measure_dict["SphEst"][1])

ph2_constants=[X, popdist, voldist, vol_const, mc_its, 2]
ph3_constants=[X, popdist, voldist, vol_const, mc_its, 3]
ph4_constants=[X, popdist, voldist, vol_const, mc_its, 4]
ph5_constants=[X, popdist, voldist, vol_const, mc_its, 5]
ph6_constants=[X, popdist, voldist, vol_const, mc_its, 6]
ph7_constants=[X, popdist, voldist, vol_const, mc_its, 7]
ph10_constants=[X,popdist,voldist,vol_const,mc_its,10]
constants=ph2_constants,ph3_constants,ph4_constants
#,ph5_constants,ph6_constants,ph7_constants,ph10_constants

cycle_prior=[LogNormal(log(20),log(2)),LogNormal(log(.9),log(1.6))]
end_prior=[Uniform(3,360)]
#vol_prior=LogNormal(log(1500),log(2.2))

phases=[2,3,4]#,5,6,7,10]

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

priors,boxes=compose_priors(phases)

uds=Vector{Vector{Function}}([[],[liwi_display],[convergence_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display],[ensemble_display]])

for (pth,prior,constants,box) in zip(paths,priors,constants,boxes)
    if isfile(pth*"/ens")
        e=deserialize(pth*"/ens")
    else
        e=CMZ_Ensemble(pth,3000,obs, prior, constants, box, GMC_DEFAULTS)
    end

    converge_ensemble!(e,backup=(true,50),upper_displays=uds, lower_displays=lds, disp_rot_its=100, mc_noise=.3, converge_factor=1e-3)
end
