using CSV,DataFrames,Distributions,NGRefTools,StatsBase,Plots,GMC_NS,Serialization, Measurements, BioSimpleStochastic, Random

Random.seed!(786)

a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"
e2ph="/bench/PhD/NGS_binaries/BSS/A10/e2ph"
e3ph="/bench/PhD/NGS_binaries/BSS/A10/e3ph"
e4ph="/bench/PhD/NGS_binaries/BSS/A10/e4ph"
e5ph="/bench/PhD/NGS_binaries/BSS/A10/e5ph"

a10df=DataFrame(CSV.read(a10pth))
a10df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a10df)]

X=Vector{Float64}()
measure_dict=Dict{String,Vector{Vector{Float64}}}()

measure_dict["PopEst"]=Vector{Vector{Float64}}()
measure_dict["VolEst"]=Vector{Vector{Float64}}()

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for t_df in groupby(a10df, "Time point (d)")
    push!(X,Float64(unique(t_df."Time point (d)")[1]))
    pevec,tvec,dvec,vvec,vevec, rtvec, rplvec, rptvec, ldvec, ondvec, onlvec, oplvec, inlvec, iplvec, gclvec = [zeros(0) for i in 1:15]

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
        length([skipmissing(i_df."Retina thickness")...])>0 && length([skipmissing(i_df."IPL")...])>0 && length([skipmissing(i_df."OPL")...])>0 && length([skipmissing(i_df."RPE length")...])>0 && push!(vevec,(mean(skipmissing(i_df."Retina thickness"))-mean(skipmissing(i_df."OPL"))-mean(skipmissing(i_df."IPL")))*mean(skipmissing(i_df."RPE length"))*(π/4))
    end

    push!(measure_dict["PopEst"],pevec)
    push!(measure_dict["VolEst"],vevec)
end

obs=[(measure_dict["PopEst"][t],measure_dict["VolEst"][t]) for t in 1:length(X)]

mc_its=Int64(1.5e6)
#noise_dist=Normal(0,.01)
noise_dist=[0]
base_scale=[144.,1.]
volconst_scale=25.
time_scale=360

popdist=fit(LogNormal,measure_dict["PopEst"][1])
voldist=fit(LogNormal,measure_dict["VolEst"][1])
ph2_constants=[X, popdist, voldist, noise_dist, mc_its, 2]
ph3_constants=[X, popdist, voldist, noise_dist, mc_its, 3]

ph2_priors=[LogNormal(log(20),log(2)),Beta(2,2),LogNormal(log(20),log(2)),Beta(2,2),Uniform(3,360),LogNormal(log(6.),log(1.6))]
ph2_box=vcat(vcat(fill(base_scale,2)...),time_scale,volconst_scale)
ph2_box=hcat([4.,0.,4.,0.,3.,.1],ph2_box)
ph2_box=GMC_NS.to_unit_ball.(ph2_box,ph2_priors)

ph3_priors=[Normal(15,5),Beta(.76,4),Normal(15,5),Beta(.76,4),Normal(15,5),Beta(.76,4),Normal(360/3,30),Normal(360/3*2,30),Normal(50.,10.)]
ph3_scale=vcat(vcat(fill(base_scale,3)...),fill(time_scale,2)...,volconst_scale)


uds=Vector{Vector{Function}}([[],[liwi_display],[convergence_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display],[ensemble_display]])


for (pth,prior,constants,box) in zip([e2ph, e3ph],[ph2_priors, ph3_priors],[ph2_constants,ph3_constants],[ph2_box, ph3_scale])
    if isfile(pth*"/ens")
        e=deserialize(pth*"/ens")
    else
        e=CMZ_Ensemble(pth,500,obs, prior, constants, box, GMC_DEFAULTS)
    end

    converge_ensemble!(e,backup=(true,50),upper_displays=uds, lower_displays=lds, disp_rot_its=100)
end
