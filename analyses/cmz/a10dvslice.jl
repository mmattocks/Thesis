using BayesianLinearRegression,CSV,DataFrames,Distributions,BioSimpleStochastic,GMC_NS, Serialization

a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"

edpath="/bench/PhD/NGS_binaries/BSS/A10/ed"
evpath="/bench/PhD/NGS_binaries/BSS/A10/ev"
edvpath="/bench/PhD/NGS_binaries/BSS/A10/edv"
paths=[edpath,evpath,edvpath]

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
            length([skipmissing(i_df."CMZ Sum")...])>0 && push!(tvec,mean(skipmissing(i_df."CMZ Sum")))
            length([skipmissing(i_df."Dorsal CMZ (#)")...])>0 && mean(skipmissing(i_df."Dorsal CMZ (#)")) >0. && push!(dvec,mean(skipmissing(i_df."Dorsal CMZ (#)")))
            length([skipmissing(i_df."Ventral CMZ (#)")...])>0 && push!(vvec,mean(skipmissing(i_df."Ventral CMZ (#)")))
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

for (pth,constants,obs) in zip(paths,constants,obsset)
    if isfile(pth*"/ens")
        e=deserialize(pth*"/ens")
    else
        e=Slice_Ensemble(pth,3000,obs, prior..., constants, box..., GMC_DEFAULTS)
    end

    converge_ensemble!(e,backup=(true,50),upper_displays=uds, lower_displays=lds, disp_rot_its=100, mc_noise=.3, converge_criterion="compression", converge_factor=1.)
end
