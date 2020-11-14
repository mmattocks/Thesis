using GMC_NS, BioSimpleStochastic, ConjugatePriors, CSV, DataFrames, Distributed, Distributions, NGRefTools, ProgressMeter, Serialization

a10path="/bench/PhD/datasets/A10 measurements 2018update.csv"
a10df=DataFrame(CSV.File(a10path))

df3=a10df[a10df."Time point (d)" .== 3, :]

inddict=Dict{Tuple,Vector}()
for row in eachrow(df3)
    ind=(row.Block, row.Slide, row.Row)
    if !(ind in keys(inddict))
        inddict[ind]=[(row["Dorsal CMZ (#)"]+row["Ventral CMZ (#)"])]
    else
        push!(inddict[ind],((row["Dorsal CMZ (#)"]+row["Ventral CMZ (#)"])))
    end
end

popvec=Vector{Float64}()
for (ind, pops) in inddict
    push!(popvec,mean(pops))
end

a25path="/bench/PhD/datasets/a25.csv"
a25df=DataFrame(CSV.File(a25path))

X=Vector{Float64}()
y=Vector{Vector{Int64}}()

for timedf in groupby(a25df,"Time")
    push!(X,timedf[1,"Time"])
    dict=Dict{Tuple,Vector{Int64}}()
    for row in eachrow(timedf)
        ind=(row.Block,row.Slide,row.Row)
        if !(ind in keys(dict))
            dict[ind]=[row.PCNA,row.EdU]
        else
            dict[ind][1]+=row.PCNA
            dict[ind][2]+=row.EdU
        end
    end

    tcounts=Vector{Int64}()
    for (ind,counts) in dict
        global popvec=vcat(popvec, counts[1])
        push!(tcounts, counts[2])
    end
    push!(y,tcounts)
end

y=y[sortperm(X)]
sort!(X)

box_1_pop=[eps() ]

single_pop_NG=fit(NormalInverseGamma,log.(popvec))
r_prior=Normal(5,1.5)
tc_prior=NormalInverseGamma(3.,1.,5.,2.)
s_prior=NormalInverseGamma(5,.2,5.,5.)
sister_prior=Beta(8.,75.)

priors_1_pop=[single_pop_NG, r_prior, tc_prior, s_prior, sister_prior]
# priors_2_pop=[NormalGamma(params(single_pop_NG).*.1...), r_prior, tc_prior, sf_prior,
#             NormalGamma(params(single_pop_NG).*.9...), r_prior, tc_prior, sf_prior]
# priors_3_pop=[NormalGamma(params(single_pop_NG).*.1...), r_prior, tc_prior, sf_prior,
#               NormalGamma(params(single_pop_NG).*.45...), r_prior, tc_prior, sf_prior,
#              NormalGamma(params(single_pop_NG).*.45...), r_prior, tc_prior, sf_prior]
#  priors_4_pop=[NormalGamma(params(single_pop_NG).*.1...), r_prior, tc_prior, sf_prior,
#              NormalGamma(params(single_pop_NG).*.3...), r_prior, tc_prior, sf_prior,
#              NormalGamma(params(single_pop_NG).*.3...), r_prior, tc_prior, sf_prior,
#              NormalGamma(params(single_pop_NG).*.3...), r_prior, tc_prior, sf_prior]

prior_sets=[priors_1_pop]
ensemble_paths=["/bench/PhD/NGS_binaries/BSS/A25/1pop"] 
                # "/bench/PhD/NGS_binaries/BSS/A25/2pop",
                #  "/bench/PhD/NGS_binaries/BSS/A25/3pop",
                #  "/bench/PhD/NGS_binaries/BSS/A25/4pop"]

const pulse_time=10.5
const mc_its=Int64(5e5)
const end_time=10.5
const retain_run=false
constants=[X, pulse_time, mc_its, end_time, retain_run]

function bound_θ!(θ)
    npops=length(θ)/8
    θ[findall(θi->θi<0., θ)].=nextfloat(0.)
    for p in 1:npops
        θ[Int64(8*p)]>1. && (θ[Int64(6*p)]=1.)
    end
    return θ
end

ensembles=Vector{Thymidine_Ensemble}()
for (ps, ep) in zip(prior_sets,ensemble_paths)
    if isfile(ep*"/ens")
        push!(ensembles,deserialize(ep*"/ens"))
    else
        @info "Assembling ensemble at $ep"
        push!(ensembles,Thymidine_Ensemble(ep, 50, y, ps, constants, GMC_DEFAULTS))
    end
end

#uds=Vector{Vector{Function}}([[tuning_display],[convergence_display],[tuning_display],[liwi_display],[tuning_display]])
uds=Vector{Vector{Function}}([[tuning_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display]])

evidence=zeros(length(ensembles))
for (ne,e) in enumerate(ensembles)
    evidence[ne]=converge_ensemble!(e, τ0=.03, backup=(true,1),upper_displays=uds, lower_displays=lds, disp_rot_its=5, mc_noise=.139)
end