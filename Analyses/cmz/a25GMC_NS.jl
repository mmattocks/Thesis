using GMC_NS, CMZNicheSims, ConjugatePriors, CSV, DataFrames, Distributed, Distributions, NGRefTools, ProgressMeter, Serialization

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

sp_max_mean,sp_max_std=get_lognormal_desc(fit(LogNormal,popvec))
sp_max_var=sp_max_std^2

pop_mu_prior=Uniform(1.,sp_max_mean)
pop_std_prior=Uniform(1.,sp_max_var)
r_prior=Normal(5,1.5)
tc_prior=NormalInverseGamma(3.,1.,5.,2.)
s_prior=NormalInverseGamma(6.,2.5,1.,1.)
sister_prior=Beta(8.,75.)

priors_1_pop=[pop_mu_prior, pop_std_prior, r_prior, marginals(tc_prior)..., marginals(s_prior)..., sister_prior]
priors_2_pop=vcat(priors_1_pop,priors_1_pop)

p1_box=[1. sp_max_mean
        1. sp_max_var
        1. 10.
        1. 4.5
        eps() .5
        eps() .5]

p2_box=vcat(p1_box,p1_box)

prior_sets=[priors_1_pop,priors_2_pop]
ensemble_paths=["/bench/PhD/NGS_binaries/BSS/A25/1pop", "/bench/PhD/NGS_binaries/BSS/A25/2pop"]
const pulse_time=10.5
const mc_its=Int64(5e5)
const end_time=10.5
const retain_run=false
constants=[X, pulse_time, mc_its, end_time, retain_run]


ensembles=Vector{Thymidine_Ensemble}()
for (ps, ep, box) in zip(prior_sets,ensemble_paths,boxes)
    if isfile(ep*"/ens")
        push!(ensembles,deserialize(ep*"/ens"))
    else
        @info "Assembling ensemble at $ep"
        push!(ensembles,Thymidine_Ensemble(ep, 50, y, ps, constants, box, GMC_DEFAULTS))
    end
end

#uds=Vector{Vector{Function}}([[tuning_display],[convergence_display],[tuning_display],[liwi_display],[tuning_display]])
uds=Vector{Vector{Function}}([[tuning_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display]])

evidence=zeros(length(ensembles))
for (ne,e) in enumerate(ensembles)
    evidence[ne]=converge_ensemble!(e, τ0=.03, backup=(true,1),upper_displays=uds, lower_displays=lds, disp_rot_its=5, mc_noise=.139)
end