using GMC_NS, CMZNicheSims, ConjugatePriors, CSV, DataFrames, Distributions, Serialization, NGRefTools

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

pop_prior=fit(NormalInverseGamma,log.(popvec))

r_prior=Normal(5,1.5)
tc_prior=NormalInverseGamma(3.,1.,5.,2.)
s_prior=NormalInverseGamma(6.,2.5,1.,1.)
sister_prior=Beta(8.,75.)

priors_1_pop=[marginals(pop_prior)..., r_prior, marginals(tc_prior)..., marginals(s_prior)..., sister_prior]

p1_box=[3.3 5.8
        .07 .26
        1. 10.
        .1 10.
        eps() 2.7
        4. 20.
        .1 200.
        eps() .5]

p1_box=GMC_NS.to_unit_ball.(p1_box,priors_1_pop)

ep="/bench/PhD/NGS_binaries/CNS/A25/1pop"
const pulse_time=10.5
const mc_its=Int64(5e5)
const end_time=10.5
constants=[X, pulse_time, mc_its, end_time]

if isfile(ep*"/ens")
    e=deserialize(ep*"/ens")
else
    @info "Assembling ensemble at $ep"
    gmcd=GMC_DEFAULTS
    gmcd[1]=5
    e=Thymidine_Ensemble(ep, 100, y, priors_1_pop, constants, p1_box, gmcd)
end

uds=Vector{Vector{Function}}([[liwi_display],[convergence_display],[evidence_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display],[ensemble_display]])

converge_ensemble!(e,backup=(true,1),upper_displays=uds, lower_displays=lds, disp_rot_its=5, mc_noise=.14, converge_factor=1e-3)

