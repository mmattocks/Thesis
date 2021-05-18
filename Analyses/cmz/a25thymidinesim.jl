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

pop_dist=fit(LogNormal,popvec)

tc_prior=NormalInverseGamma(3.5,1.,5.,2.)
g1_prior=Beta(2.,2.)
s_prior=Beta(2.,5.)
sister_prior=Beta(8.,75.)

tc_prior2=NormalInverseGamma(3.75,1.,5.,2.)
pop2_prior=Beta(1.1,10.)

priors_1_pop=[marginals(tc_prior)..., g1_prior, s_prior, sister_prior]
priors_2_pop=[marginals(tc_prior2)..., g1_prior, s_prior, sister_prior,pop2_prior]

priors_2_pop=vcat(priors_1_pop,priors_2_pop)

p1_box=[2. 8.
        .1 2
        eps() 1.
        eps() 1.
        eps() .5]

p2_box=vcat(p1_box,vcat(p1_box,[eps() 1.]))

p1_box=GMC_NS.to_unit_ball.(p1_box,priors_1_pop)
p2_box=GMC_NS.to_unit_ball.(p2_box, priors_2_pop)

ep1="/bench/PhD/NGS_binaries/CNS/A25/1pop"
ep2="/bench/PhD/NGS_binaries/CNS/A25/2pop"
const pulse_time=10.5
const mc_its=Int64(1e5)
const end_time=10.5
constants=[pop_dist, X, pulse_time, mc_its, end_time]

for (ep,priors,box) in zip([ep1,ep2],[priors_1_pop,priors_2_pop],[p1_box,p2_box])
    if isfile(ep*"/ens")
        e=deserialize(ep*"/ens")
    else
        @info "Assembling ensemble at $ep"
        gmcd=GMC_DEFAULTS
        gmcd[1]=5
        gmcd[end]=100
        e=Thymidine_Ensemble(ep, 100, y, priors, constants, box, gmcd)
    end

    uds=Vector{Vector{Function}}([[liwi_display],[convergence_display],[ensemble_display]])
    lds=Vector{Vector{Function}}([[model_obs_display],[model_obs_display],[model_obs_display]])

    converge_ensemble!(e,backup=(true,1),upper_displays=uds, lower_displays=lds, disp_rot_its=5, mc_noise=.14, converge_factor=1e-3)
end

e1=deserialize(ep1*"/ens")
e2=deserialize(ep2*"/ens")

ev1=measure_evidence(e1)
ev2=measure_evidence(e2)

CMZNicheSims.print_MAP_output(e1, "/bench/PhD/Thesis/images/cmz/a25MAP.png", its=Int64(1e7))
CMZNicheSims.print_marginals(e1,"/bench/PhD/Thesis/images/cmz/a25marginals.png")