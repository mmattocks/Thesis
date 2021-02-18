using GMC_NS, CMZNicheSims, ConjugatePriors, CSV, DataFrames, Distributions, Serialization, NGRefTools

a35path="/bench/PhD/Thesis/datasets/A35.csv"
a35df=DataFrame(CSV.File(a35path))

X=Vector{Float64}()
s_popvec=Vector{Float64}()
r_popvec=Vector{Float64}()
sib_y=Vector{Vector{Int64}}()
rys_y=Vector{Vector{Int64}}()

for timedf in groupby(a35df,"Time")
    push!(X,timedf[1,"Time"])
    println(timedf[1,"Time"])
    r_dict=Dict{Tuple,Vector{Int64}}()
    s_dict=Dict{Tuple,Vector{Int64}}()
    for row in eachrow(timedf)
        if !ismissing(row.PCNA) && !ismissing(row.EdU)
            row.Treatment=="Sib" ? (dict=s_dict; y=sib_y) : (dict=r_dict; y=rys_y)
            ind=(row.Block,row.Slide,row.Row)
            if !(ind in keys(dict))
                dict[ind]=[row.PCNA,row.EdU]
            else
                dict[ind][1]+=row.PCNA
                dict[ind][2]+=row.EdU
            end
        end
    end

    rcounts=Vector{Int64}()
    scounts=Vector{Int64}()
    for (dict,ctvec,pvec) in zip([r_dict,s_dict],[rcounts,scounts],[r_popvec,s_popvec])
        for (ind,counts) in dict
            push!(pvec, counts[1])
            push!(ctvec, counts[2])
        end
    end
    
    push!(sib_y,scounts)
    push!(rys_y,rcounts)
end

sib_y=sib_y[sortperm(X)]
rys_y=rys_y[sortperm(X)]
sort!(X)

s_pop_prior=fit(NormalInverseGamma,log.(s_popvec))
r_pop_prior=fit(NormalInverseGamma,log.(r_popvec))

g1_prior=Normal(4.,2.)
tc_prior=NormalInverseGamma(3.5,1.,5.,2.)
s_prior=NormalInverseGamma(4.,2.,1.,1.)
sister_prior=Beta(8.,75.)

s_priors=[marginals(s_pop_prior)..., g1_prior, marginals(tc_prior)..., marginals(s_prior)..., sister_prior]
r_priors=[marginals(r_pop_prior)..., g1_prior, marginals(tc_prior)..., marginals(s_prior)..., sister_prior]

box=[3.3 5.8
        .07 .26
        1. 10.
        .1 10.
        eps() 2.7
        eps() 20.
        .1 200.
        eps() .5]

box=GMC_NS.to_unit_ball.(p1_box,priors_1_pop)

sp="/bench/PhD/NGS_binaries/CNS/A35/sib"
rp="/bench/PhD/NGS_binaries/CNS/A35/rys"
const pulse_time=10.5
const mc_its=Int64(1.5e5)
const end_time=10.5
constants=[X, pulse_time, mc_its, end_time]

ep=sp; priors=s_priors
#ep=rp; priors=r_priors

if isfile(ep*"/ens")
    e=deserialize(ep*"/ens")
else
    @info "Assembling ensemble at $ep"
    gmcd=GMC_DEFAULTS
    gmcd[1]=5
    e=Thymidine_Ensemble(ep, 100, y, priors, constants, box, gmcd)
end

uds=Vector{Vector{Function}}([[liwi_display],[convergence_display],[evidence_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display],[ensemble_display]])

converge_ensemble!(e,backup=(true,1),upper_displays=uds, lower_displays=lds, disp_rot_its=5, mc_noise=.14, converge_factor=1e-3)
