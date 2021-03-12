using GMC_NS, CMZNicheSims, ConjugatePriors, CSV, DataFrames, Distributions, Serialization, NGRefTools, KernelDensity, Plots, StatsPlots
import BioBackgroundModels: lps
a35path="/bench/PhD/Thesis/datasets/A35.csv"
a35df=DataFrame(CSV.File(a35path))

X=Vector{Float64}()
s_popvec=Vector{Float64}()
r_popvec=Vector{Float64}()
sib_y=Vector{Vector{Int64}}()
rys_y=Vector{Vector{Int64}}()

for timedf in groupby(a35df,"Time")
    push!(X,timedf[1,"Time"])
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

s_box=GMC_NS.to_unit_ball.(box,s_priors)
r_box=GMC_NS.to_unit_ball.(box,s_priors)

sp="/bench/PhD/NGS_binaries/CNS/A35/sib"
rp="/bench/PhD/NGS_binaries/CNS/A35/rys"
const pulse_time=10.5
const mc_its=Int64(1.5e5)
const end_time=10.5
constants=[X, pulse_time, mc_its, end_time]

for (ep, priors, bx, y) in zip([sp,rp], [s_priors, r_priors], [s_box,r_box], [sib_y, rys_y])
    if isfile(ep*"/ens")
    e=deserialize(ep*"/ens")
    else
    @info "Assembling ensemble at $ep"
    gmcd=GMC_DEFAULTS
    gmcd[1]=5
    e=Thymidine_Ensemble(ep, 100, y, priors, constants, bx, gmcd)
    end

    uds=Vector{Vector{Function}}([[liwi_display],[convergence_display],[evidence_display]])
    lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display],[ensemble_display]])

    converge_ensemble!(e,backup=(true,1),upper_displays=uds, lower_displays=lds, disp_rot_its=5, mc_noise=.14, converge_factor=1e-3)
end

se=deserialize(sp*"/ens")
re=deserialize(rp*"/ens")

measure_evidence(se)
measure_evidence(re)

CMZNicheSims.print_MAP_output(se,"/bench/PhD/Thesis/images/rys/a35sibMAP.png", its=Int64(1e7))
CMZNicheSims.print_MAP_output(re,"/bench/PhD/Thesis/images/rys/a35rysMAP.png", its=Int64(1e7))

lt=length(se.priors)
swt=zeros(0); rwt=zeros(0)
sθ=[zeros(0) for param in 1:lt]
rθ=[zeros(0) for param in 1:lt]

for (e,wtvec,θvecs) in zip([se,re],[swt,rwt],[sθ,rθ])
    for (n,rec) in enumerate(e.posterior_samples)
        m=deserialize(rec.path)
        for i in 1:length(θvecs)
            push!(θvecs[i],m.θ[i])
        end
        push!(wtvec,e.log_Liwi[n+1])
    end

    for rec in e.models
        m=deserialize(rec.path)
        for i in 1:length(θvecs)
            push!(θvecs[i],m.θ[i])
        end
        push!(wtvec,(lps(m.log_Li,lps(e.log_Xi[end],-log(length(e.models))))))
    end
end

deleteat!(sθ,3); deleteat!(sθ,4)
deleteat!(rθ,3); deleteat!(rθ,4)

n_pops=1; plotrows=5; param_names=["LogNormal Tc μ", "LogNormal Tc σ²", "G1 Fraction", "S Fraction", "Sister Shift Fraction"]

maps=deserialize(se.models[findmax([m.log_Li for m in se.models])[2]].path)
mapr=deserialize(re.models[findmax([m.log_Li for m in re.models])[2]].path)

plots=Vector{Plots.Plot}()
for (n,(sθvec,rθvec)) in enumerate(zip(sθ,rθ))
    n==3 && (n=4)
    p_label=param_names[n]
    occursin("Fraction",p_label) ? (xls=[0,1]) : (xls=[quantile(se.priors[n],.01),quantile(se.priors[n],.99)])


    sθkde=kde(sθvec,weights=exp.(swt))
    rθkde=kde(rθvec,weights=exp.(rwt))
    n==1 ? (lblpos=:right) : (lblpos=:none)
    plt=StatsPlots.plot(se.priors[n], color=:darkmagenta, fill=true, fillalpha=.5, label="Prior", xlabel=p_label, ylabel="Density", xlims=xls)
    plot!(sθkde.x,sθkde.density, color=:green, fill=true, fillalpha=.5, label="Sib Posterior", legend=lblpos)
    plot!([maps.θ[n],maps.θ[n]],[0,maximum(sθkde.density)],color=:darkgreen, label="sib MAP", legend=lblpos)
    plot!(rθkde.x,rθkde.density, color=:orange, fill=true, fillalpha=.5, label="Rys Posterior", legend=lblpos)
    plot!([mapr.θ[n],mapr.θ[n]],[0,maximum(rθkde.density)],color=:darkred, label="Rys MAP", legend=lblpos)

    push!(plots,plt)
    n_pops>1 && n==5 && push!(plots,plot())
end

combined=plot(plots..., layout=grid(3,n_pops),size=(600*n_pops,800))

savefig(combined, "/bench/PhD/Thesis/images/rys/a35marginals.png")