using CSV,DataFrames,Distributions,GMC_NS,Measurements,Serialization

pth="/bench/PhD/datasets/A38.csv"

a38df=DataFrame(CSV.read(pth))
a38df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a38df)]
a38df.Group[1:10].="sib"
a38df.Group[11:22].="rys"

sib_measure_dict=Dict{String,Vector{Float64}}()
rys_measure_dict=Dict{String,Vector{Float64}}()

for dict in [sib_measure_dict,rys_measure_dict]
    dict["DCMZ"]=Vector{Float64}()
    dict["VCMZ"]=Vector{Float64}()
    dict["TCMZ"]=Vector{Float64}()
    dict["NucPop"]=Vector{Float64}()
    dict["Sphericity"]=Vector{Float64}()
    dict["Volume"]=Vector{Float64}()
end

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for (gdf,dict) in zip(groupby(a38df, "Group"),[sib_measure_dict,rys_measure_dict])
    for row in eachrow(gdf)
        if row."D/V"=="D"
            push!(dict["DCMZ"],row."PCNA +ve")
            push!(dict["TCMZ"],row."Total PCNA")
            push!(dict["NucPop"],row."Total nuclei")
            push!(dict["Sphericity"],row."Sphericity")
            push!(dict["Volume"],row."Volume")
        else
            push!(dict["VCMZ"],row."PCNA +ve")
            dict["Sphericity"][end]+=row."Sphericity";dict["Sphericity"][end]/=2
            dict["Volume"][end]+=row."Volume";dict["Volume"][end]/=2
        end
    end
end

for md in [sib_measure_dict,rys_measure_dict]
    md["NucPer"]=md["NucPop"]./md["TCMZ"]
end


gmc=GMC_DEFAULTS
gmc[1]=15
gmc[2]=1.e-15
gmcdir="/bench/PhD/NGS_binaries/GMC_NS/A38"


uds=Vector{Vector{Function}}([[convergence_display],[evidence_display],[info_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display]])

n_models=150

evdict=Dict{String,Measurement}()

totalmsvec=["TCMZ","NucPer","DCMZ","VCMZ","Volume","Sphericity"]
popmsvec=["TCMZ","NucPer","DCMZ","VCMZ"]
normmsvec=["Volume","Sphericity"]

for ms in totalmsvec
    evdict[ms*"_Combined"]=measurement(0,0)
    evdict[ms*"_Separate"]=measurement(0,0)
end

max_μ=300.;min_μ=0.
max_λ=.3;min_λ=0.
prior=[Uniform(min_μ,max_μ),Uniform(min_λ,max_λ)]
box=[min_μ max_μ;min_λ max_λ]

s_max_μ=1.;s_min_μ=0.
s_max_λ=5e5;s_min_λ=0.
s_prior=[Uniform(s_min_μ,s_max_μ),Uniform(s_min_λ,s_max_λ)]
s_box=[s_min_μ s_max_μ;s_min_λ s_max_λ]

for ms in popmsvec
    c_ens=gmcdir*"/c_"*ms
    c_obs=vcat(sib_measure_dict[ms],rys_measure_dict[ms])
    if isfile(c_ens*"/ens")
        e=deserialize(c_ens*"/ens")
    else
        e=LogNormal_Ensemble(c_ens,n_models, c_obs, prior, box, gmc...)
    end

    evdict[ms*"_Combined"]+=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_factor=1e-6)

    for (msd,p) in zip([sib_measure_dict,rys_measure_dict],["/s_","/r_"])
        ens=gmcdir*p*ms
        obs=msd[ms]
        if isfile(ens*"/ens")
            e=deserialize(ens*"/ens")
        else
            e=LogNormal_Ensemble(ens,n_models, obs, prior, box, gmc...)
        end
        evdict[ms*"_Separate"]+=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_factor=1e-6)
    end
end

for ms in normmsvec
    c_ens=gmcdir*"/c_"*ms
    c_obs=vcat(sib_measure_dict[ms],rys_measure_dict[ms])
    if isfile(c_ens*"/ens")
        e=deserialize(c_ens*"/ens")
    else
        if ms=="Sphericity"
            e=Normal_Ensemble(c_ens,n_models, c_obs, s_prior, s_box, gmc...)
        else
            e=Normal_Ensemble(c_ens,n_models, c_obs, prior, box, gmc...)
        end
    end

    evdict[ms*"_Combined"]+=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_factor=1e-6)

    for (msd,p) in zip([sib_measure_dict,rys_measure_dict],["/s_","/r_"])
        ens=gmcdir*p*ms
        obs=msd[ms]
        if isfile(ens*"/ens")
            e=deserialize(ens*"/ens")
        else
            if ms=="Sphericity"
                e=Normal_Ensemble(ens,n_models, obs, s_prior, s_box, gmc...)
            else
                e=Normal_Ensemble(ens,n_models, obs, prior, box, gmc...)
            end
        end

        evdict[ms*"_Separate"]+=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_factor=1e-6)
    end
end

for ms in totalmsvec
    sev=evdict[ms*"_Separate"]
    cev=evdict[ms*"_Combined"]
    ratio=sev-cev
    sig=ratio.val/ratio.err

    println("$ms & $(round(sev,digits=3)) & $(round(cev,digits=3)) & $(round(ratio,digits=3)) & $(round(sig,digits=1)) \\\\ \\hline")
end