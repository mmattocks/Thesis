using CSV,DataFrames,Distributions,Measurements,GMC_NS, NGRefTools, Plots,Measures, Serialization
gr()
default(legendfont = 10, guidefont = (12,:bold), tickfont = 10, guide = "x")


pth="/bench/PhD/Thesis/datasets/a37.csv"

a37df=DataFrame(CSV.read(pth))

atg_measure_dict=Dict{String,Vector{Float64}}()
spl_measure_dict=Dict{String,Vector{Float64}}()
ctl_measure_dict=Dict{String,Vector{Float64}}()
uninj_measure_dict=Dict{String,Vector{Float64}}()

msd=[atg_measure_dict,spl_measure_dict,ctl_measure_dict,uninj_measure_dict]

for dict in msd
    dict["DCMZ"]=Vector{Float64}()
    dict["VCMZ"]=Vector{Float64}()
    dict["TCMZ"]=Vector{Float64}()
    dict["NucPop"]=Vector{Float64}()
    dict["Sphericity"]=Vector{Float64}()
    dict["Volume"]=Vector{Float64}()
end

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for (gdf,dict) in zip(groupby(a37df, "Group"),msd)
    for row in eachrow(gdf)
        if row."D/V"=="D"
            push!(dict["DCMZ"],row."PCNA +ve")
            !ismissing(row."Total PCNA") && push!(dict["TCMZ"],row."Total PCNA")
            !ismissing(row."Total nuclei") && push!(dict["NucPop"],row."Total nuclei")
            push!(dict["Sphericity"],row."Sphericity")
            push!(dict["Volume"],row."Volume")
        else
            push!(dict["VCMZ"],row."PCNA +ve")
            dict["Sphericity"][end]+=row."Sphericity";dict["Sphericity"][end]/=2
            dict["Volume"][end]+=row."Volume";dict["Volume"][end]/=2
        end
    end
end


for md in msd
    md["NucPer"]=md["NucPop"]./md["TCMZ"]
end

gmc=GMC_DEFAULTS
gmc[1]=5
gmc[2]=1.e-15
gmcdir="/bench/PhD/NGS_binaries/BSS/A37"


uds=Vector{Vector{Function}}([[convergence_display],[evidence_display],[info_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display]])

n_models=50

max_μ=300.;min_μ=0.
max_λ=1.;min_λ=1e-5
prior=[Uniform(min_μ,max_μ),Uniform(min_λ,max_λ)]
box=[min_μ max_μ;min_λ max_λ]

np_prior=[Uniform(min_μ,5000),Uniform(min_λ,max_λ)]
np_box=[min_μ 5000;min_λ max_λ]

s_max_μ=1.;s_min_μ=0.
s_max_λ=5e5;s_min_λ=1e3
s_prior=[Uniform(s_min_μ,s_max_μ),Uniform(s_min_λ,s_max_λ)]
s_box=[s_min_μ s_max_μ;s_min_λ s_max_λ]

totalmsvec=["TCMZ","NucPop","NucPer","DCMZ","VCMZ","Volume","Sphericity"]
popmsvec=["TCMZ","NucPop","NucPer","DCMZ","VCMZ"]
normmsvec=["Volume","Sphericity"]

for (md,morph) in zip([atg_measure_dict,spl_measure_dict],["ATG","Spl"])
    evdict=Dict{String,Measurement}()

    for ms in totalmsvec
        evdict[ms*"_Combined"]=measurement(0,0)
        evdict[ms*"_Separate"]=measurement(0,0)
    end


    for ms in popmsvec
        ms=="NucPop" ? (pr=np_prior; bx=np_box) : (pr=prior; bx=box)


        c_ens=gmcdir*"/" * morph * "_c_"*ms
        c_obs=vcat(md[ms],ctl_measure_dict[ms])
        if isfile(c_ens*"/ens")
            e=deserialize(c_ens*"/ens")
        else
            e=LogNormal_Ensemble(c_ens,n_models, c_obs, pr, bx, gmc...)
        end

        evdict[ms*"_Combined"]+=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)

        for (msd,p) in zip([md,ctl_measure_dict],["/"*morph*"_","/ctl_"])
            ens=gmcdir*p*ms
            obs=msd[ms]
            if isfile(ens*"/ens")
                e=deserialize(ens*"/ens")
            else
                e=LogNormal_Ensemble(ens,n_models, obs, pr, bx, gmc...)
            end
            evdict[ms*"_Separate"]+=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)
        end
    end

    for ms in normmsvec
        c_ens=gmcdir*"/" * morph * "_c_"*ms
        c_obs=vcat(md[ms],ctl_measure_dict[ms])
        if isfile(c_ens*"/ens")
            e=deserialize(c_ens*"/ens")
        else
            if ms=="Sphericity"
                e=Normal_Ensemble(c_ens,n_models, c_obs, s_prior, s_box, gmc...)
            else
                e=Normal_Ensemble(c_ens,n_models, c_obs, prior, box, gmc...)
            end
        end

        evdict[ms*"_Combined"]+=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)

        for (msd,p) in zip([md,ctl_measure_dict],["/"*morph*"_","/ctl_"])
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

            evdict[ms*"_Separate"]+=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)
        end
    end

    for ms in totalmsvec
        sev=evdict[ms*"_Separate"]
        cev=evdict[ms*"_Combined"]
        ratio=sev-cev
        sig=ratio.val/ratio.err

        println("$morph $ms & $(round(sev,digits=3)) & $(round(cev,digits=3)) & $(round(ratio,digits=3)) & $(round(sig,digits=1)) \\\\ \\hline")
    end
end

colvec=[:darkmagenta,:darkorange,:green]

np_plt=plot_logn_MTDist([atg_measure_dict["NucPer"],spl_measure_dict["NucPer"],ctl_measure_dict["NucPer"]],colvec,[:circle,:diamond,:dtriangle],["ATG Mo.","Spl Mo.","Ctrl Mo."],"Mean Sectional Population per CMZ Cell","Posterior mean lh.")

savefig(np_plt, "/bench/PhD/Thesis/images/rys/morphnuclei.png")

