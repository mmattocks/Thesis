using CSV,DataFrames,Distributions,GMC_NS,Measurements,Serialization

a38pth="/bench/PhD/datasets/A38.csv"
a49pth="/bench/PhD/datasets/A49.csv"

a38df=DataFrame(CSV.read(a38pth))
a38df.Group[1:10].="sib"
a38df.Group[11:22].="rys"
a49df=DataFrame(CSV.read(a49pth))
a49df=a49df[2:end,:]

a38df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a38df)]

sib_measure_dict=Dict{String,Vector{Vector{Float64}}}()
rys_measure_dict=Dict{String,Vector{Vector{Float64}}}()

XCMZ=Vector{Float64}()
XEdU=Vector{Float64}()

for dict in [sib_measure_dict,rys_measure_dict]
    dict["DCMZ"]=Vector{Float64}()
    dict["VCMZ"]=Vector{Float64}()
    dict["TCMZ"]=Vector{Float64}()
    dict["EdU"]=Vector{Float64}()
end

for tdf in groupby(a49df,"Age")
    age=unique(tdf."Age")[1]
    push!(XCMZ,age);push!(XEdU,age)
    xidx=length(XCMZ)
    for sdf in groupby(tdf, "Rys/Sib")
        sr=unique(sdf."Rys/Sib")[1]
        for mdf in groupby(sdf,"D/V/P/WE")
            if unique(mdf."D/V/P/WE")[1] == "D"
                if sr=="sib"
                    push!(sib_measure_dict["DCMZ"],mdf."PCNA +ve")
                    push!(sib_measure_dict["TCMZ"],mdf."PCNA +ve")
                    push!(sib_measure_dict["EdU"],mdf."EdU +ve, PCNA +ve")
                else
                    push!(rys_measure_dict["DCMZ"],mdf."PCNA +ve")
                    push!(rys_measure_dict["TCMZ"],mdf."PCNA +ve")
                    push!(rys_measure_dict["EdU"],mdf."EdU +ve, PCNA +ve")
                end
            elseif unique(mdf."D/V/P/WE")[1] == "V"
                if sr=="sib"
                    push!(sib_measure_dict["VCMZ"],mdf."PCNA +ve")
                    sib_measure_dict["TCMZ"][xidx].+=mdf."PCNA +ve"
                    sib_measure_dict["EdU"][xidx].+=mdf."EdU +ve, PCNA +ve"
                else
                    push!(rys_measure_dict["VCMZ"],mdf."PCNA +ve")
                    rys_measure_dict["TCMZ"][xidx].+=mdf."PCNA +ve"
                    rys_measure_dict["EdU"][xidx].+=mdf."EdU +ve, PCNA +ve"
                end
            end
        end
    end
end

push!(XCMZ,5.)
for sdf in groupby(a38df,"Group")
    xidx=length(XCMZ)
    if unique(sdf."Group")[1]=="sib"
        for mdf in groupby(sdf,"D/V")
            if unique(mdf."D/V")[1] == "D"
                push!(sib_measure_dict["DCMZ"],mdf."PCNA +ve")
                push!(sib_measure_dict["TCMZ"],mdf."PCNA +ve")
            else
                push!(sib_measure_dict["VCMZ"],mdf."PCNA +ve")
                sib_measure_dict["TCMZ"][xidx].+=mdf."PCNA +ve"
            end
        end
    else
        for mdf in groupby(sdf,"D/V")
            if unique(mdf."D/V")[1] == "D"
                push!(rys_measure_dict["DCMZ"],mdf."PCNA +ve")
                push!(rys_measure_dict["TCMZ"],mdf."PCNA +ve")
            else
                push!(rys_measure_dict["VCMZ"],mdf."PCNA +ve")
                rys_measure_dict["TCMZ"][xidx].+=mdf."PCNA +ve"
            end
        end
    end
end

perm=sortperm(XCMZ)
for ms in ["DCMZ","VCMZ","TCMZ"]
    sib_measure_dict[ms]=sib_measure_dict[ms][perm]
    rys_measure_dict[ms]=rys_measure_dict[ms][perm]
end
XCMZ=[4.,10.]

for msd in [sib_measure_dict,rys_measure_dict]
    for (k,v) in msd
        for vec in v
            if any(i->i==0,vec)
                vec[findall(i->i==0,vec)].=1e-3
            end
        end
    end
end

gmc=GMC_DEFAULTS
gmc[1]=5
gmc[2]=1.e-15

uds=Vector{Vector{Function}}([[convergence_display],[evidence_display],[info_display]])
lds=Vector{Vector{Function}}([[model_obs_display],[ensemble_display]])

n_models=50

evdict=Dict{String,Measurement}()

for ms in ["TCMZ","EdU"]
    evdict[ms*"_Combined"]=measurement(0,0)
    evdict[ms*"_Separate"]=measurement(0,0)
end

max_μ=2000.
min_μ=0.
max_λ=500.
min_λ=1e-6

prior=[Uniform(min_μ,max_μ),Uniform(min_λ,max_λ)]
box=[min_μ max_μ;min_λ max_λ]

gmcdir="/bench/PhD/NGS_binaries/BSS/A38/"

for ms in ["TCMZ","EdU"]
    for (n,x) in enumerate(XCMZ)
        c_ens=gmcdir*"/c_"*ms*string(Int64(x))
        c_obs=vcat(sib_measure_dict[ms][n],rys_measure_dict[ms][n])
        if isfile(c_ens*"/ens")
            e=deserialize(c_ens*"/ens")
        else
            e=LogNormal_Ensemble(c_ens,n_models, c_obs, prior, box, gmc...)
        end

        evdict[ms*"_Combined"]+=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)


        for (msd,p) in zip([sib_measure_dict,rys_measure_dict],["/s_","/r_"])
            ens=gmcdir*p*ms*string(Int64(x))
            obs=msd[ms][n]
            if isfile(ens*"/ens")
                e=deserialize(ens*"/ens")
            else
                e=LogNormal_Ensemble(ens,n_models, obs, prior, box, gmc...)
            end
            evdict[ms*"_Separate"]+=converge_ensemble!(e,backup=(true,10000),upper_displays=uds, lower_displays=lds, disp_rot_its=10000, converge_criterion="compression", converge_factor=1e-6)
        end
    end
end

