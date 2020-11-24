using CSV, DataFrames, Measurements, Plots, Distributions
import Plots:Plot

rand6pth="/bench/PhD/Thesis/datasets/6dpfRys-Histones random primers qPCR.csv"
rand8pth="/bench/PhD/Thesis/datasets/8dpfRys-Histones random primers qPCR.csv"
dt6pth="/bench/PhD/Thesis/datasets/6dpfRys-Histones oligo dT qPCR.csv"
dt8pth="/bench/PhD/Thesis/datasets/8dpfRys-Histones oligo dT qPCR.csv"

rand6=DataFrame(CSV.read(rand6pth))
rand8=DataFrame(CSV.read(rand8pth))
dt6=DataFrame(CSV.read(dt6pth))
dt8=DataFrame(CSV.read(dt8pth))

randdfs=[rand6,rand8]
dtdfs=[dt6,dt8]
samples=Dict{Tuple{String,Int64},Vector{String}}()

measure_dict=Dict{Tuple{String,Int64,Int64,String,String},Vector{Float64}}()

for (df,d) in zip(randdfs,[6,8])
    samples[("rand",d)]=unique(df.Sample)
    for row in eachrow(df)
        str=("rand",d,row.Plate,row.Sample,row.Detector)
        if !ismissing(row.Ct)
            str in keys(measure_dict) ? push!(measure_dict[str],row.Ct) : measure_dict[str]=[row.Ct]
        end
    end
end

for (df,d) in zip(dtdfs,[6,8])
    samples[("dt",d)]=unique(df."Sample Name")
    for row in eachrow(df)
        str=("dt",d,row.Plate,row."Sample Name",row."Target Name")
        if !ismissing(row.CT)
            str in keys(measure_dict) ? push!(measure_dict[str],row.CT) : measure_dict[str]=[row.CT]
        end
    end
end

X=[6.,8.]
rand_sib=Dict{String,Vector{Vector{Measurement}}}()
rand_rys=Dict{String,Vector{Vector{Measurement}}}()
dt_sib=Dict{String,Vector{Vector{Measurement}}}()
dt_rys=Dict{String,Vector{Vector{Measurement}}}()


for ms in ["H2A","H2B","H3","H4"]
    for d in [6,8]
        xidx=findfirst(isequal(d),X)
        for plate in [1,2,3]
            wt_hk=measure_dict[("rand",d,plate,"WT cDNA 10x","ActinB")]
            length(wt_hk)>1 ? wt_hk_m=measurement(mean(wt_hk),std(wt_hk)) : wt_hk_m=measurement(mean(wt_hk),0.)
            wt_trg=measure_dict[("rand",d,plate,"WT cDNA 10x",ms)]
            length(wt_hk)>1 ? wt_trg_m=measurement(mean(wt_trg),std(wt_trg)) : wt_trg_m=measurement(mean(wt_trg),0.)
            wt_ΔCt=wt_trg_m-wt_hk_m

            for (smpl, msdict) in zip(["Sib cDNA 10x","Rys cDNA 10x"],[rand_sib,rand_rys])
                smpl_hk=measure_dict[("rand",d,plate,smpl,"ActinB")]
                length(smpl_hk)>1 ? smpl_hk_m=measurement(mean(smpl_hk),std(smpl_hk)) : smpl_hk_m=measurement(mean(smpl_hk),0.)
                smpl_trg=measure_dict[("rand",d,plate,smpl,ms)]
                length(wt_hk)>1 ? smpl_trg_m=measurement(mean(smpl_trg),std(smpl_trg)) : smpl_trg_m=measurement(mean(smpl_trg),0.)
                smpl_ΔCt=smpl_trg_m-smpl_hk_m
                ΔΔCt=smpl_ΔCt-wt_ΔCt
                fold_change=2^-ΔΔCt
                ms in keys(msdict) ? push!(msdict[ms][xidx],fold_change) : msdict[ms]=[[fold_change],Vector{Float64}()]
            end
        end
    end
end

msvec=["H2A","H2B","H3","H4"]

for ms in msvec
    for d in [6,8]
        xidx=findfirst(isequal(d),X)
        for plate in [1,2,3]
            wt_hk=measure_dict[("dt",d,plate,"WT cDNA","ActinB")]
            length(wt_hk)>1 ? wt_hk_m=measurement(mean(wt_hk),std(wt_hk)) : wt_hk_m=measurement(mean(wt_hk),0.)
            wt_trg=measure_dict[("dt",d,plate,"WT cDNA",ms)]
            length(wt_trg)>1 ? wt_trg_m=measurement(mean(wt_trg),std(wt_trg)) : wt_trg_m=measurement(mean(wt_trg),0.)
            wt_ΔCt=wt_trg_m-wt_hk_m

            for (smpl, msdict) in zip(["Sib cDNA","Rys cDNA"],[dt_sib,dt_rys])
                smpl_hk=measure_dict[("dt",d,plate,smpl,"ActinB")]
                length(smpl_hk)>1 ? smpl_hk_m=measurement(mean(smpl_hk),std(smpl_hk)) : smpl_hk_m=measurement(mean(smpl_hk),0.)
                smpl_trg=measure_dict[("dt",d,plate,smpl,ms)]
                length(smpl_trg)>1 ? smpl_trg_m=measurement(mean(smpl_trg),std(smpl_trg)) : smpl_trg_m=measurement(mean(smpl_trg),0.)
                smpl_ΔCt=smpl_trg_m-smpl_hk_m
                ΔΔCt=smpl_ΔCt-wt_ΔCt
                fold_change=2^-ΔΔCt
                ms in keys(msdict) ? push!(msdict[ms][xidx],fold_change) : msdict[ms]=[[fold_change],Vector{Float64}()]
            end
        end
    end
end

plot_halves=Vector{Plot}()
for ((sib_dict, rys_dict),assay_name) in zip([(rand_sib,rand_rys),(dt_sib,dt_rys)],["Total transcript", "Polyadenylated Transcript"])
    subplots=Vector{Plot}()
    for (n,ms) in enumerate(msvec)
        plt=0.
        mx=0.
        mnn=Inf
        for (dict,color,symbol,leg) in zip([sib_dict, rys_dict],[:green,:darkmagenta],[:circle,:diamond],["Sib","rys"])
            means=[mean(valvec) for valvec in dict[ms]]
            stdevs=[std(valvec) for valvec in dict[ms]]
            upper=Vector{Float64}()
            mns=[mn.val for mn in means]
            stds=[sd.val for sd in stdevs]
            lower=Vector{Float64}()

            for (mn,sd) in zip(mns,stds)
                dist=Normal(mn,sd)
                push!(upper,quantile(dist,.975)-mn)
                push!(lower,mn-quantile(dist,.025))
            end

            xs=vcat([[X[n] for i in 1:length(dict[ms][n])] for n in 1:length(X)]...)
            ys=vcat([[m.val for m in dict[ms][n]] for n in 1:length(X)]...)

            mx=max(mx,maximum(vcat(upper,ys)))
            mnn=min(mnn,minimum(vcat(lower,ys)))

            if dict===sib_dict
                if n!=4 && n!=1
                    plt=scatter(xs,ys, marker=symbol, color=color, markerstrokecolor=color, markersize=3, label=leg*" data", showaxis=:y, xticks=X, xformatter=_->"", ylabel="Fold change",xlims=[5.5,8.5],legend=:none)
                elseif n==1
                    plt=scatter(xs,ys, marker=symbol, color=color, markerstrokecolor=color, markersize=3, label=leg*" data", showaxis=:y, xticks=X, xformatter=_->"", ylabel="Fold change",xlims=[5.5,8.5],legend=:topleft, title=assay_name)
                else 
                    plt=scatter(xs,ys, marker=symbol, color=color, markerstrokecolor=color, markersize=3, label=leg*" data", xticks=X, ylabel="Fold change", xlabel="Age (dpf)", xlims=[5.5,8.5],legend=:none)
                end
            else
                scatter!(plt, xs,ys, marker=symbol, color=color, markerstrokecolor=color, markersize=3, label=leg)
            end

            plot!(plt, X, mns, ribbon=(lower,upper), color=color, label=leg*" mean")
        end
        plot!(X, [1. for i in 1:length(X)], linewidth=4., style=:dot, color=:black, label="WT std.")
        annotate!([(5.75,mnn+.25*(mx-mnn),Plots.text(ms,16))])
        println(ms)
        println("max: $mx min: $mnn")
        println(mnn+.25*(mx-mnn))
        push!(subplots,plt)
    end

    assay_combined=plot(subplots...,layout=grid(4,1),link=:x)
    push!(plot_halves,assay_combined)
end

combined=plot(plot_halves...,layout=grid(1,2),size=(900,900))
savefig(combined,"/bench/PhD/Thesis/images/rys/qPCR.png")

stdmat=zeros(4,2,2,2)

for ((sib_dict, rys_dict),assay_idx) in zip([(rand_sib,rand_rys),(dt_sib,dt_rys)],[1,2])
    for (n,ms) in enumerate(msvec)
        sib_means=0.
        for dict in [sib_dict, rys_dict]
            means=[measurement(mean(valvec).val,std(valvec).val) for valvec in dict[ms]]
            if dict===sib_dict
                sib_means=means
            else
                sdifs=[measurement(mean(valvec).val,std(valvec).val) for valvec in dict[ms]].-sib_means
                stds=[v.val/v.err for v in sdifs]
                stdmat[n,assay_idx,1,:]=stds

                wtdifs=[measurement(mean(valvec).val,std(valvec).val) for valvec in dict[ms]].-1.
                stds=[v.val/v.err for v in wtdifs]
                stdmat[n,assay_idx,2,:]=stds
            end
        end
    end
end

for (a,assay) in enumerate(["Total","polyA"])
    for (m,ms) in enumerate(msvec)
        for (d,age) in enumerate([6,8])
            println("$assay & $ms & $age & $(round(stdmat[m,a,1,d],digits=2)) & $(round(stdmat[m,a,2,d],digits=2))\\\\ \\hline")
        end
    end
end

rr_h2a=[Normal(m,s) for (m,s) in zip(mean.(rand_rys["H2A"]),std.(rand_rys["H2A"]))]
rs_h2a=[Normal(m,s) for (m,s) in zip(mean.(rand_sib["H2A"]),std.(rand_sib["H2A"]))]
rand_h2a_joint=exp(sum(logccdf.(rr_h2a,[d.μ for d in rs_h2a])))

rr_h2b=[Normal(m,s) for (m,s) in zip(mean.(rand_rys["H2B"]),std.(rand_rys["H2B"]))]
rs_h2b=[Normal(m,s) for (m,s) in zip(mean.(rand_sib["H2B"]),std.(rand_sib["H2B"]))]
rand_h2b_joint=exp(sum(logccdf.(rr_h2b,[d.μ for d in rs_h2b])))

dtr_h2a=[Normal(m,s) for (m,s) in zip(mean.(dt_rys["H2A"]),std.(dt_rys["H2A"]))]
dts_h2a=[Normal(m,s) for (m,s) in zip(mean.(dt_sib["H2A"]),std.(dt_sib["H2A"]))]
dt_h2a_joint=exp(sum(logccdf.(dtr_h2a,[d.μ for d in dts_h2a])))

dtr_h2b=[Normal(m,s) for (m,s) in zip(mean.(dt_rys["H2B"]),std.(dt_rys["H2B"]))]
dts_h2b=[Normal(m,s) for (m,s) in zip(mean.(dt_sib["H2B"]),std.(dt_sib["H2B"]))]
dt_h2b_joint=exp(sum(logccdf.(dtr_h2b,[d.μ for d in dts_h2b])))

all_joint=rand_h2b_joint*dt_h2a_joint*dt_h2b_joint