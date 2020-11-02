using CSV,DataFrames,Distributions,NGRefTools,StatsBase,Plots
gr()

a19_1pth="/bench/PhD/datasets/A19GR1.csv"
a19_2pth="/bench/PhD/datasets/A19GR2.csv"
a19_3pth="/bench/PhD/datasets/A19GR3.csv"

gr1df=DataFrame(CSV.read(a19_1pth))
gr2df=DataFrame(CSV.read(a19_2pth))
gr3df=DataFrame(CSV.read(a19_3pth))

X=sort(unique(gr1df."Pulse age (dpf)"))
measure_dict=Dict{String,Vector{Vector{Float64}}}()
measure_dict["ONL"]=[zeros(0) for i in 1:length(X)]
measure_dict["INL"]=[zeros(0) for i in 1:length(X)]
measure_dict["GCL"]=[zeros(0) for i in 1:length(X)]
measure_dict["Isl"]=[zeros(0) for i in 1:length(X)]
measure_dict["Pax6"]=[zeros(0) for i in 1:length(X)]
measure_dict["Isl/Pax6"]=[zeros(0) for i in 1:length(X)]
measure_dict["INL-Pax6"]=[zeros(0) for i in 1:length(X)]
measure_dict["PKCB"]=[zeros(0) for i in 1:length(X)]
measure_dict["GS"]=[zeros(0) for i in 1:length(X)]
measure_dict["HM"]=[zeros(0) for i in 1:length(X)]
measure_dict["Zpr1"]=[zeros(0) for i in 1:length(X)]

function check_layerfrac(row)
    !ismissing(row."#GCL") && !ismissing(row."#INL") && !ismissing(row."#PRL") && !ismissing(row."Total") ? (return true) : (return false)
end

function gr1_dfcheck(row)
    !ismissing(row."GCL - Isl") && !ismissing(row."GCL - Pax6") && !ismissing(row."GCL - Isl/Pax6") && !ismissing(row."INL - Pax6") ? (return true) : (return false)
end

for df in [gr1df,gr2df,gr3df]
    for row in eachrow(df)
        xidx=findfirst(i->i==row."Pulse age (dpf)",X)
        if df === gr1df
            if check_layerfrac(row)
                push!(measure_dict["GCL"][xidx],row."#GCL"/row."Total")
                push!(measure_dict["INL"][xidx],row."#INL"/row."Total")
                push!(measure_dict["ONL"][xidx],row."#PRL"/row."Total")
            end
            if gr1_dfcheck(row)
                push!(measure_dict["Isl"][xidx],row."GCL - Isl"/row."#GCL")
                push!(measure_dict["Pax6"][xidx],row."GCL - Pax6"/row."#GCL")
                push!(measure_dict["Isl/Pax6"][xidx],row."GCL - Isl/Pax6"/row."#GCL")
                push!(measure_dict["INL-Pax6"][xidx],row."INL - Pax6"/row."#INL")
            end
        else
            push!(measure_dict["GCL"][xidx],row."#GCL"/row."Total")
            push!(measure_dict["INL"][xidx],row."#INL"/row."Total")
            push!(measure_dict["ONL"][xidx],row."#PRL"/row."Total")
            if df === gr2df
                push!(measure_dict["PKCB"][xidx],row."INL-PKCB"/row."#INL")
                push!(measure_dict["GS"][xidx],row."INL-GS"/row."#INL")
                push!(measure_dict["HM"][xidx],row."INL-HM"/row."#INL")
            elseif df === gr3df
                push!(measure_dict["Zpr1"][xidx],row."PRL - Zpr1"/row."#PRL")
            end
        end
    end
end