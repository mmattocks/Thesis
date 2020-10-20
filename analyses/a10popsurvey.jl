using CSV,DataFrames

a10pth="/bench/PhD/datasets/A10 measurements 2018update.csv"

a10df=DataFrame(CSV.read(a10pth))
a10df.BSR=[string(row.Block,',',row.Slide,',',row.Row) for row in eachrow(a10df)]

X=Vector{Float64}()
measure_dict=Dict{String,Vector{Vector{Float64}}}()

measure_dict["CMZ Sum"]=Vector{Vector{Float64}}()
measure_dict["Dorsal CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["Ventral CMZ (#)"]=Vector{Vector{Float64}}()
measure_dict["Retina thickness"]=Vector{Vector{Float64}}()
measure_dict["RPE length"]=Vector{Vector{Float64}}()
measure_dict["Lens diameter"]=Vector{Vector{Float64}}()
measure_dict["ON diameter"]=Vector{Vector{Float64}}()
measure_dict["ONL"]=Vector{Vector{Float64}}()
measure_dict["OPL"]=Vector{Vector{Float64}}()
measure_dict["INL"]=Vector{Vector{Float64}}()
measure_dict["IPL"]=Vector{Vector{Float64}}()
measure_dict["GCL"]=Vector{Vector{Float64}}()

#EXTRACT ALL MEASUREMENTS FROM DF TO VECTORS
for t_df in groupby(a10df, "Time point (d)")
    push!(X,Float64(unique(t_df."Time point (d)")[1]))
    tvec,dvec,vvec,rtvec,rplvec, rptvec, ldvec, ondvec, onlvec, oplvec, inlvec, iplvec, gclvec = [zeros(0) for i in 1:13]

    for i_df in groupby(t_df,"BSR")
        length([skipmissing(i_df."CMZ Sum")...])>0 && push!(tvec,mean(skipmissing(i_df."CMZ Sum")))
        length([skipmissing(i_df."Dorsal CMZ (#)")...])>0 && push!(dvec,mean(skipmissing(i_df."Dorsal CMZ (#)")))
        length([skipmissing(i_df."Ventral CMZ (#)")...])>0 && push!(vvec,mean(skipmissing(i_df."Ventral CMZ (#)")))
        length([skipmissing(i_df."Retina thickness")...])>0 && push!(rtvec,mean(skipmissing(i_df."Retina thickness")))
        length([skipmissing(i_df."RPE length")...])>0 && push!(rptvec,mean(skipmissing(i_df."RPE length")))
        length([skipmissing(i_df."Lens diameter")...])>0 && push!(ldvec,mean(skipmissing(i_df."Lens diameter")))
        length([skipmissing(i_df."ON diameter")...])>0 && push!(ondvec,mean(skipmissing(i_df."ON diameter")))
        length([skipmissing(i_df.ONL)...])>0 && push!(onlvec,mean(skipmissing(i_df.ONL)))   
        length([skipmissing(i_df.OPL)...])>0 && push!(oplvec,mean(skipmissing(i_df.OPL)))  
        length([skipmissing(i_df.INL)...])>0 && push!(inlvec,mean(skipmissing(i_df.INL)))
        length([skipmissing(i_df.IPL)...])>0 && push!(iplvec,mean(skipmissing(i_df.IPL)))  
        length([skipmissing(i_df.GCL)...])>0 && push!(gclvec,mean(skipmissing(i_df.GCL)))   
    end

    push!(measure_dict["CMZ Sum"],tvec)
    push!(measure_dict["Dorsal CMZ (#)"],dvec)
    push!(measure_dict["Ventral CMZ (#)"],vvec)
    push!(measure_dict["Retina thickness"],rtvec)
    push!(measure_dict["RPE length"],rptvec)
    push!(measure_dict["Lens diameter"],ldvec)
    push!(measure_dict["ON diameter"],ondvec)
    push!(measure_dict["ONL"],onlvec)
    push!(measure_dict["OPL"],oplvec)
    push!(measure_dict["INL"],inlvec)
    push!(measure_dict["IPL"],iplvec)
    push!(measure_dict["GCL"],gclvec)
end

