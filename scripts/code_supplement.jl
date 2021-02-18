function crawl_minted(apdx_pth::String, path_labels::Vector{Pair{String,String}}, types::Dict{String,String}, filelabels::Dict{String,String})
    io = open(apdx_pth, "w")
    println(io, "\\chapter{Code Appendix}")
    println(io, "\\section{Code and Output Archive}")
    println(io, "\\label{sec:archive}")

    println(io, "\\setmonofont[Contextuals={Alternate}]{Fira Code}")
    for (path,label) in path_labels
        minty_walk(io, path, label, types, filelabels)
    end
    close(io) 
end

function minty_walk(io, path, label, types, filelabels)
    println(io)
    println(io, "\\section{\\protect\\path{$(basename(path))}}")
    println(io, "\\label{$label}")
    for (root, dirs, files) in walkdir(path)
        for file in files
            if occursin('.',file)
                ext=file[findlast(isequal('.'),file):end]
                if ext in keys(types)
                    type=types[ext]
                    println(io, "\\subsection{\\protect\\path{$(joinpath(root[length(path)+1:end],file))}}")
                    if file in keys(filelabels)
                        println(io, "\\label{$(filelabels[file])}")
                    end
                    println(io, "\\inputminted[breaklines, mathescape, linenos,   numbersep=5pt, frame=lines, framesep=2mm]{$type}{$(joinpath(root,file))}")
                end
            end
        end
    end
end

ca_pth="/bench/PhD/Thesis/chapters/PTIII/code.tex"
pls=[
    "/srv/git/SMME" => "SMMEcode",
    "/srv/git/NGRefTools" => "NGRefTools",
    "/srv/git/GMC_NS" => "GMCNScode",
    "/srv/git/CMZNicheSims" => "CMZNScode",
    "/srv/git/BioBackgroundModels" => "BBMcode",
    "/srv/git/BioMotifInference" => "BMIcode",
    "/srv/git/AWSWrangler" => "AWSWrangler",
    "/bench/PhD/Thesis/Analyses" => "analysiscode"
    ]

types=Dict(
    ".hpp" => "cpp",
    ".cpp" => "cpp",
    ".jl" => "julia",
    ".py" => "python"
)

filelabels=Dict(
    "He_output_plot.py" => "ssec:He_output_plot",
    "Wan_output_plot.py" => "ssec:Wan_output_plot",
    "Mitotic_rate_plot.py" => "ssec:Mitotic_rate_plot",
    "CMZ_lh.jl" => "ssec:CMZlh",
    "lens_model.jl" => "ssec:slicelensmodel",
    "a10popsurvey.jl"=>"ssec:a10popsurvey",
    "a10correlation.jl"=>"ssec:a10correlation",
    "a10periodisation.jl" => "ssec:a10periodisation",
    "a10dvratio.jl" => "ssec:a10dvratio",
    "a10dvslice.jl"=>"ssec:a10dvslice",
    "a10nvln.jl"=>"ssec:a10nvln",
    "a19lineagetrace.jl"=>"ssec:a19lineagetrace",
    "a25dvlinreg.jl" => "ssec:a25dvlinreg",
    "a27linreg.jl" => "ssec:a27linreg",
    "a27GMC_NS.jl" => "ssec:a27GMC_NS",
    "BBM_sample_prep.jl" => "ssec:BBMsampleprep",
    "BBM_refinement_prep.jl" => "ssec:BBMrefinementprep",
    "BBM_survey.jl" => "ssec:BBMsurvey",
    "BBM_survey_refinement.jl" => "ssec:BBMsurveyrefinement",
    "BBM_survey_analysis.jl" => "ssec:BBMsurveyanalysis",
    "BBM_refinement_analysis.jl" => "ssec:BBMrefinementanalysis",
    "IPM_likelihood.jl" => "ssec:IPMlikelihood",
    "ryspont.jl" => "ssec:ryspont",
    "rysp_GMC_NS.jl" => "ssec:rysp_GMC_NS",
    "a38.jl" => "ssec:a38",
    "a38GMC_NS.jl" => "ssec:a38GMC_NS",
    "caspase.jl"=>"ssec:caspasescript",
    "qPCR.jl"=>"ssec:qPCR",
    "occupancy.py"=>"ssec:occupancy",
    "position_overlap.jl"=>"ssec:position_overlap",
    "dif_pos_learner.jl" => "ssec:dif_pos_learner"
)

crawl_minted(ca_pth,pls,types, filelabels)