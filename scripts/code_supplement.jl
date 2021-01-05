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
    "CMZ_lh.jl" => "ssec:CMZlh",
    "a10periodisation.jl" => "ssec:a10periodisation",
    "BBM_sample_prep.jl" => "ssec:BBMsampleprep",
    "BBM_refinement_prep.jl" => "ssec:BBMrefinementprep",
    "BBM_survey.jl" => "ssec:BBMsurvey",
    "BBM_survey_refinement.jl" => "ssec:BBMsurveyrefinement",
    "BBM_survey_analysis.jl" => "ssec:BBMsurveyanalysis",
    "BBM_refinement_analysis.jl" => "ssec:BBMrefinementanalysis"
)

crawl_minted(ca_pth,pls,types, filelabels)