function crawl_minted(apdx_pth::String, path_labels::Vector{Pair{String,Vector{String}}}, types::Dict{String,String}, filelabels::Dict{String,String})
    io = open(apdx_pth, "w")
    println(io, "\\chapter{Code Appendix}")
    println(io, "\\label{chap:code}")

    println(io, "\\setmonofont[Contextuals={Alternate}]{Fira Code}")
    for (path,labels) in path_labels
        label,gitlabel=labels
        minty_walk(io, path, label, gitlabel, types, filelabels)
    end
    close(io) 
end

function minty_walk(io, path, label, gitlabel, types, filelabels)
    println(io)
    println(io, "\\section{\\protect\\path{$(basename(path))}}")
    println(io, "Github repository: \\url{$gitlabel}")
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
    "/srv/git/SMME" => ["SMMEcode","https://github.com/mmattocks/SMME"],
    "/srv/git/NGRefTools" => ["NGRefTools","https://github.com/mmattocks/NGRefTools.jl"],
    "/srv/git/GMC_NS" => ["GMCNScode", "https://github.com/mmattocks/GMC_NS.jl"],
    "/srv/git/CMZNicheSims" => ["CMZNScode", "https://github.com/mmattocks/CMZNicheSims.jl"],
    "/srv/git/BioBackgroundModels" => ["BBMcode","https://github.com/mmattocks/BioBackgroundModels.jl"],
    "/srv/git/BioMotifInference" => ["BMIcode", "https://github.com/mmattocks/BioMotifInference.jl"],
    "/srv/git/AWSWrangler" => ["AWSWrangler", "https://github.com/mmattocks/AWSWrangler.jl"],
    "/bench/PhD/Thesis/Analyses" => ["analysiscode", "https://github.com/mmattocks/Thesis/Analyses"]
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
    "SPSA_fixture.py" => "ssec:SPSAfix",
    "CMZ_lh.jl" => "ssec:CMZlh",
    "lens_model.jl" => "ssec:slicelensmodel",
    "a10popsurvey.jl"=>"ssec:a10popsurvey",
    "a10correlation.jl"=>"ssec:a10correlation",
    "a10periodisation.jl" => "ssec:a10periodisation",
    "a10dvratio.jl" => "ssec:a10dvratio",
    "a10dvslice.jl"=>"ssec:a10dvslice",
    "a10dvdecayslice.jl"=>"ssec:a10dvdecayslice",
    "a10nvln.jl"=>"ssec:a10nvln",
    "a19lineagetrace.jl"=>"ssec:a19lineagetrace",
    "a25dvlinreg.jl" => "ssec:a25dvlinreg",
    "a25thymidinesim.jl"=>"ssec:a25thymidinesim",
    "a27linreg.jl" => "ssec:a27linreg",
    "a27GMC_NS.jl" => "ssec:a27GMC_NS",
    "a35thymidinesim.jl"=>"ssec:a35thymidinesim",
    "BBM_sample_prep.jl" => "ssec:BBMsampleprep",
    "BBM_refinement_prep.jl" => "ssec:BBMrefinementprep",
    "BBM_survey.jl" => "ssec:BBMsurvey",
    "BBM_survey_refinement.jl" => "ssec:BBMsurveyrefinement",
    "BBM_survey_analysis.jl" => "ssec:BBMsurveyanalysis",
    "BBM_refinement_analysis.jl" => "ssec:BBMrefinementanalysis",
    "baum-welch.jl" => "ssec:BaumWelch",
    "churbanov.jl" => "ssec:Churbanov",
    "IPM_likelihood.jl" => "ssec:IPMlikelihood",
    "ryspont.jl" => "ssec:ryspont",
    "rysp_GMC_NS.jl" => "ssec:rysp_GMC_NS",
    "a38.jl" => "ssec:a38",
    "a38GMC_NS.jl" => "ssec:a38GMC_NS",
    "caspase.jl"=>"ssec:caspasescript",
    "qPCR.jl"=>"ssec:qPCR",
    "occupancy.py"=>"ssec:occupancy",
    "position_overlap.jl"=>"ssec:position_overlap",
    "dif_pos_learner.jl" => "ssec:dif_pos_learner",
    "spike_recovery.jl" => "ssec:spike_recovery",
    "dif_pos_sample_prep.jl"=>"ssec:dif_pos_sample_prep",
    "dif_pos_assembly.jl"=>"ssec:dif_pos_assembly"
)

crawl_minted(ca_pth,pls,types, filelabels)