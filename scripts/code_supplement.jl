function crawl_minted(apdx_pth::String, path_labels::Vector{Pair{String,String}}, types::Dict{String,String})
    io = open(apdx_pth, "w")
    println(io, "\\chapter{Code Appendix}")
    println(io, "\\setmonofont[Contextuals={Alternate}]{Fira Code}")
    for (path,label) in path_labels
        minty_walk(io, path, label, types)
    end
    close(io) 
end

function minty_walk(io, path, label, types)
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

crawl_minted(ca_pth,pls,types)