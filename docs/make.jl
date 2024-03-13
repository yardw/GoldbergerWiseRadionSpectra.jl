using GoldbergerWiseRadionSpectra
using Documenter

DocMeta.setdocmeta!(GoldbergerWiseRadionSpectra, :DocTestSetup, :(using GoldbergerWiseRadionSpectra); recursive=true)

makedocs(;
    modules=[GoldbergerWiseRadionSpectra],
    authors="Yaoduo WANG <yaoduowang@sjtu.edu.cn> and contributors",
    repo="https://github.com/yardw/GoldbergerWiseRadionSpectra.jl/blob/{commit}{path}#{line}",
    sitename="GoldbergerWiseRadionSpectra.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yardw.github.io/GoldbergerWiseRadionSpectra.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yardw/GoldbergerWiseRadionSpectra.jl",
    devbranch="main",
)
