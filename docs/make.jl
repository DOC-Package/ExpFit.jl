using ExpFit
using Documenter

DocMeta.setdocmeta!(ExpFit, :DocTestSetup, :(using ExpFit); recursive=true)

makedocs(;
    modules=[ExpFit],
    authors="Hideaki Takahashi <hide.kyoto.1020@gmail.com> and contributors",
    sitename="ExpFit.jl",
    format=Documenter.HTML(;
        canonical="https://DOC-package.github.io/ExpFit.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DOC-package/ExpFit.jl",
    devbranch="main",
)
