using ExpFit
using Documenter

DocMeta.setdocmeta!(ExpFit, :DocTestSetup, :(using ExpFit); recursive=true)

makedocs(;
    modules=[ExpFit],
    authors="Hideaki Takahashi <takahashi.hideaki.w33@kyoto-u.jp> and contributors",
    sitename="ExpFit.jl",
    format=Documenter.HTML(;
        canonical="https://doc-package.github.io/ExpFit.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DOC-Package/ExpFit.jl",
    devbranch="main",
)
