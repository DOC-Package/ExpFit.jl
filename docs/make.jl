using ExpFit
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(ExpFit, :DocTestSetup, :(using ExpFit); recursive=true)

#bib = CitationBibliography(
#    joinpath(@__DIR__, "src", "refs.bib");
#    style=:numeric
#)

makedocs(;
    modules=[ExpFit],
    authors="Hideaki Takahashi <takahashi.hideaki.w33@kyoto-u.jp> and contributors",
    sitename="ExpFit.jl",
    format=Documenter.HTML(;
        canonical="https://doc-package.github.io/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Exponential Fitting" => "expfit.md",
        "Exponential Reduction" => "expred.md",
    ],
    draft=get(ENV, "CI", "false") == "false",
    #plugins=[bib]
)

deploydocs(;
    repo="github.com/DOC-Package/ExpFit.jl",
    devbranch="main",
)
