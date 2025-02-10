using QFiND
using CairoMakie
using LaTeXStrings 

function plot_res(t::AbstractVector{<:Real},
                  exact::AbstractVector{<:Number},
                  approx::AbstractVector{<:Number},
                  err::AbstractVector{<:Number})
    # Style parameters (matching Python style)
    ls1    = 20      # label font size
    ls2    = 15      # tick label font size
    lw1    = 2.5     # line width for main lines
    lw2    = 2.5     # line width for reference lines
    color1 = :orangered
    color2 = :royalblue
    color3 = :black

    # Create a figure with 2 rows and 1 column.
    fig = Figure(size = (800, 700))

    # Top axis: Bath correlation function (BCF) plot.
    ax1 = Axis(fig[1, 1],
        xlabel = "",
        ylabel = L"C(t)",
        xlabelsize = ls2,
        ylabelsize = ls2
    )

    # Bottom axis: Error plot.
    ax2 = Axis(fig[2, 1],
        xlabel = L"t (\mathrm{fs})",
        ylabel = L"\delta C(t)",
        xlabelsize = ls2,
        ylabelsize = ls2
    )

    # Plot on the top axis: Approximate BCF (real and imaginary parts)
    lines!(ax1, t, real.(approx),
        label = L"\mathrm{Re}\,C(t)",
        color = color1,
        linewidth = lw1
    )
    lines!(ax1, t, imag.(approx),
        label = L"\mathrm{Im}\,C(t)",
        color = color2,
        linewidth = lw1
    )

    # Plot on the top axis: Reference (exact) BCF as dashed lines.
    lines!(ax1, t, real.(exact),
        label = L"\text{Reference}",
        color = color3,
        linestyle = :dash,
        linewidth = lw2
    )
    lines!(ax1, t, imag.(exact),
        color = color3,
        linestyle = :dash,
        linewidth = lw2
    )

    # Plot on the bottom axis: Error (real and imaginary parts)
    lines!(ax2, t, real.(err),
        color = color1,
        linewidth = lw1
    )
    lines!(ax2, t, imag.(err),
        color = color2,
        linewidth = lw1
    )

    # Add a legend to the top axis at the top-right position.
    axislegend(ax1, position = :rt, labelsize = ls2)

    # Save the figure as a PNG file.
    save("result.png", fig)

    return fig
end
