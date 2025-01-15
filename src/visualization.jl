function plot_msf(msf_exp, computed_msf)
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, msf_exp; label="Experimental MSF")
    lines!(ax, computed_msf; label="Computed MSF")
    axislegend(ax)
    ylims!(ax, 0, 0.2)
    return fig
end

