"""
@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Objective: Plot 3D family of orbits with a colorbar using Plotly
"""
import matplotlib
import plotly.graph_objs as go
import plotly.io as pio

pio.renderers.default = "browser"


def plot_orbits(
    mu,
    targeted_po_orbits,
    colourby,
    cb_label="JC",
    colourmap="plasma",
    title=None,
    data_trace=None,
    save=False,
):
    """
    Plot 3D family of orbits with a colorbar
    Dhruv Jain, 2 March 2022

    Parameters
    ----------
    mu :  float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    targeted_po_orbits : List of dictionary
        Targeted Periodic Orbit Information - Dictionary
    colourby : list of float
        Quantity by which to color orbit family members
    cb_label : string, optional
        Quantity name for colorbar. The default is 'JC'.
    colourmap : colourmap, optional
        The default is 'plasma'.
    title : string, optional
        Title of plot, will also be used as file name if save==True. The default is None.
    data_trace : dictionay of plotly.graph_objs object, optional
        To add additional traces
    save : boolean, optional
        Choose if save the plot or not. The default is False.
    """

    layout = go.Layout(autosize=True)
    layout.scene = dict(
        xaxis=dict(
            title="x [nd]",
            tickfont=dict(size=14),
        ),
        yaxis=dict(
            title="y [nd]",
            tickfont=dict(size=14),
        ),
        zaxis=dict(
            title="z [nd]",
            tickfont=dict(size=14),
        ),
        aspectmode="data",
        aspectratio={"x": 1, "y": 1, "z": 1},
    )

    colour_range = [min(colourby), max(colourby)]
    cmap = matplotlib.cm.get_cmap(colourmap)
    norm = matplotlib.colors.Normalize(vmin=colour_range[0], vmax=colour_range[1])
    colours = []
    for i in colourby:
        colours.append(cmap(norm(i)))

    fig = go.Figure(layout=layout)

    for i in range(len(targeted_po_orbits)):
        fig.add_trace(
            go.Scatter3d(
                x=targeted_po_orbits[i]["states"][:, 0],
                y=targeted_po_orbits[i]["states"][:, 1],
                z=targeted_po_orbits[i]["states"][:, 2],
                mode="lines",
                hovertext=f"c={colourby[i]:.3f}",
                line=go.scatter3d.Line(color=f"rgba{colours[i]}", width=5),
            )
        )

    colorbar_trace = go.Scatter3d(
        x=[None],
        y=[None],
        z=[None],
        mode="markers",
        marker=dict(
            colorscale=colourmap,
            showscale=True,
            cmin=colour_range[0],
            cmax=colour_range[1],
            colorbar=dict(
                title=dict(text=cb_label, font=dict(size=18)),
                thickness=10,
                outlinewidth=0,
                tickfont=dict(size=14),
            ),
        ),
        hoverinfo="none",
    )
    fig.add_trace(colorbar_trace)
    if data_trace != None:
        for i in range(len(data_trace)):
            fig.add_trace(data_trace[i])

    fig.update_layout(
        title=title,
        title_x=0.5,
        font=dict(size=18),
        template="plotly_dark",
        showlegend=False,
    )

    if save == True:
        fig.write_html(title + ".html")
    else:
        fig.show()

    return None
