import plotly
import plotly.graph_objs as go
import plotly.figure_factory as ff
from sklearn.decomposition import PCA
import pandas as pd




def draw_PCA_with_envfit_arrow(raw_ord, fit_result:pd.DataFrame,base_length=3):

    fig = go.Figure()

    fig.add_scatter(x=raw_ord.iloc[:, 0],
                    y=raw_ord.iloc[:, 1],
                    text=raw_ord.index,
                    marker=dict(size=20, opacity=0.7),
                    mode="markers", )
    fit_result = fit_result.loc[fit_result.iloc[:,1]<=0.05,:]
    for idx,row in fit_result.iterrows() :
        scale = row["r2"] * base_length
        source,end = row["Source"]*scale,row["End"]*scale
        fig.add_scatter(x=[0,source],
                        y=[0,end],
                        mode="lines",
                        showlegend=False
                        )
        fig.add_scatter(x=[source],y=[end],
                        mode="markers+text",
                        text=[idx],
                        showlegend=False,

                        textposition="middle right")
    fig.layout.hovermode = "closest"
    fig.layout.font.size = 15
    return fig