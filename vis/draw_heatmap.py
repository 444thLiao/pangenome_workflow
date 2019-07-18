import click
import plotly.figure_factory as ff
import plotly.graph_objs as go
from plotly.subplots import make_subplots


@click.command()
def cli():
    fig = create_heatmap()


def create_heatmap():
    pass


import pandas as pd

left_dist = pd.read_csv('/home/liaoth/data2/project/shenzhen_Acinetobacter/pipelines_190619/summary_output/pairwise_mash/pairwise_mash.dist',
                        index_col=0)
main_matrix = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/pipelines_190619/all_roary_o/gene_presence_absence.Rtab",
                          sep='\t',
                          index_col=0).T
accessory_matrix = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/pipelines_190619/summary_all_roary_o/samples2annotate.csv",
                               index_col=0)
create_heatmap(left_dist,
               main_matrix,
               )


def create_heatmap(left_dist,
                   main_matrix,
                   up_dist=None,
                   accessory_matrix=None,
                   height=1500,
                   width=None,
                   main_matrix_text=None,
                   accessory_matrix_text=None,
                   return_matrix=False):
    fig = make_subplots(2, 3,
                        specs=[[None, None, {}],
                               [{}, {}, {}]],
                        shared_yaxes=True,
                        shared_xaxes=True)
    # This is the format of your plot grid:
    #     (empty)          (empty)      [ (1,3) x3,y1 ]
    # [ (2,1) x1,y2 ]  [ (2,2) x2,y2 ]  [ (2,3) x3,y2 ]
    sub_df = main_matrix.copy()
    _sub_df = None
    # create top dendrogram
    if up_dist is not None:
        up_dendro = ff.create_dendrogram(up_dist.values,
                                         orientation='bottom',
                                         labels=up_dist.index)

        sub_df = main_matrix.reindex(columns=up_dendro.layout.xaxis.ticktext)
        main_matrix_text = main_matrix_text.reindex(columns=up_dendro.layout.xaxis.ticktext) if main_matrix_text is not None else None
        # coordinate to up dend
        for i in up_dendro.data:
            i.showlegend = False
            fig.append_trace(i, 1, 3)
    # Create left Dendrogram
    if type(left_dist) == go.Figure:
        side_dendro = left_dist
    else:
        side_dendro = ff.create_dendrogram(left_dist.values,
                                           orientation='right',
                                           labels=left_dist.index)
    for i in side_dendro.data:
        i.showlegend = False
        fig.append_trace(i, 2, 1)
    sub_df = sub_df.reindex(side_dendro.layout.yaxis.ticktext)
    main_matrix_text = main_matrix_text.reindex(side_dendro.layout.yaxis.ticktext) if main_matrix_text is not None else None
    # coordinate to left dend
    ############################################################
    # draw (2,3) heatmap
    heatmap = go.Heatmap(
        x=sub_df.columns,
        y=sub_df.index,
        z=sub_df.values,
        text=main_matrix_text.values if main_matrix_text is not None else None,
        colorscale='Earth',
        reversescale=True,
        showscale=False
    )
    if up_dist is not None:
        heatmap['x'] = up_dendro.layout.xaxis.tickvals
    heatmap['y'] = side_dendro.layout.yaxis.tickvals
    # fig.layout.xaxis3.ticktext = sub_df.columns
    # fig.layout.xaxis3.tickangle = 30
    fig.append_trace(heatmap, 2, 3)

    ############################################################
    if accessory_matrix is not None:
        _sub_df = accessory_matrix.reindex(side_dendro.layout.yaxis.ticktext)
        accessory_matrix_text = accessory_matrix_text.reindex(side_dendro.layout.yaxis.ticktext) if accessory_matrix_text is not None else None
        # import pdb;pdb.set_trace()
        heatmap = go.Heatmap(
            x=_sub_df.columns,
            y=_sub_df.index,
            z=_sub_df.values,
            text=accessory_matrix_text.values if accessory_matrix_text is not None else None,
            hoverinfo='all',
            colorscale='Earth',
            opacity=0.9,
            reversescale=True,
            showscale=False
        )
        heatmap['y'] = side_dendro.layout.yaxis.tickvals
        fig.append_trace(heatmap, 2, 2)
    ############################################################

    if accessory_matrix is not None:
        fig.layout.xaxis1.domain = [0, 0.05]
        fig.layout.xaxis2.domain = [0.05, 0.1]
        fig.layout.xaxis3.domain = [0.1, 1]
    else:
        fig.layout.xaxis1.domain = [0, 0.17]
        fig.layout.xaxis2.domain = [0.0, 0.0]
        fig.layout.xaxis3.domain = [0.2, 1]
    if up_dist is not None:
        fig.layout.yaxis1.domain = [0.9, 1]
        fig.layout.yaxis2.domain = [0.0, 0.9]
    else:
        fig.layout.yaxis2.domain = [0.0, 0]
        fig.layout.yaxis1.domain = [0.0, 1]
    fig.layout.yaxis2.ticktext = side_dendro.layout.yaxis.ticktext
    fig.layout.yaxis2.tickvals = side_dendro.layout.yaxis.tickvals
    if up_dist is not None:
        fig.layout.xaxis3.ticktext = up_dendro.layout.xaxis.ticktext
        fig.layout.xaxis3.tickvals = up_dendro.layout.xaxis.tickvals
    # fig.layout.xaxis.anchor = 'x2'
    # fig.layout.margin.b = 250
    fig.layout.width = main_matrix.shape[1] * 13 if width is None else width
    fig.layout.height = height
    fig.layout.hovermode = 'closest'
    fig.layout.xaxis1.zeroline = fig.layout.xaxis2.zeroline = fig.layout.xaxis3.zeroline = False
    fig.layout.yaxis1.zeroline = fig.layout.yaxis2.zeroline = False
    fig.layout.xaxis1.showgrid = False
    fig.layout.yaxis2.showgrid = False
    fig.layout.xaxis1.showgrid = fig.layout.xaxis2.showgrid = fig.layout.xaxis3.showgrid = False
    fig.layout.yaxis1.showgrid = False
    fig.layout.xaxis1.showticklabels = False
    fig.layout.xaxis3.showticklabels = False

    fig.layout.xaxis1.visible = False
    fig.layout.xaxis3.visible = False
    fig.layout.xaxis2.visible = False

    fig.layout.yaxis1.showticklabels = False
    if return_matrix:
        if main_matrix_text is not None or accessory_matrix_text is not None:
            return fig, main_matrix_text, accessory_matrix_text
        else:
            return fig, sub_df, _sub_df
    else:
        return fig


if __name__ == '__main__':
    cli()
