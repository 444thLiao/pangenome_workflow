# require plotly >= 4.0.0

import click
import pandas as pd
import plotly.express as px
from skbio.stats.ordination import pcoa


@click.command()
@click.option("-i", "--input", "dist_file")
@click.option("-m", "--metadata", "metadata")
@click.option("-col", "col_name")
@click.option("-o", "output")
@click.option("-w", "width", default=1920, required=False)
@click.option("-h", "height", default=1080, required=False)
def cli(dist_file, metadata, col_name, output,
        width, height):
    dist = pd.read_csv(dist_file, index_col=0)
    metadata = pd.read_csv(metadata, index_col=0,dtype=str)
    columns = '\t'.join(map(str, metadata.columns))
    if col_name not in metadata.columns:
        raise Exception(f"requested '{col_name}' is not in the given metadata table."
                        f"validated columns are {columns}")
    col = metadata.loc[:, col_name]
    fig = draw_PCoA(dist, col)
    fig.layout.width = width
    fig.layout.height = height
    fig.write_html(file=output, auto_open=True)


def draw_PCoA(dist, col):
    """
    :param dist:
    :param col: should be series or a order list
    :return:
    """
    pcoa_result = pcoa(dist)
    ord_data = pcoa_result.samples
    ord_data.index = dist.index
    ord_data.loc[:, "sample ID"] = ord_data.index
    if isinstance(col, pd.Series):
        col = col.loc[dist.index]
    elif isinstance(col, list):
        pass
    ord_data.loc[:, 'annotated'] = col
    ord_data.loc[:, 'size'] = 20
    if len(set(col)) > 10:
        color_discrete_sequence = px.colors.qualitative.Dark24
    else:
        color_discrete_sequence = px.colors.qualitative.Plotly

    fig = px.scatter(ord_data,
                     x="PC1",
                     y="PC2",
                     size='size',
                     opacity=0.6,
                     color="annotated",
                     marginal_x="violin",
                     marginal_y="violin",
                     hover_data=["sample ID", 'annotated'],
                     # text="sample ID",
                     color_discrete_sequence=color_discrete_sequence
                     )
    fig.layout.hovermode = "closest"
    fig.layout.xaxis.title = "NMDS1({:.2f}%)".format(pcoa_result.proportion_explained[0] * 100)
    fig.layout.yaxis.title = "NMDS2({:.2f}%)".format(pcoa_result.proportion_explained[1] * 100)
    fig.layout.font.size = 15
    fig.layout.font.size = 15
    return fig


if __name__ == '__main__':
    cli()

    "python3 /home/liaoth/data2/project/genome_pipelines/vis/draw_pcoa.py -i /home/liaoth/data2/project/shenzhen_Acinetobacter/pipelines_190619/summary_output/pairwise_mash/pairwise_mash.dist -m /home/liaoth/data2/project/shenzhen_Acinetobacter/pipelines_190619/summary_output/species_annotated.csv -col 'annotated organism' -o /home/liaoth/Desktop/test.html"
    # dist = pd.read_csv('/home/liaoth/data2/project/shenzhen_Acinetobacter/pipelines_190619/summary_output/pairwise_mash/pairwise_mash.dist',
    #                     index_col=0)
    # species_annoated = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/pipelines_190619/summary_output/species_annotated.csv",
    #                                index_col=0)
    #
    # pcoa_result = pcoa(dist)
    # ord_data = pcoa_result.samples
    # ord_data.index = dist.index
    # ord_data.loc[:,'species'] = species_annoated.loc[ord_data.index,'annotated organism']
    # ord_data.loc[:,'size'] = 20
    # fig = px.scatter(ord_data,
    #                  x="PC1",
    #                  y="PC2",
    #                  size='size',
    #                  opacity=0.6,
    #                  color="species",
    #                  )
    # fig.layout.hovermode = "closest"
    # fig.layout.xaxis.title = "NMDS1({:.2f}%)".format(pcoa_result.proportion_explained[0] * 100)
    # fig.layout.yaxis.title = "NMDS2({:.2f}%)".format(pcoa_result.proportion_explained[1] * 100)
    # fig.layout.font.size = 15
    # fig.write_html('./temp_plot.html',auto_open=True)
