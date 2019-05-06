from vis.heatmap import main, total_df
from project_specific import params
import pandas as pd
import plotly

data = pd.read_csv('/home/liaoth/project/shenzhen_Acinetobacter/ad_analysis/whole_region_roary/heatmap_data.csv', index_col=0)
data = data.replace('-', pd.np.nan)
data = data.iloc[:-1, :]
roary_data_col = data.columns

metadata = total_df.loc[:, params.res_cols]
metadata = metadata.reindex(data.index)

map_dict = {}
map_dict.update(params.res2fun)
map_dict.update(params.vf2fun)

total_data = pd.concat([data, metadata], axis=1)
fig = main(total_data, roary_data_col,
           filter_sample=True,
           tree_p="/home/liaoth/project/shenzhen_Acinetobacter/roary_o/core_gene.newick",
           accessory_cols=params.res_cols,
           up_cluster=None,
           width=3000
           )

plotly.offline.plot(fig,
                    filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/whole_region_roary/21_with_metadata.html",
                    auto_open=False)
if __name__ == '__main__':
    pass
