
import xgboost
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.preprocessing import LabelEncoder
import pandas as pd
from sklearn.metrics import confusion_matrix

data = pd.read_csv('/home/liaoth/data2/project/shenzhen_actineto/roary_o/gene_presence_absence.csv',index_col=0)
sub_df = data.loc[:,['34960', '35082', '35989', '37502', '37706', '39054', '39058', '39292',
       '39502', '39528', '40067', '40116', '40283', '40615', '41194', '41296',
       '41833', '42121', '42268', '42349', '43458',]]
sub_df = sub_df.loc[~pd.isna(sub_df).all(1),:]
num_df = sub_df.applymap(lambda x:abs(int(bool(pd.isnull(x)))-1) ).T

metadata = pd.read_csv("/home/liaoth/data2/project/shenzhen_actineto/汇总结果/metadata.csv",index_col=0)
sub_meta = metadata.loc[num_df.index,'source']

le = LabelEncoder()
rfs = RandomForestClassifier(n_estimators=2,max_depth=2)
rfs.fit(num_df,le.fit_transform(sub_meta))

print(confusion_matrix(le.fit_transform(sub_meta), rfs.predict(num_df)))


dbc = GradientBoostingClassifier(n_estimators=2,max_depth=2)
dbc.fit(num_df,le.fit_transform(sub_meta))
print(confusion_matrix(le.fit_transform(sub_meta), dbc.predict(num_df)))