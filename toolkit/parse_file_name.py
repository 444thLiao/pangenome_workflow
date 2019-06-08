"""
Basic function for parse input sample name to formatted format in order to import into pipelines.
:type Parse_Function
"""
import os

import pandas as pd


class fileparser():
    def __init__(self, filename):
        filename = os.path.abspath(filename)

        self.df = pd.read_csv(filename, sep='\t', index_col=0, dtype=str)
        self.cols, self.df = validate_df(self.df)

    def get_attr(self, col):
        if col == self.df.index.name:
            return list(self.df.index)
        if col not in self.cols:
            raise Exception("attr %s not in input df" % col)
        else:

            return self.df[col].to_dict()

    def get_full_PE_info(self):

        PE_rows = self.get_PE_rows()
        if PE_rows.shape[1] > 2:
            other_info = PE_rows.to_dict(orient='index')
        else:
            other_info = None
        return other_info

    def get_PE_rows(self):
        input_df = self.df
        PE_rows = input_df.loc[(~input_df.iloc[:, 0].isna()) & (~input_df.iloc[:, 1].isna()), :]
        # the rows which R1 and R2 are both not null
        return PE_rows

    def get_PE_info(self):
        PE_rows = self.get_PE_rows()
        pairreads = tuple(zip(PE_rows.index,
                              PE_rows.iloc[:, 0],
                              PE_rows.iloc[:, 1], ))
        return pairreads

    def get_SE_info(self):
        # except the PE_rows
        input_df = self.df
        PE_rows = self.get_PE_rows()
        Single_rows = input_df.loc[~input_df.index.isin(PE_rows.index), :]
        singlereads = tuple(zip(Single_rows.index,
                                Single_rows.iloc[:, 0]))
        return singlereads


def validate_df(df):
    template_file = os.path.join(os.path.dirname(__file__),
                                 "data_input.template")
    columns_values = open(template_file).read().strip('\n').split('\t')

    if set(df.columns) != set(columns_values):
        raise Exception("INPUT file has unknown header")

    if df.loc[df.index.drop_duplicates(), :].shape[0] != df.shape[0]:
        raise Exception("sample name contains duplicates")
    both_null = df.loc[(df.iloc[:, 0].isna()) & (df.iloc[:, 1].isna()), :]
    if both_null.shape[0] != 0:
        raise Exception("Some rows doesn't have R1 and R2. ")
    return columns_values, df
