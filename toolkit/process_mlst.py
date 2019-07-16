import os
import sys

import pandas as pd
from tqdm import tqdm
from pandas.errors import EmptyDataError
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from pipelines.soft_db_path import mlst_db


def parse_mlst(infiles, sample_name=None):
    collect_df = []
    for infile in infiles:
        indf = pd.read_csv(infile, sep='\t', index_col=None, header=None)
        header = ["raw_path", "Scheme", "ST"]
        header += ["Locus%s" % (_ + 1)
                   for _ in range(indf.shape[1] - len(header))]
        indf.columns = header
        collect_df.append(indf)
    final_df = pd.concat(collect_df, axis=0)
    if sample_name is None:
        final_df.index = range(final_df.shape[0])
    else:
        final_df.index = [sample_name + '.%s' % _
                          for _ in range(final_df.shape[0])]
    return final_df


def extract_gene(g_list):
    # input a row of info output by mlst
    # filterout not locus part with `list(val[[_ for _ in val.index if 'locus' in _.lower()]])`
    gene_list = []
    st_list = []
    for info in g_list:
        gene = info.split('(')[0]
        st = info.split('(')[1].strip(')')
        if ',' in st:
            st = st.split(',')
        else:
            st = [st]
        new_st = []
        for _st in st:
            if '~' in _st:
                new_st.append((_st.strip('~'), '~'))
            elif '?' in _st:
                new_st.append((_st.strip('?'), '?'))
            else:
                new_st.append((_st, ''))
        gene_list.append(gene)
        st_list.append(new_st)
    return gene_list, st_list


def construct_st(num_list, final=[], status=''):
    # input a row of info output by `extract_gene`
    """
    recursive function for get ST which from different combination of genes
    st_list = [[('1', '')],
               [('3', '')],
               [('3', ''), ('189', '')],
               [('2', '')],
               [('2', '~')],
               [('96', '')],
               [('3', '')]]
    construct_st(st_list) == [(['1', '3', '3', '2', '2', '96', '3'], '~'),
                              (['1', '3', '3', '2', '2', '22', '3'], '~'),
                              (['1', '3', '189', '2', '2', '96', '3'], ''),
                              (['1', '3', '189', '2', '2', '22', '3'], '')]
    st_list = [[('1', '')],
               [('3', '')],
               [('3', '~'), ('189', '')],
               [('2', '')],
               [('2', '')],
               [('96', ''), ('22', '')],
               [('3', '')]]
    construct_st(st_list) == [(['1', '3', '1', '3', '3', '2', '2', '96', '3'], '~'),
                              (['1', '3', '1', '3', '3', '2', '2', '22', '3'], '~'),
                              (['1', '3', '1', '3', '189', '2', '2', '96', '3'], ''),
                              (['1', '3', '1', '3', '189', '2', '2', '22', '3'], '')]

    st_list =  [[('40', '~')],
                 [('52', '')],
                 [('32', '~')],
                 [('43', '~')],
                 [('4', '')],
                 [('278', '')],
                 [('3', '?')]]
    construct_st(st_list) == (['40', '52', '32', '43', '4', '278', '3'], '~~~?')
    :param st_list:
    :param final:
    :param status:
    :return:
    """
    st_list = num_list[::]
    if len(st_list) == 0:
        _final = final[::]
        _status = status[::]
        final = []
        status = ''
        return _final, _status
    val = st_list.pop(0)
    if len(val) == 1:
        final.append(val[0][0])
        status += val[0][1]
        return construct_st(st_list, final, status)
    else:
        differt_path = []
        for _val in val:
            _final = final[::]
            _final.append(_val[0])
            # print(_val,_final,st_list)
            _status = status[::]
            _status += _val[1]
            after_ = construct_st(st_list[::], _final, _status)
            if type(after_) == tuple:
                differt_path.append(after_)
            else:
                differt_path += after_
        return differt_path


def format_ST_dict(sub_df):
    ST_dict = sub_df.T.to_dict('list')
    ST_dict = {','.join(map(str, val)): k for k, val in ST_dict.items()}
    return ST_dict


def redefine_mlst(mlst_df: pd.DataFrame, scheme, db=mlst_db):
    db_file = os.path.join(db, scheme, scheme + '.txt')
    try:
        ST_df = pd.read_csv(db_file, sep='\t')
    except EmptyDataError:
        return mlst_df
    if mlst_df.shape[0] == 0:
        return mlst_df
    for idx, val in tqdm(mlst_df.iterrows()):
        gene_cols = list(val[[_
                              for _ in val.index
                              if 'locus' in _.lower()]])
        # get columns contain gene names
        gene_list, extract_st = extract_gene(gene_cols)
        st_list = construct_st(extract_st, [], '')
        st_list = [st_list] if type(st_list) == tuple else st_list

        sub_df = ST_df.loc[:, gene_list]
        sub_df.index = ST_df.loc[:, 'ST']
        ST_dict = format_ST_dict(sub_df)
        # get mapping dict of ST

        scheme_source = gene_list[0].split('_')[0]
        # Oxf or Pas... for abaumannii
        for count, st_s in enumerate(st_list):
            assert len(st_s[0]) == len(gene_list)
            ST = ST_dict.get(','.join(map(str, st_s[0])), '-')
            if ST == '-':
                final_ST = 'NEW ST'
            else:
                final_ST = 'ST%s' % ST
            mlst_df.loc[idx, '%s_ST.%s' % (scheme_source, (count + 1))] = final_ST
    # remove duplicated columns
    may_dup_cols = [_ for _ in mlst_df.columns if '_ST.' in _]
    if may_dup_cols:
        necessary_cols = list(mlst_df.columns[:mlst_df.columns.get_loc(may_dup_cols[0])])
        sub_df = mlst_df.loc[:,may_dup_cols].fillna('?')
        unique_cols = list(sub_df.sum(0).drop_duplicates().index)
        new_mlst_df = mlst_df.loc[:,necessary_cols+unique_cols]
        new_mlst_df.columns = necessary_cols + ['%s_ST.%s' % (scheme_source,_)
                                                for _ in range(1,len(unique_cols)+1)]

    return mlst_df


def main(mlst_df):
    output_mlst = {}
    # init a empty dict
    schemes = set(mlst_df.Scheme)
    scheme_df_dict = {}
    for scheme in schemes:
        scheme_df_dict[scheme] = mlst_df.loc[mlst_df.Scheme == scheme, :]
    # get all used scheme & init a dict which contain multiple df corresponding scheme

    scheme_df_dict = {k: redefine_mlst(_df, scheme=k)
                      for k, _df in scheme_df_dict.items()}

    collect_df = []
    for scheme, _df in scheme_df_dict.items():
        _df.index = [str(_).split('.')[0] for _ in _df.index]
        ST_columns = [_ for _ in _df.columns if "_ST" in _]
        collect_df.append(_df.loc[:, ST_columns])
        output_mlst[scheme] = _df
    merged_df = pd.concat(collect_df, axis=1)
    return output_mlst, merged_df


if __name__ == '__main__':
    import argparse
    from pipelines.tasks import valid_path

    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--infile", help='mlst output df by merged_mlst')
    parse.add_argument("-o", "--outdir", help='')
    parse.add_argument("-p", "--prefix", help='')
    args = parse.parse_args()
    mlst_file = os.path.abspath(args.infile)
    odir = os.path.abspath(args.outdir)
    prefix = args.prefix
    valid_path(mlst_file, check_size=1)
    valid_path(odir, check_odir=1)

    mlst_df = pd.read_csv(mlst_file,
                          sep=',',
                          index_col=0)  # todo: auto detect the sep
    output_mlst = main(mlst_df)
    if prefix is None:
        prefix = os.path.basename(mlst_file).split('.')[0]
    for scheme, odf in output_mlst.items():
        with open(os.path.join(odir,
                               "{sn}_{scheme}_mlst.txt".format(
                                   sn=prefix,
                                   scheme=scheme)),
                  'w') as f1:
            odf.to_csv(f1, index=0)
