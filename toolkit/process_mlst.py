import os
import sys

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from pipelines.constant_str import mlst_db


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
    ST_df = pd.read_csv(db_file, sep='\t')
    for idx, val in tqdm(mlst_df.iterrows()):
        gene_cols = list(val[[_ for _ in val.index if 'locus' in _.lower()]])
        gene_list, extract_st = extract_gene(gene_cols)
        st_list = construct_st(extract_st, [], '')
        st_list = [st_list] if type(st_list) == tuple else st_list
        sub_df = ST_df.loc[:, gene_list]
        sub_df.index = ST_df.loc[:, 'ST']
        ST_dict = format_ST_dict(sub_df)

        scheme_source = gene_list[0].split('_')[0]
        for count, st_s in enumerate(st_list):
            assert len(st_s[0]) == len(gene_list)
            ST = ST_dict.get(','.join(map(str, st_s[0])), '-')
            if ST == '-':
                final_ST = 'NEW ST'
            else:
                final_ST = 'ST%s' % ST
            mlst_df.loc[idx, '%s_ST.%s' % (scheme_source, (count + 1))] = final_ST
    return mlst_df


def main(mlst_df):
    output_mlst = {}
    schemes = set([col.split('.')[-1] for col in mlst_df.columns if '.' in col])
    scheme_df_dict = {}
    sub_df = mlst_df.loc[:, [col for col in mlst_df.columns if '.' not in col]]
    scheme = sub_df.loc[sub_df.index[0], [_ for _ in sub_df.columns if 'Scheme' in _]][0]
    scheme_df_dict[scheme] = sub_df
    for scheme_idx in schemes:
        if scheme_idx:
            sub_df = mlst_df.loc[:, [col for col in mlst_df.columns if '.' in col]]
            scheme = sub_df.loc[sub_df.index[0], [_ for _ in sub_df.columns if 'Scheme.' + str(scheme_idx) in _]][0]
            scheme_df_dict[scheme] = sub_df
    scheme_df_dict = {k: redefine_mlst(_df, scheme=k) for k, _df in scheme_df_dict.items()}

    for scheme, _df in scheme_df_dict.items():
        output_mlst[scheme] = _df
        # _df.to_csv(output_mlst.format(scheme=scheme), sep=',', index=True)
    return output_mlst

if __name__ == '__main__':
    import argparse
    from pipelines.tasks import valid_path

    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--infile", help='mlst raw output df')
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
    for scheme,odf  in output_mlst.items():
        with open(os.path.join(odir,"{sn}_{scheme}_mlst.txt".format(sn=prefix,
                                                                    scheme=scheme)),'w') as f1:
            odf.to_csv(f1,index=1)
    # todo: test is it useful for mlst only which without pandoo
