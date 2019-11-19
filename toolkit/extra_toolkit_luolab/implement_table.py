from subprocess import check_call
from glob import glob
from os.path import *
from tqdm import tqdm
import multiprocessing as mp

columns = ['ID',
           'Quality',
           'Collaborator',
           'assemble time',
           'sample time',
           '16S_genus',
           '16S_species',
           '16S_similarity',
           'additional_16S_genus',
           'additional_16S_species',
           'additional_16S_similarity',
           'Completeness',
           'Contamination',
           'Heterogeneity',
           'Genome size (bp)',
           'contigs',
           'N50',
           'GC',
           'Coding density',
           'Genes',
           'depth']

data_input = './data_input_SZCCT_modified.tab'
Collaborator = "Tao Jinjin"
assemble_time = '2019.11'
sample_time = '-'
###
annotate_16s = './matchResults.csv'

import pandas as pd

_df = pd.read_csv(data_input, sep='\t', index_col=0)
all_ids = list(_df.index)


def run(cmd):
    check_call(cmd, shell=1)


params = []
for fdir in tqdm(glob('./pipelines_o/assembly_o/regular/*/*.fa')):
    fdir = dirname(fdir)
    name = basename(fdir)
    if name in all_ids:
        cmd = "/home-user/thliao/anaconda3/envs/py2/bin/checkm taxonomy_wf domain Bacteria -t 60 -x .fa {fdir} ./tmp/{name}".format(fdir=fdir, name=name)
        if not exists("./tmp/{name}".format(name=name)):
            check_call(cmd, shell=1,
                       )
            params.append(cmd)
# with mp.Pool(processes=2) as tp:
#     list(tp.imap(run,tqdm(params)))
for f in tqdm(glob('./tmp/*/*.ms')):
    fdir = dirname(f)
    name = basename(dirname(f))
    ofile = './tmp/{name}.qa'.format(name=name)
    cmd = "/home-user/thliao/anaconda3/envs/py2/bin/checkm qa {msfile} {fdir} -f {ofile} --tab_table -o 2".format(msfile=f,
                                                                                                                  fdir=fdir,
                                                                                                                  ofile=ofile)
    try:
        if not exists(ofile):
            check_call(cmd,
                       shell=1, )
    except:
        pass
# 16s annoated
annotate_16s = pd.read_csv(annotate_16s, sep='\t')
derep_a = annotate_16s.drop_duplicates('sequence ID')
derep_a.index = [_.split('_')[0] for _ in derep_a.loc[:, 'sequence ID']]

final_df = pd.DataFrame(columns=columns)

final_df.loc[:, 'ID'] = all_ids
final_df.loc[:, 'Collaborator'] = Collaborator
final_df.loc[:, 'assemble time'] = assemble_time
final_df.loc[:, 'sample time'] = sample_time
final_df.index = all_ids
for gid in all_ids:
    if gid not in derep_a.index:
        final_df.loc[gid, 'Quality'] = 'not-qualified'
        continue
    sub_df = derep_a.loc[gid]
    if len(sub_df.shape) == 1:
        t_hit = sub_df['hitTaxonomy']
        genus = t_hit.split(';;')[-2].split(';')[0]
        species = t_hit.split(';;')[-2].split(';')[1]
        final_df.loc[gid, ['16S_genus',
                           '16S_species',
                           '16S_similarity', ]] = [genus, species, sub_df['similarity']]
    else:
        sub_df = sub_df.sort_values('similarity', ascending=False).iloc[:2, :]
        _df = sub_df.iloc[0, :]
        t_hit = _df['hitTaxonomy']
        genus = t_hit.split(';;')[-2].split(';')[0]
        species = t_hit.split(';;')[-2].split(';')[1]
        final_df.loc[gid, ['16S_genus',
                           '16S_species',
                           '16S_similarity', ]] = [genus, species, _df['similarity']]
        _df = sub_df.iloc[1, :]
        t_hit = _df['hitTaxonomy']
        genus = t_hit.split(';;')[-2].split(';')[0]
        species = t_hit.split(';;')[-2].split(';')[1]
        final_df.loc[gid, ['additional_16S_genus',
                           'additional_16S_species',
                           'additional_16S_similarity', ]] = [genus, species, _df['similarity']]
    try:
        sub_df = pd.read_csv('./tmp/{name}.qa'.format(name=gid), sep='\t')
        final_df.loc[gid, ['Completeness',
                           'Contamination',
                           'Heterogeneity',
                           'Genome size (bp)',
                           'contigs',
                           'N50',
                           'GC',
                           'Coding density',
                           'Genes', ]] = list(sub_df.loc[0, ['Completeness',
                                                             'Contamination',
                                                             'Strain heterogeneity',
                                                             'Genome size (bp)',
                                                             '# contigs',
                                                             'N50 (contigs)',
                                                             'GC',
                                                             'Coding density',
                                                             '# predicted genes'
                                                             ]])
    except:
        pass
    tmp = open('pipelines_o/assembly_o/regular/{name}/shovill.log'.format(name=gid))
    row = [row for row in tmp if 'Estimated sequencing depth:' in row]
    depth = row[0].split(' ')[-2].strip('\n')
    depth = '%sx' % depth
    final_df.loc[gid, 'depth'] = depth

for gid, row in final_df.iterrows():
    if (not pd.isna(row['additional_16S_genus'])) and (row['additional_16S_genus'] != row['16S_genus']):
        print(gid)
        row['Quality'] = 'not-qualified'
    elif row['Completeness'] <= 95 or row['Contamination'] >= 5:
        # print(gid)
        row['Quality'] = 'not-qualified'
    else:
        row['Quality'] = 'qualified'
