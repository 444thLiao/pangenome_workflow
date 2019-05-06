import os
import re
from collections import defaultdict
from glob import glob
from subprocess import check_call

import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from toolkit.utils import get_locus2group

samples2locus = defaultdict(list)
for ffn_pth in tqdm(glob('/home/liaoth/project/shenzhen_actineto/prokka_o/*/*.ffn')):  # server
    for db in ['card', 'argannot', 'ncbi', 'vfdb', 'vfdb_full', 'resfinder', 'plasmidfinder', 'victors']:
        cmd = "abricate {infile} --threads 40 --db {db} --mincov 80 --quiet> {ofile}"
        ofile = os.path.join(os.path.dirname(ffn_pth), 'abricate_{db}.tab'.format(db=db))
        new_cmd = cmd.format(infile=ffn_pth,
                             db=db,
                             ofile=ofile)
        # check_call(new_cmd, shell=True)
        ############################################################
    samples_name = os.path.basename(os.path.dirname(ffn_pth))
    samples2locus[samples_name] = [_.id for _ in SeqIO.parse(ffn_pth,format='fasta')]

def gene_process(gene):
    gene = gene.rsplit('.', 1)[0]
    if 'full' in gene:
        gene = gene.split('full_')[-1]
    if gene.endswith('_1') or gene.endswith('_2'):
        gene = re.split('_1', gene)[0]
        gene = re.split('_2', gene)[0]
    if gene.startswith('('):
        gene = gene.split(')')[1]
    if gene.startswith('bla'):
        gene = gene.split('bla')[1]
    return gene

locus2annotate = {}
annotate2db = {}
for tab in tqdm(glob("/home/liaoth/project/shenzhen_actineto/prokka_o/*/abricate_*.tab")):
    df = pd.read_csv(tab, sep='\t')
    seq_id = list(df.loc[:, "SEQUENCE"])
    genes = [gene_process(_) for _ in df.loc[:, 'GENE']]
    db = os.path.basename(tab).split('_')[1].split('.')[0]
    annotate2db.update({g:db for g in genes})
    # only compare among one sample.
    prepare_update_dict = dict(zip(seq_id, genes))
    for k in seq_id:
        if k in locus2annotate.keys():
            ori = locus2annotate[k]
            aft = prepare_update_dict[k]
            if len(ori) < len(aft):
                prepare_update_dict.pop(k)
    locus2annotate.update(prepare_update_dict)

roary_dir = '/home/liaoth/project/shenzhen_actineto/65_roary_o'
locus2group = get_locus2group(roary_dir)
group2annotate = defaultdict(list)
for locus, annotate in locus2annotate.items():
    formatted_annotate = gene_process(annotate)
    if locus not in locus2group:
        # it mean the roary dir is not full
        print(locus)
        pass
    else:
        group2annotate[locus2group[locus]].append(formatted_annotate)

# for unify similary annotate genes
rename_genes = defaultdict(list)
# check for, exam one group match too many annotated genes
for g, a in group2annotate.items():
    if len(set(a)) > 1:
        unify_name = min(set(a),key=lambda x:len(x))
        for _ in set(a):
            rename_genes[_] = unify_name

samples2genes = defaultdict(lambda: defaultdict(lambda:0))

for sample in samples2locus.keys():
    for locus in samples2locus[sample]:
        if locus in locus2annotate:
            annotate = locus2annotate[locus]
            formatted_annotate = rename_genes.get(annotate,annotate)
            samples2genes[sample][formatted_annotate] += 1

for locus,annotate in locus2annotate.items():
    locus2annotate[locus] = rename_genes.get(annotate,annotate)
locus2annotate_df = pd.DataFrame.from_dict({0:locus2annotate})
locus2annotate_df.loc[:,'db'] = [annotate2db[_] for _ in locus2annotate_df.loc[:,0]]
locus2annotate_df.to_csv('/home/liaoth/project/shenzhen_actineto/all_65_annotated_locus2annotate.tab',sep='\t',index=1)

abricate_result = pd.DataFrame.from_dict(samples2genes,orient='index')
abricate_result.loc['db',:] = [annotate2db[_] for _ in abricate_result.columns]
abricate_result.to_csv('/home/liaoth/project/shenzhen_actineto/abricate_result.tab',sep='\t')
