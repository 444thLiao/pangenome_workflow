import itertools
import os
from glob import glob

import pandas as pd
import plotly
import plotly.graph_objs as go
from BCBio import GFF
from Bio import AlignIO
import plotly.io as pio

group2 = ['39292', '39502', '41296', '34960', '35082', '39054', '40067',
          '42349', '43458', '37706', '37502', '35989', '39528']

roary_dir = "/home/liaoth/data2/project/shenzhen_Acinetobacter/group2_roary_o_test"

def hamming(s1, s2):
    """Return the hamming distance between 2 DNA sequences"""
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2)) + abs(len(s1) - len(s2))

gene_diversity = {}
for aln in glob(os.path.join(roary_dir,"pan_genome_sequences","*.aln")):
    gene_name = os.path.basename(aln).split('.fa')[0]
    aln_obj = AlignIO.read(aln,format='fasta')
    length_aln = len(aln_obj)
    num_seq = len(aln_obj[:,0])
    xi = 1/num_seq
    xj = 1/num_seq
    count = 0
    for i,j in itertools.combinations(range(num_seq),2):
        seqi = aln_obj[i,:]
        seqj = aln_obj[j,:]
        if str(seqi.seq) !=str(seqj.seq):
            pi_ij = hamming(seqi,seqj)/length_aln
        else:
            pi_ij = 0
        count += xi*xj*pi_ij
    gene_diversity[gene_name] = count

locus2annotate_df = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/all_65_annotated_locus2annotate.tab",sep='\t',index_col=0)
locus2annotate = dict(zip(locus2annotate_df.index,locus2annotate_df.loc[:,'gene']))

odir = "/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/detect_sweeps"
from toolkit.utils import  get_locus2group
locus2group = get_locus2group(roary_dir)
group2locus = {g: l for l, g in locus2group.items()}  # it will remove redundancy, it does not matter
df = pd.read_csv('/home/liaoth/data2/project/shenzhen_Acinetobacter/group2_roary_o_test/gene_presence_absence.csv',index_col=0)
df = df.iloc[:,13:]
df = df.T
core_genes = df.columns[(~df.isna()).all(0)]
############################################################
fig = go.Figure()
fig.add_scatter(x=[_ for _ in core_genes if _ in gene_diversity],
                y=[gene_diversity[_] for _ in core_genes if _ in gene_diversity],
                mode='markers')
fig.layout.width = 1000
fig.layout.height = 1000
fig.layout.font.size= 15
fig.layout.yaxis.title = "Intra-group gene diversity"
pio.write_image(fig,os.path.join(odir,'overview.png'))
plotly.offline.plot(fig,filename=os.path.join(odir,'overview.html'),auto_open=False)
############################################################
# high diversity genes
threshold = 5
high_diversity_genes = [_ for _ in core_genes if _ in gene_diversity and gene_diversity[_] > threshold]
for g in high_diversity_genes:
    locus = group2locus[g]
    if locus in locus2annotate:
        print(g,locus2annotate[locus])

############################################################
def get_neighbour_genes(centre_group,num_neighbor=20):
    order_subset_locus = {}
    for sn in group2:
        gff_file = '/home/liaoth/data2/project/shenzhen_Acinetobacter/prokka_html/prokka_gff/{sn}.gff'.format(sn=sn)
        gff_obj = GFF.parse(gff_file)
        locus_name = df.loc[sn,centre_group]
        try:
            contig = [_ for _ in gff_obj if locus_name in [record.id for record in _.features]][0]
        except:
            import pdb;pdb.set_trace()
        order_locus_id = [_.id for _ in contig.features if locus2group.get(_.id,'') in core_genes]
        idx = order_locus_id.index(locus_name)
        if idx >= num_neighbor:
            start = idx-num_neighbor
        else:
            start = 0
        if idx+num_neighbor >= len(order_locus_id):
            end = len(order_locus_id)
        else:
            end = idx+num_neighbor
        # print(idx,start, end,len(order_locus_id))
        order_subset_locus[sn] = [locus2group.get(_) for _ in order_locus_id[start:end]]
    return order_subset_locus

#### group specific group_92
# order_subset_locus = get_neighbour_genes('group_92')
# for k,v in order_subset_locus.items():
#     # if k == '41296' or k == '39054':
#     print(k,'\t'.join([_ if _ is not None else 'None' for _ in v[::-1]  ]))
#     else:
#         print(k, '\t'.join([_ if _ is not None else 'None' for _ in v]))

# for k,v in order_subset_locus.items():
#     if k == '41296' or k == '39054':
#         print(k,'\t'.join([str(gene_diversity[_]) if _ is not None else 'None' for _ in v[::-1]  ]))
#     else:
#         print(k, '\t'.join([str(gene_diversity[_]) if _ is not None else 'None' for _ in v]))
target_gene = 'group_92'  # bap

order_subset_locus = get_neighbour_genes(target_gene,20)


fig = go.Figure()
for sample in order_subset_locus.keys():
    xs = [_ for _ in order_subset_locus[sample] if _ is not None]
    ys = [gene_diversity[_] for _ in xs if _ is not None]
    fig.add_scatter(x=xs,
                    y=ys,
                    mode='markers',
                    name=sample,
                    showlegend=False,
                    marker=dict(size=40,
                                color='#1DACE2'))
fig.add_scatter(x=[target_gene],
                y=[gene_diversity[target_gene]],
                showlegend=False,
                mode='markers',
                marker=dict(size=55,
                            color='red')
                )
plotly.offline.plot(fig)
