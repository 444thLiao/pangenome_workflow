from glob import glob
from tqdm import tqdm
from subprocess import check_call
import os
import pandas as pd
from BCBio import GFF
from Bio import SeqIO
from subprocess import check_output
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from toolkit.utils import get_locus2group, run_cmd, valid_path
from toolkit.get_gene_info import get_gff, add_fea2gff
from toolkit.utils import run_cmd, get_locus2group, get_length_fasta, valid_path
from toolkit.get_gene_info import get_gff, add_fea4plasmid, get_gff_pth
from toolkit.process_plasmid import get_gene_in_plasmids

def get_accessory(roary_dir,abricate_file,prokka_o):
    from toolkit.process_IS import get_IS_CDS, get_locus2group
    locus2group = get_locus2group(roary_dir)

    locus2annotate_df = pd.read_csv(abricate_file, sep=',', index_col=0)
    locus2annotate = locus2annotate_df.loc[:, 'gene'].to_dict()

    prepare_dict = {}
    for sn in map(str, locus2annotate_df['sample'].unique()):
        gff_pth = os.path.join(prokka_o, "{sn}", '{sn}.gff').format(sn=sn)
        gff_db = get_gff(gff_pth, mode='db')
        gff_obj = get_gff(gff_pth, mode='bcbio')
        # remove all existing features
        for record in gff_obj.values():
            record.features.clear()
            record.annotations['sequence-region'].clear()
            record.annotations['sequence-region'] = "%s %s %s" % (record.id, 1, len(record))
        prepare_dict[sn] = (gff_db, gff_obj)

    return locus2group,locus2annotate,prepare_dict

def write_gff(phigaro_tab_pth, ori_gff):
    ori_gff = ori_gff.path
    phigaro_tab_pth = phigaro_tab_pth.path
    phigaro_tab = pd.read_csv(phigaro_tab_pth, sep='\t', index_col=0)
    # header: scaffold	begin	end	transposable	taxonomy	vog
    gff_obj = get_gff(ori_gff, mode='bcbio')
    for record in gff_obj.values():
        record.features.clear()
        record.annotations['sequence-region'].clear()
        record.annotations['sequence-region'] = "%s %s %s" % (record.id, 1, len(record))
    count = 0
    for contig, vals in phigaro_tab.iterrows():
        record = gff_obj[contig]
        add_fea2gff(record,
                    vals["begin"],
                    vals["end"],
                    ID="prophage_{:0>5}".format(count),
                    type="Prophage",
                    source="phigaro",
                    phage_tax=vals["taxonomy"],
                    transposable=vals["transposable"],
                    vogID=vals["vog"]
                    )
        count += 1
    with open(phigaro_tab_pth.replace('.out',".gff"),'w') as f1:
        GFF.write(list(gff_obj.values()),f1)



def main(phage_dict, phage_genes, phage_records, full_records,
         indir, prokka_dir, odir):

    only_p = os.path.join(odir, "onlyphage")
    valid_path([only_p], check_odir=1)

    samples_name = phage_dict.keys()

    summary_df = pd.DataFrame(
        columns=['total contigs',
                 'total CDS',
                 'total length',
                 'contigs belong to plasmid',
                 'CDS belong to plasmid',
                 'total length of plasmids',
                 'ratio of plasmid'])
    for sample_name in samples_name:
        contig_pth = os.path.join(indir,
                                  'regular',
                                  sample_name,
                                  'contigs.fa')
        num_contigs = int(check_output("grep -c '^>' %s " % contig_pth, shell=True))
        length_contigs = sum(get_length_fasta(contig_pth).values())

        gff_pth = get_gff_pth(prokka_dir, sample_name)
        # avoid missing formatted prokka dir.

        num_CDS = int(check_output("grep -c 'CDS' %s " % gff_pth, shell=True))

        if os.path.getsize(plasmidcontig_pth) == 0:
            length_plasmidscontigs = num_contigs4plasmid = 0
        else:
            length_plasmidscontigs = sum(get_length_fasta(plasmidcontig_pth).values())

        _sub_df = pd.DataFrame(columns=summary_df.columns,
                               index=[sample_name],
                               data=[[num_contigs,
                                      num_CDS,
                                      length_contigs,
                                      num_contigs4plasmid,
                                      len(phage_genes.get(sample_name, [])),
                                      length_plasmidscontigs,
                                      length_plasmidscontigs / length_contigs
                                      ]])
        summary_df = summary_df.append(_sub_df)
        if phage_records[sample_name]:
            with open(os.path.join(only_p, "%s_phage.fa" % sample_name), 'w') as f1:
                SeqIO.write(phage_records[sample_name], f1, format='fasta')
            with open(os.path.join(only_p, "%s_phage.gff" % sample_name), 'w') as f1:
                GFF.write(phage_records[sample_name], f1)

            with open(os.path.join(fullwithannotated, "%s_phageannotated.gff" % sample_name), 'w') as f1:
                GFF.write(full_records[sample_name], f1)
    summary_pth = os.path.join(odir, "phage_summary.csv")
    summary_df.to_csv(summary_pth, sep=',')



if __name__ == '__main__':
    import argparse

    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--indir", help='')
    parse.add_argument("-o", "--outdir", help='')

    args = parse.parse_args()
    main(os.path.abspath(args.indir),
         os.path.abspath(args.outdir))