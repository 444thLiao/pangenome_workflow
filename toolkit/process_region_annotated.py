import copy
import json
import os
import sys
from collections import defaultdict
from subprocess import check_output

import pandas as pd
from BCBio import GFF
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from toolkit.get_gene_info import add_fea2gff
from toolkit.utils import get_length_fasta, valid_path, get_locus2group
from toolkit.get_gene_info import get_gff, get_gff_pth


def get_accessory_obj(roary_dir,
                      abricate_file,
                      prokka_o):
    # For parse unify data into a new annotated gff
    # we need below data(3 differnt type)
    # 1. roary_dir (*clustered info* for minimize the locus ID)
    # 2. abricate_file  (*locus2annotated* for annotated the `extracted` locus)
    # 3. identical prokka_o and corresponding gff (If this prokka is not the same, the locus ID will be different)
    locus2group = get_locus2group(roary_dir)  #

    locus2annotate_df = pd.read_csv(abricate_file, sep=',', index_col=0)
    locus2annotate = locus2annotate_df.loc[:, 'gene'].to_dict()

    sample2gff = {}
    for sn in map(str, locus2annotate_df['sample'].unique()):
        gff_pth = os.path.join(prokka_o,
                               "{sn}",
                               '{sn}.gff').format(sn=sn)
        # load gff twice for two different purpose
        # 1. db mode for query locus
        # 2. bcbio mode for rewrite the gff file
        gff_db = get_gff(gff_pth, mode='db')
        empty_gff_obj = get_gff(gff_pth, mode='bcbio')
        # remove all existing features
        for record in empty_gff_obj.values():
            record.features.clear()
            record.annotations['sequence-region'].clear()
            # record.annotations['sequence-region'] = "%s %s %s" % (record.id,
            #                                                       1,
            #                                                       len(record))
        gff_obj = get_gff(gff_pth, mode='bcbio')
        sample2gff[sn] = (gff_db, empty_gff_obj, gff_obj)

    return locus2group, locus2annotate, sample2gff


def write_new_gff(unify_regions_pth,
                  sample2gff,
                  source="software"):
    # sample2gff must be empty
    # do not subset the gff file, we just remove all existing features
    # and add regions base annotation into it
    # return a list of data
    unify_regions = pd.read_csv(unify_regions_pth, sep='\t', index_col=0, dtype=str)
    new_sample2gff = copy.deepcopy(sample2gff)

    for region_ID, vals in unify_regions.iterrows():
        region = vals["region"]
        sample = str(vals["sample"])
        other_info = vals["other info"]
        contig = region.split(':')[0]
        start, end = map(int,
                         region.split(':')[1].split('-'))

        empty_gff_obj = new_sample2gff[sample]
        record = empty_gff_obj[contig]

        formatted_ID = region_ID
        if type(other_info) == str:
            extra_annotated = json.loads(other_info)
        else:
            extra_annotated = {}

        add_fea2gff(record,
                    start,
                    end,
                    ID=formatted_ID if formatted_ID else record.id,
                    type=source.capitalize(),
                    source=source,
                    **extra_annotated
                    )
    sample2gff_records = {sample: list(gff_dict.values())
                          for sample, gff_dict in new_sample2gff.items()}
    return sample2gff_records


def cut_old_gff(unify_regions_pth,
                sample2gff,
                source,
                locus2annotate=None):
    # sample2gff must be annotated
    # cut the full length gff into small region of gff
    # and annotated these region (it could not see the neighbour gene)

    unify_regions = pd.read_csv(unify_regions_pth, sep='\t', index_col=0, dtype=str)
    annotated_sample2gff = copy.deepcopy(sample2gff)
    sample2gff_records = {}
    for sample in unify_regions["sample"].unique():
        sample = str(sample)
        sub_regions_df = unify_regions.loc[unify_regions["sample"] == sample, :]
        annotated_contig2record = annotated_sample2gff[sample]
        sample2gff_records[sample] = []
        for region_ID, vals in sub_regions_df.iterrows():
            region = vals["region"]
            other_info = vals["other info"]
            contig = region.split(':')[0]
            start, end = map(int,
                             region.split(':')[1].split('-'))
            start, end = sorted([start, end])  # todo: make sure the order is right
            record = annotated_contig2record[contig]
            subrecord = record[start:end]

            if type(other_info) == str:
                extra_annotated = json.loads(other_info)
            else:
                extra_annotated = {}
            # add a full length annotated
            if source != "full":
                add_fea2gff(subrecord,
                            0,
                            len(subrecord),
                            ID=region_ID,
                            type=source.capitalize(),
                            source=source,
                            **extra_annotated
                            )
            if locus2annotate is not None:
                for cds in subrecord.features:
                    locus_id = cds.id  # locus ID
                    annotated = locus2annotate.get(locus_id, locus_id)
                    cds.qualifiers["ID"] = cds.qualifiers['locus_tag'] = [annotated]

            sample2gff_records[sample].append(subrecord)
    return sample2gff_records


def summary_into_matrix(sample2pieces, unique_by=None):
    # from annotated sub_gff, we could summarized a matrix for ML or other statistic analysis
    samples2genes_among_regions = defaultdict(lambda: defaultdict(int))
    for sample, records in sample2pieces.items():
        for fea in [fea for record in records for fea in record.features]:
            if fea.type == 'CDS':
                if unique_by is None:
                    summarized_id = fea.id
                else:
                    summarized_id = fea.qualifiers[unique_by][0]
                samples2genes_among_regions[sample][summarized_id] += 1
    samples2annotated_df = pd.DataFrame.from_dict(samples2genes_among_regions, orient='index')
    samples2annotated_df = samples2annotated_df.fillna(0)
    return samples2annotated_df


def summary_statistic(ori_sample2gff,
                      region_sample2gff,
                      name):
    summary_df = pd.DataFrame(
        columns=['total contigs',
                 'total CDS',
                 'total length',
                 "num of %s" % name,
                 'num of contig have %s' % name,
                 'num CDS covered by %s' % name,
                 'total length of %s' % name
                 ])

    for sample, ori_records in ori_sample2gff.items():
        num_contig = len(ori_records.keys())
        num_CDS = sum([len(record.features) for record in ori_records.values()])
        total_length = sum([len(record) for record in ori_records.values()])

        if sample in region_sample2gff:
            region_records = region_sample2gff[sample]
            num_regions = len(region_records)
            num_contig_r = len(set([region_record.id
                                    for region_record in region_records]))
            num_cds_r = sum([len(region_record.features) - 1  # it have a full region features add by `cut_old_gff`
                             for region_record in region_records])
            total_length_r = sum([len(region_record)
                                  for region_record in region_records])
        else:
            num_regions = num_contig_r = num_cds_r = total_length_r = 0
            
        summary_df = summary_df.append(pd.DataFrame([[num_contig,
                                                      num_CDS,
                                                      total_length,
                                                      num_regions,
                                                      num_contig_r,
                                                      num_cds_r,
                                                      total_length_r]],
                                                    index=[sample],
                                                    columns=summary_df.columns))

    return summary_df


# def write_gff(phigaro_tab_pth, ori_gff):
#     ori_gff = ori_gff.path
#     phigaro_tab_pth = phigaro_tab_pth.path
#     phigaro_tab = pd.read_csv(phigaro_tab_pth, sep='\t', index_col=0)
#     # header: scaffold	begin	end	transposable	taxonomy	vog
#     gff_obj = get_gff(ori_gff, mode='bcbio')
#     for record in gff_obj.values():
#         record.features.clear()
#         record.annotations['sequence-region'].clear()
#         record.annotations['sequence-region'] = "%s %s %s" % (record.id, 1, len(record))
#     count = 0
#     for contig, vals in phigaro_tab.iterrows():
#         record = gff_obj[contig]
#         add_fea2gff(record,
#                     vals["begin"],
#                     vals["end"],
#                     ID="prophage_{:0>5}".format(count),
#                     type="Prophage",
#                     source="phigaro",
#                     phage_tax=vals["taxonomy"],
#                     transposable=vals["transposable"],
#                     vogID=vals["vog"]
#                     )
#         count += 1
#     with open(phigaro_tab_pth.replace('.out', ".gff"), 'w') as f1:
#         GFF.write(list(gff_obj.values()), f1)


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
