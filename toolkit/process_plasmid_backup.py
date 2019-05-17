import os
import re
import sys
from collections import defaultdict
from glob import glob
from subprocess import check_output
import pandas as pd
from BCBio import GFF
from Bio import SeqIO
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from toolkit.utils import run_cmd, get_locus2group, get_length_fasta, valid_path
from toolkit.get_gene_info import get_gff, add_fea4plasmid, get_gff_pth


def align_plasmid(indir):
    # parse the plasmid result into a unify format of data
    # unify format:
    # {sample_name : [(chrom:start-end, other info),
    #                 (chrom2:start2-end2, other info)]
    #                 }
    aligned_info = {}
    for contig in tqdm(glob(os.path.join(indir,
                                         'regular',
                                         '*',
                                         'contigs.fa'))):
        # For each contig, align to corresponding assembly fasta to get matched positions.
        run_cmd("bwa index %s" % contig, dry_run=False, log_file='/tmp/tmp.log')
        plasmid_contig = contig.replace('/regular/',
                                        '/plasmidsSpades/') + 'sta'
        osam = plasmid_contig.replace('.fasta','.sam')
        # regular is 'contigs.fa', plasmids must is 'contigs.fasta'
        if not os.path.isfile(plasmid_contig):
            # it may the SE data, so it should not process plasmid
            continue
        sample_name = os.path.basename(os.path.dirname(contig))
        if os.path.getsize(plasmid_contig) > 0:
            # it may not find plasmid fasta
            if not os.path.isfile(osam):
                sam_content = check_output("bwa mem -x intractg {regular_one} {plasmid_one} > {ofile}".format(regular_one=contig,
                                                                                               plasmid_one=plasmid_contig,
                                                                                                         ofile=osam),
                                      shell=True,
                                      executable='/usr/bin/zsh')
                sam_content = sam_content.decode('utf-8')
            else:
                sam_content = open(osam,'r').read()
            validated_row_of_sam = [_ for _ in sam_content.split('\n') if (not _.startswith('@')) and (_) and (not _.startswith('*'))]

            matched_contig_names = [re.findall("component_([0-9]+)$",
                                    row.split('\t')[0])[0]
                                    for row in validated_row_of_sam]  # parse the contig name of plasmid ot plasmid id
            match_region = ["%s:%s-%s" % (row.split('\t')[2],
                                          int(row.split('\t')[3]) - 1,  # convert 0-coord
                                          int(row.split('\t')[3]) - 1 + len(row.split('\t')[9]))
                            for row in validated_row_of_sam]
            aligned_info[sample_name] = list(zip(match_region,
                                                 matched_contig_names ))

    return aligned_info

def summary_plasmid(indir):
    # summary the unify format data of plasmid result into
    # {sample : {plasmid_ID: [region1,region2]}}
    aligned_info = align_plasmid(indir)
    sample2plasmid_summary = {}
    for sample_name, regions_info in aligned_info.items():
        plasmid_count_dict = defaultdict(list)
        for info in regions_info:
            region,other_info = info
            if not region.startswith('*'):
                # it mean "A QNAME ‘*’ indicates the information is unavailable"
                plasmid_count_dict[other_info].append(region)
        sample2plasmid_summary[sample_name] = plasmid_count_dict
    return sample2plasmid_summary



def get_gene_in_plasmids(plasmids_dict, locus2group, prokka_dir):
    plasmids_genes = {}
    plasmids_records = {}
    full_records = {}
    for sn, p2contig in plasmids_dict.items():
        gff_pth = get_gff_pth(prokka_dir, sn)
        # avoid missing formatted prokka dir.
        gff_dict = get_gff(gff_pth, mode='bcbio')
        # remove all features, make sure the output file only contains detected plasmid
        for record in gff_dict.values():
            record.features.clear()
            record.annotations['sequence-region'].clear()
            record.annotations['sequence-region'] = "%s %s %s" % (record.id, 1, len(record))

        plasmids_genes[sn] = {}
        plasmids_records[sn] = []
        full_records[sn] = []
        for plasmid, p2regions in p2contig.items():
            plasmid_chrome = []
            count = 0
            for region in p2regions:
                contig = region.split(':')[0]
                start, end = map(int, region.split(':')[1].split('-'))
                record = gff_dict[contig]

                subrecord = record[start:end]
                subrecord.id = 'plasmid%s_piece%s' % (str(plasmid),
                                                      str(count))
                subrecord.name = subrecord.description = ''
                add_fea4plasmid(record, start, end, subrecord.id)
                for cds in subrecord.features:
                    cds_id = cds.id
                    group = locus2group.get(cds_id, "removed")
                    if group != 'removed':
                        # todo: process these lost genes
                        plasmids_genes[sn][group] = 1
                        cds.id = group
                        cds.qualifiers["ID"] = cds.qualifiers['locus_tag'] = [group]
                plasmid_chrome.append(subrecord)
                count += 1

            plasmids_records[sn] += plasmid_chrome
        plasmids_records[sn] = list(sorted(plasmids_records[sn],
                                           key=lambda record: tuple(map(int, record.id.strip('plasmid').split('_piece')))))
        full_records[sn] = list(gff_dict.values())
    return plasmids_genes, plasmids_records, full_records


def main(indir, roary_dir, prokka_dir, odir):
    only_p = os.path.join(odir,
                          "onlyplasmid")
    fullwithannotated = os.path.join(odir,
                                     "fullWithAnnotated")
    valid_path([only_p, fullwithannotated], check_odir=1)

    plasmids_dict = get_plasmids(indir)
    locus2group = get_locus2group(roary_dir)
    plasmids_genes, plasmids_records, full_records = get_gene_in_plasmids(plasmids_dict,
                                                                          locus2group,
                                                                          prokka_dir)
    samples_name = plasmids_dict.keys()

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
        plasmidcontig_pth = os.path.join(indir,
                                         'plasmidsSpades',
                                         sample_name,
                                         'contigs.fa')

        if os.path.getsize(plasmidcontig_pth) == 0:
            length_plasmidscontigs = num_contigs4plasmid = 0
        else:
            length_plasmidscontigs = sum(get_length_fasta(plasmidcontig_pth).values())
            num_contigs4plasmid = sum([len(_) for _ in plasmids_dict[sample_name].values()])
        _sub_df = pd.DataFrame(columns=summary_df.columns,
                               index=[sample_name],
                               data=[[num_contigs,
                                      num_CDS,
                                      length_contigs,
                                      num_contigs4plasmid,
                                      len(plasmids_genes.get(sample_name, [])),
                                      length_plasmidscontigs,
                                      length_plasmidscontigs / length_contigs
                                      ]])
        summary_df = summary_df.append(_sub_df)
        if plasmids_records[sample_name]:
            with open(os.path.join(only_p, "%s_plasmid.fa" % sample_name), 'w') as f1:
                SeqIO.write(plasmids_records[sample_name], f1, format='fasta')
            with open(os.path.join(only_p, "%s_plasmid.gff" % sample_name), 'w') as f1:
                GFF.write(plasmids_records[sample_name], f1)

            with open(os.path.join(fullwithannotated, "%s_plasmidannotated.gff" % sample_name), 'w') as f1:
                GFF.write(full_records[sample_name], f1)
    summary_pth = os.path.join(odir, "plasmid_summary.csv")
    summary_df.to_csv(summary_pth, sep=',')


if __name__ == '__main__':
    import argparse

    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--indir", help='')
    parse.add_argument("-r", "--roary_dir", help='')
    parse.add_argument("-p", "--prokka_dir", help='')
    parse.add_argument("-o", "--outdir", help='')

    args = parse.parse_args()
    indir = os.path.abspath(args.indir)
    odir = os.path.abspath(args.outdir)
    roary_dir = os.path.abspath(args.roary_dir)
    prokka_dir = os.path.abspath(args.prokka_dir)

    os.makedirs(odir, exist_ok=True)
    main(indir, roary_dir, prokka_dir, odir)
