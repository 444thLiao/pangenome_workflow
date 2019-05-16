import os
import random
import string
import sys
from glob import glob
from subprocess import check_call
from collections import defaultdict
import pandas as pd
from Bio import SeqIO, Phylo
from .get_gene_info import get_gff,add_fea2gff

def validate_table(df):
    if df.loc[df.index.drop_duplicates(),:].shape[0] != df.shape[0]:
        raise Exception("sample name contains duplicates")
    both_null = df.loc[(df.iloc[:,0].isna()) & (df.iloc[:,1].isna()),:]
    if both_null.shape[0] != 0:
        raise Exception("Some rows doesn't have R1 and R2. ")

def run_cmd(cmd, dry_run=False, log_file=None, **kwargs):
    outstream = None
    if type(log_file) == str:
        if os.path.isfile(log_file):
            if os.path.getsize(log_file) == 0:
                outstream = open(log_file, 'a')
        if outstream is None:
            outstream = open(log_file, 'w')
    elif log_file is None:
        outstream = sys.stdout
    else:
        outstream = log_file

    print(cmd, file=outstream)
    outstream.flush()
    if not dry_run:
        check_call(cmd,
                   shell=True,
                   executable="/usr/bin/zsh",
                   stdout=outstream,
                   stderr=outstream,
                   **kwargs)
        outstream.flush()


def valid_path(in_pth,
               check_size=False,
               check_dir=False,
               check_glob=False,
               check_odir=False,
               check_ofile=False):
    if type(in_pth) == str:
        in_pths = [in_pth]
    else:
        in_pths = in_pth[::]
    for in_pth in in_pths:
        in_pth = os.path.abspath(os.path.realpath(in_pth))
        if in_pth is None:
            continue
        if check_glob:
            query_list = glob(in_pth)
            if not query_list:
                raise Exception('Error because of input file pattern %s' % in_pth)
        if check_dir:
            if not os.path.isdir(in_pth):
                raise Exception("Error because %s doesn't exist" % in_pth)
        if check_size:
            if os.path.getsize(in_pth) <= 0:
                raise Exception("Error because %s does not contain content." % in_pth)
        if check_odir:
            if not os.path.isdir(in_pth):
                os.makedirs(in_pth, exist_ok=True)
        if check_ofile:
            odir_file = os.path.dirname(in_pth)
            if not os.path.isdir(odir_file):
                os.makedirs(odir_file, exist_ok=True)
    return True


def batch_ln(source_dir, target_dir, suffix=None,
             dry_run=False,
             log_file=None):
    ln_cmd = "ln -s {in_file} -t %s" % target_dir
    if suffix is not None:
        valid_path(os.path.join(source_dir, '*.' + suffix),
                   check_glob=True)
        in_files = glob(os.path.join(source_dir, '*.' + suffix))
    else:
        in_files = glob(os.path.join(source_dir, '*'))
    for in_file in in_files:
        run_cmd(ln_cmd.format(in_file=in_file),
                dry_run=dry_run, log_file=log_file)
        log_file.flush()


def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_letters + string.digits
    return ''.join(random.choice(letters) for _ in range(stringLength))


def get_length_fasta(fasta):
    fh = SeqIO.parse(fasta, format='fasta')
    return {_.id: len(_.seq) for _ in fh}


def get_locus2group(roary_dir):
    protein_cluster_info = os.path.join(roary_dir, 'clustered_proteins')
    locus2group = dict()
    for row in open(protein_cluster_info).readlines():
        row = row.strip('\n')
        cluster_group, remained_locus = row.split(': ')
        for locus in remained_locus.split('\t'):
            locus2group[locus] = cluster_group
    return locus2group


def get_fasta_by_ID(roary_dir, seqid, output_fasta=None):
    # mainly for fetch reference sequence of group
    fa = SeqIO.parse(os.path.join(roary_dir, 'pan_genome_reference.fa'), format='fasta')
    results = []
    for record in fa:
        if record.description.endswith(seqid):
            results.append(record)
    if output_fasta is not None:
        SeqIO.write(results, open(output_fasta, 'w'), format='fasta')
    else:
        return results


def get_tree(tree_pth, rooted=False):
    t = Phylo.read(tree_pth, 'newick')
    if rooted is False:
        return t
    elif rooted == 'midpoint':
        t.root_at_midpoint()
        return t
    elif rooted in [_.name for _ in t.get_terminals()]:
        rooted_node = [_ for _ in t.get_terminals() if _.name == rooted]
        t.root_with_outgroup(rooted_node[0])
        return t
    return t


def construct_pandoo_table(samples_name, odir):
    """ Only for this pipelines, mostly file structure is fixed."""
    isolates_df = pd.DataFrame()
    for sn in samples_name:
        isolates_df = isolates_df.append(pd.DataFrame(data=[[sn,
                                                             os.path.join(odir, "assembly_o", "regular", sn, "contigs.fasta"),
                                                             os.path.join(odir, "cleandata", sn + '_R1.clean.fq.gz'),
                                                             os.path.join(odir, "cleandata", sn + '_R2.clean.fq.gz')
                                                             ]], index=[0]))
    return isolates_df


def group_specific(df, group_dict, threshold=0):
    """
    :param df: 1/0 dataframe instead of np.nan dataframe
    :param group_dict: {groupA:[s1,s2],groupB:[s3,s4]}
    :return:
    """
    # df = df.loc[:,[idx for idx,v in df.iteritems() if pd.api.types.is_numeric_dtype(v)]]
    group_specific_cols = {}
    for g1 in group_dict.keys():
        g2 = set(group_dict.keys()) - {g1}
        g1_vals = group_dict[g1]
        g2_vals = [v for g in g2 for v in group_dict[g]]

        g1_specific = df.columns[(df.loc[g1_vals, :].fillna(0).sum(0) >= len(g1_vals) - threshold) &
                                 (df.loc[g2_vals, :].fillna(0).sum(0) <= threshold)]
        group_specific_cols[g1] = list(g1_specific)
    return group_specific_cols


def get_group_dict(df, col):
    series = df.loc[:, col]
    filtered_series = series[~pd.isna(series)]
    group_dict = {v: [] for v in set(filtered_series)}
    for idx, v in filtered_series.items():
        group_dict[v].append(idx)
    return group_dict


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


def parse_source(source, count, other_info=None):
    if source == "plasmid":
        formatted_ID = 'plasmid%s_piece%s' % (str(other_info),
                                              str(count))
        extra_annotated = {}
    elif source == "phigaro":
        formatted_ID = "phage_{:0>5}".format(count)
        extra_annotated = other_info
    elif source == "isescan":
        formatted_ID = "phage_{:0>5}".format(count)
        extra_annotated = other_info
    elif source == "full":
        formatted_ID = ''
        extra_annotated = {}
    else:
        raise Exception
    return formatted_ID, extra_annotated

def get_new_gff(unify_regions,
                empty_gff_obj,
                source="software"):
    # do not subset the gff file, we just remove all existing features
    # and add regions base annotation into it
    # return a list of data
    new_gff_dict = empty_gff_obj.copy(deep=True)

    count = 1
    for region, other_info in unify_regions:
        contig = region.split(':')[0]
        start,end = map(int,
                        region.split(':')[1].split('-'))
        record = new_gff_dict[contig]

        formatted_ID, extra_annotated = parse_source(source, count, other_info)
        add_fea2gff(record,
                    start,
                    end,
                    ID=formatted_ID if formatted_ID else record.id,
                    type=source.capitalize(),
                    source= source,
                    **extra_annotated
                    )
        count += 1
    return list(new_gff_dict.values())

def alter_gff(gff_obj,
                locus2group,
                locus2annotate,
                source,
              unify_regions=None):
    # cut the full length gff into small region of gff
    # and annotated these region (it could not see the neighbour gene)

    count = 1
    collected_records = []
    if unify_regions is None:
        unify_regions = [("%s:%s-%s" % (contig,0,len(record)),
                          None)
                         for contig,record in gff_obj.items()]
        # create whole chrome regions

    for region, other_info in unify_regions:
        contig = region.split(':')[0]
        start, end = map(int,
                         region.split(':')[1].split('-'))

        start, end = sorted([start, end])  #todo: make sure the order is right

        record = gff_obj[contig]
        subrecord = record[start:end]

        formatted_ID, extra_annotated = parse_source(source, count, other_info)
        subrecord.id = formatted_ID
        subrecord.name = subrecord.description = ''  # remove the description

        for cds in subrecord.features:
            locus_id = cds.id  # locus ID
            if locus_id in locus2annotate.keys():
                group = locus2annotate[locus_id]
            else:
                group = locus2group.get(locus_id, locus_id)

            cds.qualifiers["ID"] = cds.qualifiers['locus_tag'] = [group]
        count += 1
        collected_records.append(subrecord)
    return collected_records


def summary_into_matrix(sample2pieces):
    # from annotated sub_gff, we could summarized a matrix for ML or other statistic analysis
    samples2genes_among_regions = defaultdict(lambda: defaultdict(int) )
    for sample, sub_gffs in sample2pieces.items():
        for fea in [fea for sub_gff in sub_gffs for fea in sub_gff.features]:
            locus_id = fea.id
            samples2genes_among_regions[sample][locus_id] += 1

    return samples2genes_among_regions

