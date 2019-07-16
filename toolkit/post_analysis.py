import sys
from os.path import dirname, realpath, join

sys.path.insert(0, dirname(dirname(__file__)))
from BCBio import GFF
from toolkit.batch_abricate import get_abricate_df
from toolkit.process_region_annotated import *
from toolkit.utils import valid_path
from pipelines import constant_str as constant
import click

@click.command()
@click.argument("-r","--roary_d","roary_dir")
@click.argument("-o","--output_d","output_dir")
@click.argument("-p","--prokka_d","prokka_o")
@click.argument("-a","--abricate_f","abricate_file")
def cli(roary_dir, output_dir, prokka_o, abricate_file):
    post_analysis(roary_dir,
                  output_dir,
                  prokka_o,
                  abricate_file)


def post_analysis(roary_dir, output_dir, prokka_o=None, abricate_file=None):
    roary_dir = realpath(roary_dir)
    root_output_dir = dirname(roary_dir)
    summary_dir = join(root_output_dir, constant.summary_dir)
    ######################
    output_dir = realpath(output_dir)

    if prokka_o is None:
        prokka_o = join(root_output_dir, 'prokka_o')
    abricate_odir = join(root_output_dir, 'abricate_result')
    abricate_result = pd.DataFrame()
    if abricate_file is None:
        print("summarizing the abricate output")
        abricate_file, abricate_result = get_abricate_df(prokka_dir=prokka_o,
                                                         roary_dir=roary_dir,
                                                         abricate_odir=abricate_odir)
        # abricate_file = workflow_task.input()["abricate"].path
    valid_path(output_dir, check_odir=1)
    ############################################################
    # For summary the region based annotated result
    # 1. write info to empty gff for each sample
    # 2. extract each region with annotated info for each sample
    # 3. summary into matrix
    # 4. summary statistic info
    # prepare the accessory obj
    locus2group, locus2annotate, sample2gff = get_accessory_obj(roary_dir,
                                                                abricate_file,
                                                                prokka_o)
    # get locus annotation(abricate)/group(roary) from different files.
    # sample2gff contains three objs:
    #   gff_db(for query), empty_gff_obj(gff obj but removed features), gff_obj(complete gff obj)
    empty_sample2gff = {sn: vals[1]
                        for sn, vals in sample2gff.items()}
    ori_sample2gff = {sn: vals[2]
                      for sn, vals in sample2gff.items()}

    merged_locus2annotate = locus2group.copy()
    merged_locus2annotate.update(locus2annotate)
    # use annotate to overlap group info.
    # merged them
    ############
    summary_task_source = ["phigaro",
                           "plasmidSpades+bwa",
                           "isescan"]
    names = ["Prophage", "Plasmid", "IS"]
    ofiles = [join(summary_dir,
                   'phage_summary.tab'),
              join(summary_dir,
                   "plasmid_summary.tab"),
              join(summary_dir,
                   'IS_summary.tab')
              ]
    annotated_sample2gff = copy.deepcopy(ori_sample2gff)
    for record in [record
                   for contig2record in annotated_sample2gff.values()
                   for record in contig2record.values()]:
        for fea in record.features:
            if fea.type == 'CDS':
                locus_id = fea.id
                annotated = merged_locus2annotate.get(locus_id, locus_id)
                fea.qualifiers["ID"] = fea.qualifiers['locus_tag'] = [annotated]
                fea.id = annotated
    ###########
    # main part for generate multiple GFF files
    print("generating the GFF files...")
    for ofile, source, name in zip(ofiles,
                                   summary_task_source,
                                   names):
        print(source,name,ofile)
        pth = ofile

        annotated_sample2gff_records = write_new_gff(pth,
                                                     empty_sample2gff,
                                                     source=source)
        subset_sample2gff_records = cut_old_gff(pth,
                                                annotated_sample2gff,
                                                source=source)
        samples2annotated_df = summary_into_matrix(subset_sample2gff_records,
                                                   unique_by=None)
        summary_df = summary_statistic(ori_sample2gff,
                                       subset_sample2gff_records,
                                       name
                                       )
        # output
        full_gff_with_region_dir = os.path.join(output_dir,
                                                "gff_with_%s" % name)
        subset_gff_with_annotated_dir = os.path.join(output_dir,
                                                     "gff_of_%s" % name)
        annotated_gff_odir = os.path.join(output_dir,
                                          "annotated_gff")
        valid_path([full_gff_with_region_dir,
                    subset_gff_with_annotated_dir,
                    annotated_gff_odir
                    ], check_odir=1)
        for sn, records in annotated_sample2gff_records.items():
            filename = "%s.gff" % sn
            with open(os.path.join(full_gff_with_region_dir,
                                   filename), 'w') as f1:
                GFF.write(records, f1)
            with open(os.path.join(subset_gff_with_annotated_dir,
                                   filename), 'w') as f1:
                GFF.write(subset_sample2gff_records[sn], f1)
            with open(os.path.join(annotated_gff_odir,
                                   filename), 'w') as f1:
                GFF.write(list(annotated_sample2gff[sn].values()), f1)
        with open(os.path.join(output_dir, "%s_annotated_matrix.csv" % name), 'w') as f1:
            samples2annotated_df.to_csv(f1, sep=',', index=1)
        with open(os.path.join(output_dir, "%s_statistic.csv" % name), 'w') as f1:
            summary_df.to_csv(f1, sep=',', index=1)

    ############################################################
    # mlst
    mlst_ofile = join(summary_dir,
                      "mlst_all.csv")
    pandoo_df = pd.read_csv(mlst_ofile, index_col=0)
    type_ST = list(pandoo_df.columns)
    type_ST = set([_.split('.')[0] for _ in type_ST])

    for scheme in type_ST:
        with open(os.path.join(output_dir,
                               "%s_mlst.csv" % scheme.replace('_ST', '')), 'w') as f1:
            pandoo_df.loc[:, pandoo_df.columns.str.startswith(scheme)].to_csv(f1, index=1)
    ############################################################
    # abricate
    if type(abricate_file) == str:
        abricate_ofile = abricate_file
        abricate_dir = os.path.dirname(abricate_ofile)
        os.system("cp %s %s" % (abricate_ofile,
                                output_dir))
        os.system("cp %s %s" % (os.path.join(abricate_dir,
                                             "samples2annotate.csv"),
                                output_dir))
    else:
        abricate_ofile = join(output_dir, 'locus2annotate.csv')
        abricate_ofile2 = join(output_dir, "samples2annotate.csv")
        abricate_file.to_csv(abricate_ofile, index=1)
        abricate_result.to_csv(abricate_ofile2, index=1)

    abricate_gff_dir = os.path.join(output_dir,
                                    "annotated_gff_simplified")
    valid_path([abricate_gff_dir], check_odir=1)
    os.system("cp %s %s" % (os.path.join(abricate_odir,
                                         "*",
                                         "*.gff"),
                            abricate_gff_dir))

if __name__ == '__main__':
    cli()