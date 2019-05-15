import os

import gffutils
from BCBio import GFF

tmp_dir = '/tmp'


def get_gff(gff_fn, mode='db'):
    if mode == 'db':
        dbfn = os.path.join(tmp_dir,
                            os.path.basename(gff_fn).rsplit('.', 1)[0] + '.gffdb')
        if os.path.isfile(dbfn):
            os.remove(dbfn)
            fn = gffutils.create_db(gff_fn, dbfn=dbfn, merge_strategy='merge')
        else:
            fn = gffutils.create_db(gff_fn, dbfn=dbfn, merge_strategy='merge')
        return fn
    elif mode == 'bcbio':
        gff_obj = GFF.parse(gff_fn)
        record_dict = {_.id: _ for _ in gff_obj}
        return record_dict


def get_gene_with_regin(gff_fn, regions):
    gff_f = get_gff(gff_fn)
    all_genes = []
    for region in regions:
        for cds in gff_f.region(region=(region), completely_within=True):
            if "ID" in cds.attributes.keys():
                all_genes.append(cds["ID"][0])
            else:
                all_genes.append(cds["note"][0])

    os.system("rm %s/*.gffdb" % tmp_dir)
    return all_genes
