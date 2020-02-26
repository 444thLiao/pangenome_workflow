from Bio import SeqIO
from glob import glob
from os.path import *
from tqdm import tqdm
from collections import defaultdict

#ffffff
def extract_16s(indir):
    name2_rrna = defaultdict(list)
    for f in tqdm(glob(join(indir, '*', '*.gbk'))):
        name = basename(dirname(f))
        if not exists(f):
            print('no sequence in ', name)
            continue
        records = SeqIO.parse(f, format='genbank')

        for r in records:
            # each r is a contig
            rrna = [_ for _ in r.features if _.type == 'rRNA']
            if rrna:
                rrna_16S = [_
                            for _ in rrna
                            if "16S ribosomal RNA" in _.qualifiers['product']]
                if len(rrna_16S) == 1 and name2_rrna.get(name) is not None:
                    print('multiple 16s in ', name)
                    _locus_tag = rrna_16S[0].qualifiers['locus_tag']
                    _record = rrna_16S[0].extract(r)
                    _record.description = _record.id = _record.name
                    _record.id = _locus_tag[0]
                    name2_rrna[name].append(_record)
                elif len(rrna_16S) == 1:
                    _locus_tag = rrna_16S[0].qualifiers['locus_tag']
                    _record = rrna_16S[0].extract(r)
                    _record.description = _record.id = _record.name
                    _record.id = _locus_tag[0]
                    name2_rrna[name].append(_record)
                elif len(rrna_16S) > 1:
                    # multiple 16S within one contig
                    for _16s in rrna_16S:
                        _locus_tag = _16s.qualifiers['locus_tag']
                        _record = _16s.extract(r)
                        _record.description = _record.id = _record.name
                        _record.id = _locus_tag[0]
                        name2_rrna[name].append(_record)
                    print('multiple 16s(single contig) in ', name)
        if name2_rrna.get(name) is None:
            print('no 16s in ', name)
    return name2_rrna


if __name__ == '__main__':
    import sys

    assert len(sys.argv) == 3


    # extract_16s.py ./pipelines_o/prokka_o /home-user/thliao/data/jjtao_20191113/16S.fasta
    def process_path(path):
        if not '/' in path:
            path = './' + path
        return abspath(path)


    #     indir = './pipelines_o/prokka_o'
    # mode = 'prokka'
    # ofile = '/home-user/thliao/data/jjtao_20191113/16S.fasta'
    indir = process_path(sys.argv[1])
    ofile = process_path(sys.argv[2])
    name2_rrna = extract_16s(indir)
    records = [_ for k, v in name2_rrna.items() for _ in v]
    with open(ofile, 'w') as f1:
        SeqIO.write(records, f1, format='fasta-2line')
