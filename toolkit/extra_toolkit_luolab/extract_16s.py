from Bio import SeqIO
from glob import glob
from os.path import *
from tqdm import tqdm
from collections import defaultdict
indir = './pipelines_o/prokka_o'
mode = 'prokka'


name2_rrna = defaultdict(list)
for f in tqdm(glob(join(indir,'*','*.gbk'))):
    records = SeqIO.parse(f,format='genbank')
    name = basename(dirname(f))
    for r in records:
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
                for _16s in rrna_16S:
                    _locus_tag = _16s.qualifiers['locus_tag']
                    _record = _16s.extract(r)
                    _record.description = _record.id = _record.name
                    _record.id = _locus_tag[0]
                    name2_rrna[name].append(_record)
                print('multiple 16s in ', name)
    if name2_rrna.get(name) is None:
        print('no 16s in ',name)


ofile = '/home-user/thliao/data/jjtao_20191113/16S.fasta'
records = [_ for k,v in name2_rrna.items() for _ in v]
with open(ofile,'w') as f1:
    SeqIO.write(records,f1,format='fasta-2line')