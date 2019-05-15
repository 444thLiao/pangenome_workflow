import os
import time
from glob import glob

import pandas as pd
from Bio import Entrez
from tqdm import tqdm

Entrez.email = 'l0404th@gmail.com'


def fetch(nid, output, start=None, end=None, max_try=50):
    record = ''
    count = 0
    while not record:
        try:
            if not start or not end:
                record = Entrez.efetch(db='nuccore', id=nid, retmode='text', rettype='fasta').read()
            else:
                record = Entrez.efetch(db='nuccore', id=nid, retmode='text', rettype='fasta', seq_start=start,
                                       seq_stop=end).read()
            record = record.split('\n')
            record = record[0] + '\n' + ''.join(record[1:])
        except:
            time.sleep(1)
        count += 1
        if count >= max_try:
            break
    if os.path.isfile(output):
        if '(' in output:
            output_n = output.partition('(')[0]
            output_number = int(output.split('(')[-1].strip(')'))
            output = output_n + '(%s)' % str(output_number + 1)
        else:
            output = output + '(2)'
    with open(output, 'w') as f1:
        f1.write(record)


if __name__ == '__main__':

    for dfn in glob("/home/liaoth/data2/project/shenzhen_Acinetobacter/plasmids/*.txt"):
        data = pd.read_csv(dfn, index_col=0, sep='\t')
        odir = os.path.join(os.path.dirname(dfn),
                            os.path.basename(dfn).split('_plasmids')[0])
        os.makedirs(odir, exist_ok=True)
        for idx, vals in tqdm(data.iterrows(), total=data.shape[0]):
            refseq_id = vals["RefSeq"] if vals["RefSeq"] != '-' else vals["INSDC"]
            if not os.path.isfile(os.path.join(odir, '%s.fa' % refseq_id)):
                fetch(refseq_id, os.path.join(odir, '%s.fa' % refseq_id))
