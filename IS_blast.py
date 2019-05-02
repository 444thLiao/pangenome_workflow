import os
from subprocess import check_output
from glob import glob
import pandas as pd
import io
cmd = "blastn -query {fasta} -db /home/liaoth/data2/project/shenzhen_actineto/ISfinder/ISfinder-sequences/IS.fna -outfmt 6"



def process_df(df):

    return df
file2df = {}
for fa in glob("/home/liaoth/data2/project/shenzhen_actineto/ready_prokka/*.fasta"):
    n_cmd = cmd.format(fasta=fa)
    output = check_output(n_cmd,shell=True)
    output = output.decode('utf-8')
    #header = []
    data = pd.read_csv(io.StringIO(output),sep='\t',header=None)
    data = process_df(data)
    base_name = os.path.basename(fa).split('.')[0]
    file2df[base_name] = data